#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use DBI;
use Data::Dumper;

our %opt;	# Getopt options

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script builds the list of genes in all breaks in a separate table.
	
	-h        : this (help) message
	
	Input:
	-d <path> : path to the sqlite3 database
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hd:', \%opt) or usage();
	usage() if $opt{h};
	usage("Sqlite3 database needed (-d)") unless $opt{d};
}

sub get_genomes
{
	my ($dbh) = @_;
	
	my $query = 'SELECT DISTINCT sp FROM orthos O, genes P WHERE O.pid1=P.pid';
	my $genomes = $dbh->selectcol_arrayref($query) or die("'$query': $!");
	
	return $genomes;
}

sub get_all_breaks
{
	my ($dbh) = @_;

	my $select = 'SELECT * FROM breaks_all';
	my $query = $select;
	
	my $breaks_list = $dbh->selectall_arrayref($query, { Slice => {} });

	# Put in a hash of species/genome_parts
	my %breaks = ();
	foreach my $b (@$breaks_list) {
		push @{ $breaks{ $b->{ 'sp1' } }{ $b->{ 'sp2' } }}, $b;
	}
	
	return \%breaks;
}

sub copy {
	my ($ref) = @_;
	
	my $clone = {};
	if (ref($ref) eq 'HASH') {
		foreach my $k (keys %$ref) {
			$clone->{ $k } = $ref->{$k};
		}
	}
	return $clone;
}

sub get_all_genes
{
	my ($dbh, $genomes) = @_;
	
	my $genes = {};
	print STDERR "Get genes";
	foreach my $sp (@$genomes) {
		print STDERR "\t$sp";
		my $select = 'SELECT pid, pnum_all, loc_start, loc_end, loc_length, product, feat, paralogs, paralogs_n, sp, sp1, sp2, strand, oid, pid2 AS ortho FROM genes G LEFT JOIN orthos_all O ON (G.pid=O.pid1 AND G.sp=O.sp1) WHERE sp=?';
		my $query = $select;
		
		my $genes_list = $dbh->selectall_arrayref($query, { Slice => {} }, $sp);

		# Put in a hash of species
		my ($n, $s) = (0,0);
		my $g =@$genes_list ;
		
		# Be careful: use ALL only if sp2->[n] is not defined
		foreach my $g (@$genes_list) {
			# Keep all genes with orthos with sp2 in a separate file
			if ( $g->{ 'oid' } and $g->{ 'sp2' } ) {
				$genes->{ $sp }->{ $g->{ 'sp2' } }->[ $g->{ 'pnum_all' } ] = $g;
				$s++;
			}
			my $gall = copy($g);
			delete $gall->{ 'oid' };
			delete $gall->{ 'sp2' };
			delete $gall->{ 'ortho' };
			# Still, keep all genes for sp1 (but get rid of everything from any sp2)
			$genes->{ $sp }->{ 'ALL' }->[ $g->{ 'pnum_all' } ] = $gall;
			$n++;
		}
	}
	print STDERR "\n";
	return $genes;
}

sub get_genes_list
{
	my ($breaks, $genes, $genome1, $genome2) = @_;
	my $breaks_genes = {};

	foreach my $b ( @$breaks ) {
		# Take all genes in genome 1
		for(my $i = $b->{ 'pnum_all_left1' }+1 ; $i < $b->{ 'pnum_all_right1' } ; $i++) {
			my $gene = $genes->{ $genome1 }->{ $genome2 }->[ $i ];
			if (not $gene) {
				$gene = $genes->{ $genome1 }->{ 'ALL' }->[ $i ];
			}
			$gene->{ 'sp1' } = $genome1;
			$gene->{ 'sp2' } = $genome2;
			$gene->{ 'side' } = 1;
			push @{ $breaks_genes->{ $b->{ 'breakid' } }->{ 1 } }, $gene;
			#print STDERR "1\t$b->{'breakid'}\t$gene->{'pid'}\t$gene->{'pnum_all'}\n";
		}
		# Take all genes in strand 2
		my $from = $b->{ 'pnum_all_left2' };
		my $to = $b->{ 'pnum_all_right2' };
		# Reversed order
		if ($from > $to) {
			my $tmp = $from;
			$from = $to;
			$to = $tmp;
		}
		for(my $i = $from + 1 ; $i < $to ; $i++) {
			my $gene = $genes->{ $genome2 }->{ $genome1 }->[ $i ];
			if (not $gene) {
				$gene = $genes->{ $genome2 }->{ 'ALL' }->[ $i ];
			}
			$gene->{ 'sp1' } = $genome1;
			$gene->{ 'sp2' } = $genome2;
			$gene->{ 'side' } = 2;
			push @{ $breaks_genes->{ $b->{ 'breakid' } }->{ 2 } }, $gene;
			#print STDERR "2\t$b->{'breakid'}\t$gene->{'pid'}\t$gene->{'pnum_all'}\n";
		}
	}
	
	return $breaks_genes;
}

sub add_genes_list
{
	my ($dbh, $genes_list) = @_;
	
	my @valorder = qw(sp1 sp2 breakid side pid oid ortho ortho_in paralogs paralogs_n pnum_all loc_start loc_end loc_length product feat strand);
	my $valorder_names = join(", ", @valorder);
	my $placeholders = join(",", ('?')x($#valorder+1));
	$dbh->do( 'BEGIN' );
	my $query = "INSERT INTO breaks_genes($valorder_names) VALUES ($placeholders)";
	my $sth = $dbh->prepare( $query );
	foreach my $breakid ( keys %$genes_list ) {
		foreach my $side (1, 2) {
			foreach my $gene (@{ $genes_list->{ $breakid }->{ $side } }) {
				$gene->{'breakid'} = $breakid;
				my @vals = ();
				foreach my $name (@valorder) {
					push @vals, $gene->{ $name };
				}
				$sth->execute( @vals ) or die( Dumper( $gene ) . " in $breakid");
			}
		}
	}
	$dbh->do( 'COMMIT' );
}

sub check_orthos
{
	my ($genes_list) = @_;

	foreach my $breakid ( keys %$genes_list ) {
		my %list = ();
		# List all pid / orthos couples
		foreach my $side (1, 2) {
			foreach my $gene (@{ $genes_list->{ $breakid }->{ $side } }) {
				$list{ $side }{ $gene->{ 'pid' } } = $gene->{ 'ortho' };
			}
		}
		# See if the orthos are in the other side
		foreach my $side (1, 2) {
			my $other_side = ($side == 1) ? 2 : 1;
			foreach my $gene (@{ $genes_list->{ $breakid }->{ $side } }) {
				my $pid = $gene->{ 'pid' };
				my $ortho = $list{ $side }{ $pid };

				# Found in the same break on the other side?
				if ( $ortho and $list{ $other_side }{ $ortho } ) {
					$gene->{ 'ortho_in' } = 1;
				}
				
			}
		}
	}
	return $genes_list;

}

##################################
# MAIN
init();

# Connect
my $dbh = DBI->connect( "dbi:SQLite:dbname=$opt{d}" );

# Retrieve the breaks
my $breaks_list = get_all_breaks($dbh);
my @genomes = sort keys %$breaks_list;

# Retrieve all the genes/orthos for each genome
my $genes = get_all_genes($dbh, \@genomes);

# Create new table
$dbh->do("DROP TABLE IF EXISTS breaks_genes");
my $create_table = <<REQ;
	CREATE TABLE breaks_genes(
		brgid INTEGER PRIMARY KEY,
		sp1 TEXT,
		sp2 TEXT,
		breakid INT REFERENCES breaks(breakid) ON DELETE CASCADE,
		pid INT REFERENCES genes(pid),
		ortho INT REFERENCES genes(pid),
		ortho_in INT DEFAULT 0,
		paralogs_n INT DEFAULT 0,
		paralogs TEXT,
		side INT,
		oid TEXT REFERENCES orthos(oid),
		pnum_all INT,
		loc_start INT,
		loc_end INT,
		loc_length INT,
		product TEXT,
		feat TEXT,
		strand INT,
		UNIQUE (breakid, pid)
	);
REQ
$dbh->do( $create_table );
$dbh->do( "CREATE INDEX idx_brg_breakside ON breaks_genes(breakid, side)" );

print STDERR "Process each break genes list... \n";
foreach my $genome1 (@genomes) {
	my $noo = @{ $genes->{ $genome1 }->{ 'ALL'} };
	print STDERR "$genome1\t$noo\n";
	foreach my $genome2 (@genomes) {
		next if ( $genome1 eq $genome2 );
		# Now retrieve each break genes list
		my $breaks = $breaks_list->{ $genome1 }->{ $genome2 };
		my $genes_list = get_genes_list($breaks, $genes, $genome1, $genome2);
		# Check orthologs inside
		$genes_list = check_orthos($genes_list);
		# Add them to the new table
		add_genes_list( $dbh, $genes_list );
	}
}

__END__

