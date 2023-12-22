#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Data::Dumper;
use Getopt::Std;
use DBI;
use List::Util qw(min max);
use Statistics::Basic qw(median);

our %opt;	# Getopt options

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script orders parts (e.g. contigs) of a genome with another.
	
	-h        : this (help) message
	
	Input:
	-d <path> : path to the sqlite3 database to create
	
	Automatic
	-a        : automatic orders, based on the closest genomes
	or manual
	-m <str>  : name of the model genome
	-s <str>  : name of the genome to align
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hd:am:s:', \%opt) or usage();
	usage() if $opt{h};
	usage("Sqlite3 database needed (-d)") unless $opt{d};
	usage("Mode needed: auto or manual (-a|-m -s)") unless $opt{a} xor ($opt{m} and $opt{s});
	usage("Model genome needed (-m)") if not $opt{a} and not $opt{m};
	usage("Name of the genome to align needed (-s)") if not $opt{a} and not $opt{s};
}

sub get_genes
{
	my ($dbh, $sp) = @_;
	
	my $query = 'SELECT * FROM genes WHERE sp=?';
	my $genes = $dbh->selectall_arrayref($query, { Slice => {} }, ($sp)) or die("'$query': $!");
	#print STDERR "Genes in $sp:\t" . (@$genes+0) . "\n";
	
	return $genes;
}

sub get_orthos
{
	my ($dbh, $sp1, $sp2) = @_;
	
	my $query = 'SELECT * FROM orthos_all WHERE sp1=? AND sp2=?';
	my $orthos = $dbh->selectall_hashref($query, 'pid1', undef, ($sp1, $sp2)) or die("'$query': $!");
	#print STDERR "Orthologs from $sp1 to $sp2:\t" . (keys %$orthos) . "\n";
	
	return $orthos;
}

sub partition
{
	my ($draft) = @_;

	my $parts = {};

	foreach my $gene (@$draft) {
		push @{ $parts->{ $gene->{'gpart'} } }, $gene;
	}
	#print STDERR "Parts:\t" . (keys %$parts) . "\n";
	
	return $parts;
}

sub get_parts_data
{
	my ($draft_parts, $model) = @_;
	my %medians = ();
	my %directions = ();
	my $n = 0;
	foreach my $part_name (sort keys %$draft_parts) {
		my $part = $draft_parts->{ $part_name };
		my @pnums = ();
		my $cumul = 0;
		my $last_pnum2 = undef;
		
		# Read each gene
		foreach my $gene (@$part) {
			my $pid = $gene->{ 'pid' };

			# Does this gene has an ortholog?
			if ( $model->{ $pid } ) {
				my $pnum2 = $model->{ $pid }->{ 'pnum_all2' };
				# Add to the list
				push @pnums, $pnum2;
				$n++;
				# Also, scan to get the direction
				if (defined($last_pnum2)) {
					my $direction = $pnum2 - $last_pnum2;
					$direction = $direction / abs($direction);	# 1 or -1
					$cumul += $direction;
					$last_pnum2 = $pnum2;
				} else {
					$last_pnum2 = $pnum2;
				}
			}
		}
		if (@pnums > 2) {
			if (
				(@pnums < 50 and abs(max(@pnums) - min(@pnums)) < 200)
				or
				abs($cumul) > 20) {
				$medians{ $part_name } = median(@pnums);
			} else {
				#print STDERR "NOT $part_name\n";
				$medians{ $part_name } = 99999999;
			}
			$directions{ $part_name } = $cumul;
			#	print STDERR "\t$part_name\t$medians{ $part_name }\n";
		} elsif (@pnums > 0) {
			$medians{ $part_name } = median(@pnums);
			$directions{ $part_name } = $cumul > 0 ? $cumul : 0;
		} else {
			$medians{ $part_name } = 99999999;	# max value to put it at the end
			$directions{ $part_name } = 0;
			#print STDERR "No ortholog for part $part_name\n";
		}
	}
	if ($n == 0) {
		#print STDERR "Really nothing\n";
		#print STDERR join("\n\t", keys %$model);
		#print STDERR "\n";
	}
	
	my %stats = ('medians' => \%medians, 'directions' => \%directions);
	return \%stats;
}

sub order_parts
{
	my ($dbh, $draft_parts, $stats) = @_;
	
	my @order = sort { $stats->{'medians'}->{$a} <=> $stats->{'medians'}->{$b} or $a cmp $b } keys %$draft_parts;
	#@order = sort keys %$draft_parts;
	
	my $pnum = 0;
	$dbh->do('BEGIN');
	foreach my $part_name (@order) {
		#print STDERR "$part_name\t$stats->{'directions'}->{$part_name}\t$stats->{'medians'}->{ $part_name }\n";
		my $sth = $dbh->prepare_cached('UPDATE genes SET pnum_display=? WHERE pid=?');

		my @gene_order = @{ $draft_parts->{ $part_name } };
		@gene_order = reverse(@gene_order) if $stats->{'directions'}->{$part_name} < 0;
		foreach my $gene (@gene_order) {
			$pnum++;
			$sth->execute($pnum, $gene->{'pid'});
		}
	}
	$dbh->do('COMMIT');
}

sub order_gparts
{
	my ($dbh, $sp) = @_;
	$dbh->do('DELETE FROM genome_parts WHERE sp=?', undef, ($sp));
	my $populate_genome_parts_table = <<"REQ";
	INSERT INTO genome_parts (sp, gpart, min, max)
	SELECT sp, gpart, min(pnum_display), max(pnum_display) FROM genes WHERE sp=? GROUP BY gpart ORDER BY sp, min(pnum_display)
REQ
	$dbh->do($populate_genome_parts_table, undef, ($sp));
}

sub regenerate_blocks_all
{
	my ($dbh) = @_;
	$dbh->do("DROP TABLE IF EXISTS blocks_all");
	$dbh->do("CREATE TABLE blocks_all AS SELECT * FROM blocks_all_view");
	$dbh->do("CREATE INDEX idx_blocks_all_all ON blocks_all(direction, gpart2, gpart1, sp2, sp1)");
	$dbh->do("CREATE INDEX idx_blocks_all_block ON blocks_all(blockid)");
}

sub regenerate_breaks_all
{
	my ($dbh) = @_;
	$dbh->do('DROP TABLE IF EXISTS breaks_all');
	$dbh->do("CREATE TABLE breaks_all AS SELECT * FROM breaks_all_view");
	$dbh->do("CREATE INDEX idx_breaks_all_sp ON breaks_all(sp1, sp2)");
	$dbh->do("CREATE INDEX idx_breaks_all_break ON breaks_all(breakid)");
}

sub get_genomes_parts
{
	my ($dbh) = @_;
	my %genomes = ();
	my $query = "SELECT sp, count(gpart) AS count FROM genome_parts GROUP BY sp";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	while(my $row = $sth->fetchrow_hashref) {
		$genomes{$row->{'sp'}} = $row->{'count'};
	}
	return \%genomes;
}

sub get_common_orthos
{
	my ($dbh) = @_;
	my %orthos = ();
	my $query = "SELECT sp1, sp2, count(*) AS count FROM orthos_all GROUP BY sp1, sp2";
	my $sth = $dbh->prepare($query);
	$sth->execute();
	while(my $row = $sth->fetchrow_hashref) {
		$orthos{ $row->{'sp1'} }{ $row->{'sp2'} } = $row->{'count'};
	}
	return \%orthos;
}

sub find_couples
{
	my ($genomes, $common) = @_;
	
	# Separate fragmented from complete genomes
	my %fragmented = ();
	my @models = ();
	foreach my $g (keys %$genomes) {
		if ($genomes->{$g} > 1) {
			$fragmented{$g} = $common->{$g};
		} else {
			push @models, $g;
		}
	}
	
	my %couples = ();
	foreach my $f (keys %fragmented) {
		#print STDERR "Finding model for $f\n";
		my $best = find_best_match($fragmented{$f}, \@models);
		$couples{$f} = $best;
	}
	
	return \%couples;
}

sub find_best_match
{
	my ($counts, $models) = @_;
	
	# Find the best model genome (max orthos)
	my $max = 0;
	my $best = '';
	foreach my $g (@$models) {
		if ($counts->{$g} > $max) {
			$max = $counts->{$g};
			$best = $g;
		}
	}
	return $best;
}

sub order_all_genomes
{
	my ($dbh) = @_;
	
	# First get the list of genomes with their number of parts
	my $genomes = get_genomes_parts($dbh);
	# Second get the number of common orthos between each pair of genomes
	my $common = get_common_orthos($dbh);
	# Finally, find the closest complete genome for each fragmented one and use it as model
	my $couples = find_couples($genomes, $common);
	my $n = 0;
	foreach my $sample (sort keys %$couples) {
		my $model = $couples->{$sample};
		print STDERR "Order $sample based on $model\n";
		order_genome($dbh, $model, $sample);
		$n++;
	}
	print STDERR "$n genomes ordered\n";
	return $n;
}

sub order_genome
{
	my ($dbh, $model, $sample) = @_;
	
	# First retrieve the genes from the draft, and the orthos from the draft to the model
	my $draft = get_genes($dbh, $sample);
	my $orthos = get_orthos($dbh, $sample, $model);

	# Second, put all the draft genes in separate genome_part hashes
	my $draft_parts = partition($draft);

	# Then for each of these parts, compute the median of the pnum_all of the ortholog in the model
	my $stats = get_parts_data($draft_parts, $orthos);

	# Finally change all the pnum_all in genes of the draft genome to the order defined
	order_parts($dbh, $draft_parts, $stats);
	order_gparts($dbh, $sample);
	return 1;
}

##################################
# MAIN
init();

# Connect
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt{d}","","");

my $changed = 0;
if ($opt{a}) {
	$changed = order_all_genomes($dbh);
} else {
	$changed = order_genome($dbh, $opt{m}, $opt{s});
}

# Regenerate the blocks and breaks tables from views
if ($changed > 0) {
	regenerate_blocks_all($dbh);
	regenerate_breaks_all($dbh);
}
__END__

