#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use DBI;

our %opt;	# Getopt options

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script extends the synteny blocks in a syntebase database.
	
	-h        : this (help) message
	
	Input:
	-d <path> : path to the sqlite3 database
	-t <int>  : pairing tolerance: 0=no tolerance >0: size of pairs
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hd:t:', \%opt) or usage();
	usage() if $opt{h};
	usage("Sqlite3 database needed (-d)") unless $opt{d};
}

##################################
#
sub create_pairs
{
	my ($dbh, $table) = @_;

	# First create the pairs table
	# Uses a request adapted from syntebase
	print STDERR "Create $table table... \n";
	$dbh->do("DROP TABLE IF EXISTS $table");
	my $create_pairs_table = <<"REQ";
	CREATE TABLE $table (
		pairid INTEGER PRIMARY KEY,
		pair_order1 INTEGER,
		pair_order2 INTEGER,
		oid_start REFERENCES orthos(oid),
		oid_end REFERENCES orthos(oid),
		direction INTEGER,
		inblocks1 INTEGER,
		inblocks2 INTEGER
	)
REQ
	$dbh->do($create_pairs_table);
}

sub populate_pairs
{
	my ($dbh, $table, $tolerance) = @_;
	$tolerance = 0 if not defined($tolerance) or $tolerance < 0;
	print STDERR "Populate $table table (tolerance $tolerance) \n";
	my $populate_pairs_table = <<"REQ";
	INSERT INTO $table (oid_start, oid_end, direction, pair_order1, pair_order2, inblocks1, inblocks2)
	SELECT
		ostart.oid AS oid_start,
		oend.oid AS oid_end,
		(oend.pnum_order2 - ostart.pnum_order2) AS direction,
		oend.pnum_CDS1 AS pair_order1,
		oend.pnum_CDS2 AS pair_order2,
		abs(oend.pnum_CDS1 - ostart.pnum_CDS1) AS inblocks1,
		abs(oend.pnum_CDS2 - ostart.pnum_CDS2) AS inblocks2
	FROM
		orthos_all AS ostart,
		orthos_all AS oend
	WHERE
		    ostart.sp1 = oend.sp1
		AND ostart.sp2 = oend.sp2
		AND ostart.gpart1 = oend.gpart1
		AND ostart.gpart2 = oend.gpart2
	/* The juicy part: only orthos that are next to each other in both genomes */
	/* (the direction in the second genome may be different, although the order is the same!) */
	AND oend.pnum_CDS1   > ostart.pnum_CDS1
	AND oend.pnum_CDS1   < ostart.pnum_CDS1 + 2 + $tolerance
	AND oend.pnum_order1 = ostart.pnum_order1 + 1
	AND (
		(
		    oend.pnum_CDS2   > ostart.pnum_CDS2
		AND oend.pnum_CDS2   < ostart.pnum_CDS2 + 2 + $tolerance
		AND oend.pnum_order2 = ostart.pnum_order2 + 1
		)
		OR
		(
		    oend.pnum_CDS2   < ostart.pnum_CDS2
		AND oend.pnum_CDS2   > ostart.pnum_CDS2 - 2 - $tolerance
		AND oend.pnum_order2 = ostart.pnum_order2 - 1
		)
	)
	ORDER BY oend.sp1, oend.sp2, oend.pid1, oend.pid2
REQ
	$dbh->do($populate_pairs_table);
}

sub mark_notinpairs
{
	my ($dbh, $table) = @_;

	# Foreach genome couple, mark the orthologs that are not in a pair
	my $mark = <<"REQ";
	UPDATE orthos SET noblock=1
	WHERE oid IN (
		SELECT oid FROM orthos_all
		WHERE
		   oid NOT in (SELECT oid_start FROM $table)
		OR oid NOT in (SELECT oid_end   FROM $table)
	)
REQ
	$dbh->do($mark);
}

sub create_blocks
{
	my ($dbh, $table) = @_;
	print STDERR "Create $table table... \n";
	$dbh->do("DROP TABLE IF EXISTS $table");
	my $create_blocks_table = <<"REQ";
	CREATE TABLE $table (
		blockid INTEGER PRIMARY KEY,
		block_order1 INTEGER,
		block_order2 INTEGER,
		oid_start REFERENCES orthos(oid),
		oid_end REFERENCES orthos(oid),
		direction INTEGER,
		block_size INTEGER
	)
REQ
	$dbh->do($create_blocks_table);
}

sub populate_blocks
{
	my ($dbh, $pairs_table, $blocks_table) = @_;

	# Retrieve the different genomes names
	print STDERR "Populate blocks table... \n";
	my $genomes = get_genomes($dbh);
	
	foreach my $genome1 (sort @$genomes) {
		foreach my $genome2 (sort @$genomes) {
			next if $genome1 eq $genome2;
			print STDERR sprintf("%-12s\t%-12s", $genome1, $genome2);
			
			# Get the neighbours for those genomes
			my $pairs = get_pairs($dbh, $pairs_table, $genome1, $genome2);
			
			# Aggregate the neighbours to create blocks
			my $blocks = extend_blocks($pairs);
			
			# Add sort keys for blocks
			$blocks = sort_blocks($blocks);
			
			# Add those blocks in the table "syntenic_blocks"
			add_to_table($dbh, $blocks_table, $blocks, $genome1, $genome2);
			print STDERR "\n";
		}
	}
}

sub create_blocks_all
{
	my ($dbh) = @_;

	# Create a view to ease the searches
	$dbh->do("DROP VIEW IF EXISTS blocks_all_view");
	my $create_blocks = <<"ADD";
	CREATE VIEW blocks_all_view AS SELECT
		blockid,
		oid_start,
		oid_end,
		block_order1,
		block_order2,
		block_size,
		BLOCK.direction,
		S1.sp AS sp1,
		S2.sp AS sp2,
		S1.gpart AS gpart1,
		S2.gpart AS gpart2,
		S1.pid AS start1,
		E1.pid AS end1,
		S2.pid AS start2,
		E2.pid AS end2,
		
		S1.pnum_CDS AS pnum_CDS_start1,
		E1.pnum_CDS AS pnum_CDS_end1,
		S2.pnum_CDS AS pnum_CDS_start2,
		E2.pnum_CDS AS pnum_CDS_end2,
		
		S1.pnum_all AS pnum_all_start1,
		E1.pnum_all AS pnum_all_end1,
		S2.pnum_all AS pnum_all_start2,
		E2.pnum_all AS pnum_all_end2,
		
		S1.pnum_display AS pnum_display_start1,
		E1.pnum_display AS pnum_display_end1,
		S2.pnum_display AS pnum_display_start2,
		E2.pnum_display AS pnum_display_end2
		
	FROM blocks BLOCK, orthos START, orthos END, genes S1, genes E1, genes S2, genes E2 WHERE
		BLOCK.oid_start = START.oid
	AND	BLOCK.oid_end=END.oid
	AND S1.pid = START.pid1
	AND E1.pid = END.pid1
	AND S2.pid = START.pid2
	AND E2.pid = END.pid2
	ORDER BY sp1, sp2, pnum_all_start1
ADD
	$dbh->do($create_blocks);
	regenerate_blocks_all($dbh);
}

sub regenerate_blocks_all
{
	my ($dbh) = @_;
	$dbh->do("DROP TABLE IF EXISTS blocks_all");
	my $create_blocks = "CREATE TABLE blocks_all AS SELECT * FROM blocks_all_view";
	$dbh->do($create_blocks);
	$dbh->do("CREATE INDEX idx_blocks_all_all ON blocks_all(direction, gpart2, gpart1, sp2, sp1)");
	$dbh->do("CREATE INDEX idx_blocks_all_block ON blocks_all(blockid)");
}

sub get_genomes
{
	my ($dbh) = @_;
	
	my $query = 'SELECT DISTINCT sp FROM orthos O, genes P WHERE O.pid1=P.pid';
	my $genomes = $dbh->selectcol_arrayref($query) or die("'$query': $!");
	
	return $genomes;
}

sub get_pairs
{
	my ($dbh, $table, $genome1, $genome2) = @_;
	
	my $npairs = 0;
	
	print STDERR "\tPairs: ";
	my $select = 'SELECT pairid, oid_start, oid_end, direction, pair_order1, pair_order2, G1.gpart AS gpart1, G2.gpart AS gpart2';
	my $from = "FROM $table P, orthos O, genes G1, genes G2";
	my $where = 'WHERE P.oid_start=O.oid AND o.pid1=G1.pid AND o.pid2=G2.pid AND G1.sp="' . $genome1 . '" AND G2.sp="' . $genome2 . '"';
	my $order = 'ORDER BY pair_order1';
	my $query = "$select $from $where $order";
	
	my $pairs = $dbh->selectall_arrayref($query, { Slice => {} });
	$npairs += @$pairs;
	
	print STDERR sprintf("%6s", $npairs);
	
	return $pairs;
}

sub extend_blocks
{
	my ($neighbour_pairs) = @_;
	
	my $nturns = 0;
	
	print STDERR "\tBlocks: ";
	
	# Put those data in 3 hashes
	my %pairs = ();
	my %starts = ();
	my %ends = ();
	foreach my $pair (@$neighbour_pairs) {
		my $pair_id = $pair->{'pairid'};
		$pairs{ $pair_id } = $pair;
		# Check unicity?
		$starts{ $pair->{'oid_start'} }	= $pair;
		$ends{ $pair->{'oid_end'} }		= $pair;
	}
	
	# Extend each pair into blocks
	my @blocks = ();
	foreach my $id (sort keys %pairs) {
		my $pair = $pairs{ $id };
		
		# Already added to another block
		next if not defined($pair);
		
		# Start this new block
		my %block = %$pair;
		$block{'block_size'} = 2;
		delete($block{ $id });
		
		# Clean the hashes
		delete($pairs{ $id });
		delete($starts{ $block{'oid_start'} });
		delete($ends{ $block{'oid_end'} });
		
		# Extend the start to the left:
		# If we find another pair which ends with the same orthologs than the current block start,
		# then we can append it to the start
		while($ends{ $block{'oid_start'} }) {
			my $added_pair = $ends{ $block{'oid_start'} };
			$block{'oid_start'} = $added_pair->{'oid_start'};
			
			# Clean the hashes, the added pair can't be used anymore
			delete($pairs{ $added_pair->{'pairid'} });
			delete($starts{ $added_pair->{'oid_start'} });
			delete($ends{ $added_pair->{'oid_end'} });
			$block{'block_size'}++;
		}
		
		# Extend the end to the right:
		# If we find another pair which starts with the same orthologs than the current block end,
		# then we can append it to the end
		while($starts{ $block{'oid_end'} }) {
			my $added_pair = $starts{ $block{'oid_end'} };
			$block{'oid_end'} = $added_pair->{'oid_end'};
			
			# Clean the hashes, the added pair can't be used anymore
			delete($pairs{ $added_pair->{'pairid'} });
			delete($starts{ $added_pair->{'oid_start'} });
			delete($ends{ $added_pair->{'oid_end'} });
			$block{'block_size'}++;
		}
		
		# Add this block to the new list
		push @blocks, \%block;
	}
	
	my $nb = @blocks;
	print STDERR "$nb";
	
	return \@blocks;
}

sub sort_blocks
{
	my ($blocks) = @_;

	# Divide by gpart1
	my %gparts1 = ();
	foreach my $bl (@$blocks) {
		push @{ $gparts1{ $bl->{ 'gpart1' } } }, $bl;
	}
	my $new_blocks = [];
	foreach my $gp (keys %gparts1) {
		push @$new_blocks, @{ sort_blocks_gparts($gparts1{ $gp }, 1) };
	}
	$blocks = $new_blocks;

	# Divide by gpart2
	my %gparts2 = ();
	foreach my $bl (@$blocks) {
		push @{ $gparts2{ $bl->{ 'gpart2' } } }, $bl;
	}
	$new_blocks = [];
	foreach my $gp (keys %gparts2) {
		push @$new_blocks, @{ sort_blocks_gparts($gparts2{ $gp }, 2) };
	}

	return $new_blocks;
}

sub sort_blocks_gparts
{
	my ($blocks, $num) = @_;
	
	# Initial order (conserve it)
	for (my $i = 0; $i < @$blocks; $i++) {
		$blocks->[$i]->{'rank'} = $i;
	}
	
	# Order of the blocks on the genome $num
	@$blocks = sort { $a->{'pair_order'.$num} <=> $b->{'pair_order'.$num} } @$blocks;
	for (my $i = 0; $i < @$blocks; $i++) {
		$blocks->[$i]->{'block_order'.$num} = $i+1;	# Integer > 0
	}
	
	@$blocks = sort { $a->{'rank'} <=> $b->{'rank'} } @$blocks;
	
	return $blocks;
}

sub create_table
{
	my ($dbh, $table, $columns) = @_;
	
	# RECREATE THE TABLE
	$dbh->do("DROP TABLE IF EXISTS $table");
	$dbh->do("CREATE TABLE $table ($columns)");
	
	# Also create a schema?
}

sub add_to_table
{
	my ($dbh, $table, $data, $genome1, $genome2) = @_;
	
	# NB : gpart needs to be implemented
	my $gpart1 = $genome1;
	my $gpart2 = $genome2;
	
	# Add the data to the table
	#print STDERR "Add the data to the table synteny_blocks...\n";
	$dbh->do("BEGIN");
	
	my @cols = qw(oid_start oid_end direction block_size block_order1 block_order2);
	
	foreach my $n (@$data) {
		my $insert = "INSERT INTO $table (".join(', ', @cols).")";
		my @values = map { $n->{$_} } @cols;
		my $values = 'VALUES ("' . join('", "', @values) . '")';
		my $query = "$insert $values";
		$dbh->do($query);
	}
	$dbh->do("COMMIT");
}

##################################
# MAIN
init();

# Connect
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt{d}","","");

create_pairs($dbh, 'pairs');
populate_pairs($dbh, 'pairs', $opt{t});
mark_notinpairs($dbh, 'pairs');
create_blocks($dbh, 'blocks');
populate_blocks($dbh, 'pairs', 'blocks');
create_blocks_all($dbh, 'blocks', 'blocks_all');

__END__
