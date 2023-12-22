#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use DBI;
use Digest::SHA qw(sha1_hex);

our %opt;	# Getopt options
my $max_included_blocks = 0;

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script cleans superfluous breaks (leaving only the shorter ones).
	
	-h        : this (help) message
	
	Input:
	-d <path> : path to the sqlite3 database
	-b <int>  : maximum number of allowed (distant) blocks in breaks (default: $max_included_blocks)
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hd:b:', \%opt) or usage();
	usage() if $opt{h};
	usage("Sqlite3 database needed (-d)") unless $opt{d};
	$max_included_blocks = $opt{b} if defined($opt{b});
}

sub create_breaks_table
{
	my ($dbh) = @_;
	$dbh->do("DROP TABLE IF EXISTS breaks");
	my $create_breaks_table = <<"REQ";
	CREATE TABLE breaks(
		breakid INTEGER PRIMARY KEY,
		break_sum TEXT,
		left_block REFERENCES blocks(blockid),
		right_block REFERENCES blocks(blockid),
		direction INTEGER,
		opposite REFERENCES breaks(breakid) ON DELETE CASCADE,
		break_size1 INTEGER,
		break_size2 INTEGER,
		inblocks1 INTEGER,
		inblocks2 INTEGER
	)
REQ
	$dbh->do($create_breaks_table);
	$dbh->do("CREATE INDEX idx_breaks_size1 ON breaks(break_size1)");
}

sub populate_breaks_table
{
	my ($dbh) = @_;

	my $populate_breaks_table = <<"REQ";
	INSERT INTO breaks (left_block, right_block, direction, break_size1, break_size2, inblocks1, inblocks2)
	SELECT
	A.blockid AS left_block,
	B.blockid AS right_block,
	A.direction AS direction,
	abs(B.pnum_CDS_start2 - A.pnum_CDS_end2) - 1 AS break_size1,
	abs(B.pnum_CDS_start1 - A.pnum_CDS_end1) - 1 AS break_size2,
	abs(B.block_order1 - A.block_order1) - 1 AS inblocks1,
	abs(B.block_order2 - A.block_order2) - 1 AS inblocks2

	FROM blocks_all A, blocks_all B
	WHERE
	    A.sp1 = B.sp1
	AND A.sp2 = B.sp2
	AND A.gpart1 = B.gpart1
	AND A.gpart2 = B.gpart2

	AND A.direction = B.direction

	-- Consecutive in genome 1 (only 1 possible direction)
	AND (
		B.block_order1 > A.block_order1
		AND B.block_order1 < A.block_order1 + (2 + $max_included_blocks)
	)

	-- Consecutive in genome 2 (2 possible directions)
	AND
	(
	(
		/* Same direction */
		A.direction = 1
		AND B.block_order2 > A.block_order2
		AND B.block_order2 < A.block_order2 + (2 + $max_included_blocks)
	)

	OR
	(
		/* Opposite direction */
		A.direction = -1
		AND B.block_order2 < A.block_order2
		AND B.block_order2 > A.block_order2 - (2 + $max_included_blocks)
	)
	)
REQ
	$dbh->do($populate_breaks_table);
}

sub create_breaks_all
{
	my ($dbh) = @_;
	# Finally add a full view
	$dbh->do('DROP VIEW IF EXISTS breaks_all_view');
	my $create_breaks = <<"REQ";
	CREATE VIEW breaks_all_view AS
	SELECT
		breakid,
		break_sum,
		left_block,
		right_block,
		B.direction,
		opposite,
		break_size1,
		break_size2,
		inblocks1,
		inblocks2,
		oL1.sp AS sp1,
		oL2.sp AS sp2,
		oL1.gpart AS gpart1,
		oL2.gpart AS gpart2,

		oL1.pid AS left1,
		oR1.pid AS right1,
		oL2.pid AS left2,
		oR2.pid AS right2,
		
		oL1.pnum_CDS AS pnum_CDS_left1,
		oR1.pnum_CDS AS pnum_CDS_right1,
		oL2.pnum_CDS AS pnum_CDS_left2,
		oR2.pnum_CDS AS pnum_CDS_right2,
		
		oL1.pnum_all AS pnum_all_left1,
		oR1.pnum_all AS pnum_all_right1,
		oL2.pnum_all AS pnum_all_left2,
		oR2.pnum_all AS pnum_all_right2,
		
		oL1.pnum_display AS pnum_display_left1,
		oR1.pnum_display AS pnum_display_right1,
		oL2.pnum_display AS pnum_display_left2,
		oR2.pnum_display AS pnum_display_right2,
		
		L.block_size AS left_size,
		R.block_size AS right_size
		
	FROM breaks B, blocks L, blocks R, orthos OL, orthos orR, genes oL1, genes oR1, genes oL2, genes oR2
	WHERE
		B.left_block = L.blockid
	AND B.right_block = R.blockid
	AND L.oid_end = oL.oid
	AND R.oid_start = orR.oid
	AND oL.pid1 = oL1.pid
	AND orR.pid1 = oR1.pid
	AND oL.pid2 = oL2.pid
	AND orR.pid2 = oR2.pid
REQ
	$dbh->do($create_breaks);
	regenerate_breaks_all($dbh);
}

sub regenerate_breaks_all
{
	my ($dbh) = @_;
	$dbh->do('DROP TABLE IF EXISTS breaks_all');
	$dbh->do("CREATE TABLE breaks_all AS SELECT * FROM breaks_all_view");
	$dbh->do("CREATE INDEX idx_breaks_all_sp ON breaks_all(sp1, sp2)");
	$dbh->do("CREATE INDEX idx_breaks_all_break ON breaks_all(breakid)");
	$dbh->do("CREATE INDEX idx_breaks_all_break_sum ON breaks_all(break_sum)");
}

sub get_genomes
{
	my ($dbh) = @_;
	
	my $query = 'SELECT DISTINCT sp FROM orthos O, genes P WHERE O.pid1=P.pid';
	my $genomes = $dbh->selectcol_arrayref($query) or die("'$query': $!");
	
	return $genomes;
}

sub get_breaks
{
	my ($dbh, $genome1, $genome2) = @_;

	my $select = 'SELECT * FROM breaks_all';
	my $where = '';
	if ($genome1 and $genome2) {
		$where = 'WHERE sp1="' . $genome1 . '" AND sp2="' . $genome2 . '"';
	}
	my $order = 'ORDER BY left_block';
	my $query = "$select $where $order";
	
	my $breaks = $dbh->selectall_arrayref($query, { Slice => {} });
	return $breaks;
}

sub get_all_breaks
{
	my ($dbh) = @_;
	my $all = get_breaks($dbh);
	my $breaks = {};
	foreach my $br (@$all) {
		push @{ $breaks->{ $br->{'sp1'} }->{ $br->{'sp2'} } },  $br;
	}
	return $breaks;
}

sub skip_big_breaks
{
	my ($br, $to_remove) = @_;
	
	if ($br->{'break_size1'} >= 30 and $br->{'break_size2'} >= 100 or
			$br->{'break_size2'} >= 30 and $br->{'break_size1'} >= 100 or
			$br->{'break_size1'} >= 400 or
			$br->{'break_size2'} >= 400) {
		push @$to_remove, $br;
		return 1;
	} else {
		return 0;
	}
}

sub clean_all_breaks
{
	my ($dbh) = @_;
	
	# Retrieve the different genomes names
	my $genomes = get_genomes($dbh);
	
	# Get all current breaks
	my $all_breaks = get_all_breaks($dbh);
	
	# Recreate the table!
	create_breaks_table($dbh);
	
	my $ntotal = 0;
	foreach my $genome1 (sort @$genomes) {
		print STDERR "$genome1\t";
		foreach my $genome2 (sort @$genomes) {
			next if $genome1 eq $genome2;
			
			# Get the breaks for those genomes
			my $breaks = $all_breaks->{$genome1}->{$genome2};
			if (not $breaks) {
				print STDERR " (No break with $genome2) ";
				next;
			}
			my $n = @$breaks;
			
			# Clean the breaks to keep the shorter ones
			$breaks = clean_breaks($dbh, $breaks, 'left');
			$breaks = clean_breaks($dbh, $breaks, 'right');
			print STDERR ".";
			add_breaks($dbh, $breaks);
			my $nbr = @$breaks;
			$ntotal += $nbr;

			## TODO: remove overlapping breaks
		}
		print STDERR "\n";
	}
	print STDERR " $ntotal breaks found\n";
	# Recreate a breaks_all table with all remaining breaks
	create_breaks_all($dbh);
}

sub add_breaks
{
	my ($dbh, $breaks) = @_;
	$dbh->do('BEGIN');
	my $insert = $dbh->prepare("INSERT INTO breaks (left_block, right_block, direction, break_size1, break_size2, inblocks1, inblocks2, break_sum) VALUES(?, ?, ?, ?, ?, ?, ?, ?)");
	foreach my $br (@$breaks) {
		my @vals = (
			$br->{'left_block'},
			$br->{'right_block'},
			$br->{'direction'},
			$br->{'break_size1'},
			$br->{'break_size2'},
			$br->{'inblocks1'},
			$br->{'inblocks2'},
			$br->{'break_sum'},
		);
		$insert->execute(@vals);
	}
	$dbh->do('COMMIT');
}

sub clean_breaks
{
	my ($dbh, $breaks, $direction) = @_;
	
	my $ndels = 0;
	
	my $block = '';
	if ($direction and $direction eq 'left') {
		$block = 'left_block';
	} elsif ($direction and $direction eq 'right') {
		$block = 'right_block';
	} else {
		die("Wrong direction: must be either 'left' or 'right'");
	}
	
	# Put the breaks in a hash for every start
	my %starts = ();
	foreach my $br (@$breaks) {
		#next if (skip_big_breaks($br, \@to_remove));
		
		my $side = $br->{ $block };
		if (not $side) {
			foreach my $k (keys %$br) {
				print STDERR "$k : $br->{$k}\n";
			}
			die("No side in $direction : $side");
		}
		push @{ $starts{$side} }, $br;
	}
	
	# Prepare order: left-to-right or the opposite?
	my @order = $block eq 'left_block' ? sort { $a <=> $b } keys %starts : sort { $b <=> $a } keys %starts;
	
	# Then, for every start, only keep the shorter break
	my %keep = ();
	foreach my $st (@order) {
		my @brs = @{ $starts{$st} };
		if (@brs == 1) {
			$keep{$brs[0]->{'breakid'}} = $brs[0];
		}
		
		my $shorter = {};
		my $shorter_length = 9999999;
		
		# Find the shorter one
		foreach my $break (@brs) {
			my $length = $break->{'break_size1'} + $break->{'break_size2'};
			
			if ($length < $shorter_length) {
				# Put the new shorter one in place
				$shorter_length = $length;
				$shorter = $break;
			}
		}
		$keep{$shorter->{'breakid'}} = $shorter;
	}
	
	# Return the breaks that are kept
	my @to_keep = ();
	foreach my $br (sort keys %keep) {
		push @to_keep, $keep{$br};
	}
	return \@to_keep;
}

sub opposites
{
	my ($dbh) = @_;

	my $breaks = get_breaks($dbh);
	
	my %borders = ();
	foreach my $b (@$breaks) {
		my $id = join('|', $b->{'sp1'}, $b->{'sp2'}, $b->{'left1'}, $b->{'right1'});
		$borders{$id} = $b;
	}
	
	$dbh->do('PRAGMA foreign_keys = ON');
	$dbh->do('BEGIN');
	my $sth = $dbh->prepare("UPDATE breaks SET opposite=? WHERE breakid=?");
	my $delete = $dbh->prepare("DELETE FROM breaks WHERE breakid=?");
	my %opposites = ();
	foreach my $b (@$breaks) {
		my $id1 = join('|', $b->{'sp2'}, $b->{'sp1'}, $b->{'left2'}, $b->{'right2'});
		my $id2 = join('|', $b->{'sp2'}, $b->{'sp1'}, $b->{'right2'}, $b->{'left2'});
		
		my $opposite;
		if ($borders{$id1}) {
			$opposite = $borders{$id1};
		}
		elsif ($borders{$id2}) {
			$opposite = $borders{$id2};
		} else {
			print STDERR "Warning: no opposite found for $b->{'breakid'} (delete the break)\n";
			print STDERR "\t$id1\n\t$id2\n";
			$delete->execute($b->{'breakid'});
			next;
		}
		$sth->execute($opposite->{'breakid'}, $b->{'breakid'});
	}
	$dbh->do('COMMIT');
}

sub breaks_sums
{
	my ($dbh) = @_;

	my $breaks = get_breaks($dbh);
	$dbh->do('PRAGMA foreign_keys = ON');
	$dbh->do('BEGIN');
	my $sth = $dbh->prepare("UPDATE breaks SET break_sum=? WHERE breakid=?");
	foreach my $b (@$breaks) {
		# Create a checksum using the gene ids of the four extremities of the break
		# NB: not using sp1/sp2 since this name may change: this means that if there is a
		# genome in different versions, then breaks may have the same break_sum, but within
		# different comparisons
		# Hence, this breaks_sum is only useful for the user to map identical breaks between different computations
		my $sum = sha1_hex(join('|', $b->{'left1'}, $b->{'right1'}, $b->{'left2'}, $b->{'right2'}));
		$sth->execute($sum, $b->{'breakid'});
	}
	$dbh->do('COMMIT');
}

##################################
# MAIN
init();

# Connect
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt{d}","","");

# Create all breaks with a single SQL request
print STDERR "Create breaks table... \n";
create_breaks_table($dbh);
populate_breaks_table($dbh);
create_breaks_all($dbh);

# Post-processing: clean breaks
print STDERR "Cleaning breaks... \n";
clean_all_breaks($dbh);

# Compute the opposite of every break
opposites($dbh);
# Compute the break_sum (more permament identifier)
breaks_sums($dbh);
# Finish by regenerating the view
regenerate_breaks_all($dbh);


__END__
