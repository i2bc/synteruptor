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
	
	This script analyses breaks in order to rank them.
	
	-h        : this (help) message
	
	Input:
	-d <path> : path to the sqlite3 database
	-C        : also clean bad breaks (e.g. filled with non-orthologs)
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hd:C', \%opt) or usage();
	usage() if $opt{h};
	usage("Sqlite3 database needed (-d)") unless $opt{d};
}

sub create_ranking_table
{
	my ($dbh) = @_;
	$dbh->do("DROP TABLE IF EXISTS breaks_ranking");
	my $create_ranking_table = <<"REQ";
	CREATE TABLE breaks_ranking(
		breakid INT REFERENCES breaks(breakid) ON DELETE CASCADE,
		real_size1 INT DEFAULT 0,
		real_size2 INT DEFAULT 0,
		tRNA_both INT DEFAULT 0,
		tRNA_both_ext INT DEFAULT 0,
		content1 TEXT,
		content2 TEXT,
		paralogs1 INT DEFAULT 0,
		paralogs2 INT DEFAULT 0,
		delta_GC1 INT,
		delta_GC2 INT,
		UNIQUE(breakid)
	)
REQ
	$dbh->do($create_ranking_table);
}

sub create_graph_table
{
	my ($dbh) = @_;
	$dbh->do("DROP TABLE IF EXISTS breaks_graph");
	my $create_graph_table = <<"REQ";
	CREATE TABLE breaks_graph(
		graphid TEXT,
		from_name TEXT,
		to_name TEXT
	)
REQ
	$dbh->do($create_graph_table);
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
		
	my $select = 'SELECT breakid, left_block, right_block, break_size1, break_size2, direction, left1, right1, opposite, sp1, sp2';
	my $from = 'FROM breaks_all';
	my $where = '';
       	if ($genome1 and $genome2) {
		$where .= 'WHERE sp1="' . $genome1 . '" AND sp2="' . $genome2 . '"';
	}
	my $order = 'ORDER BY left_block';
	my $query = "$select $from $where $order";
	my $breaks = $dbh->selectall_arrayref($query, { Slice => {} });
	
	return $breaks;
}

sub get_breaks_hash
{
	my ($dbh, $genome1, $genome2) = @_;
	my $array = get_breaks($dbh, $genome1, $genome2);
	my $hash = {};
	foreach my $break (@$array) {
		$hash->{ $break->{'breakid'} } = $break;
	}
	return $hash;
}

sub get_break_genes
{
	my ($dbh, $breakid, $num) = @_;

	my $query = "SELECT * FROM genes WHERE pid IN (SELECT pid FROM breaks_genes WHERE breakid=? AND side=?) ORDER BY pnum_all";

	my @values = (
		$breakid,
		$num
	);
	my $bgenes = $dbh->selectall_arrayref($query, { Slice => {} }, @values);

	return $bgenes;
}

sub get_breaks_genes
{
	my ($dbh, $sp1, $sp2) = @_;

	my $query = "SELECT * FROM genes G, breaks_genes B WHERE G.pid=B.pid AND sp1=? AND sp2=?";
	my @values = (
		$sp1,
		$sp2
	);
	my $bgenes = $dbh->selectall_arrayref($query, { Slice => {} }, @values);

	my $breaks_genes = {};
	foreach my $gene (@$bgenes) {
		push @{ $breaks_genes->{ $gene->{'breakid'} }->{ $gene->{'side'} }}, $gene;
	}

	return $breaks_genes;
}

sub rank_genomes_breaks
{
	my ($dbh, $clean) = @_;
	# Retrieve the different genomes names
	my $genomes = get_genomes($dbh);
	my @breaks_to_delete = ();

	foreach my $genome1 (sort @$genomes) {
		my ($n, $d) = (0, 0);
		print STDERR "$genome1\t";
		foreach my $genome2 (sort @$genomes) {
			next if $genome1 eq $genome2;
			
			# Compute ranking scores for each break
			my ($stats, $to_delete) = rank_breaks($dbh, $genome1, $genome2, $clean);
			$n += $stats->{'ranked'};
			$d += $stats->{'deleted'};
			push @breaks_to_delete, @$to_delete;
		}
		print STDERR "$n breaks ranked\t$d bad breaks to delete\n";
	}
	print STDERR "Deleting bad breaks... ";
	delete_breaks($dbh, \@breaks_to_delete);
	print STDERR "done\n";
}

sub rank_breaks
{
	my ($dbh, $sp1, $sp2, $clean) = @_;
	
	print STDERR ".";
	my $breaks = get_breaks_hash($dbh, $sp1, $sp2);
	my $breaks_genes = get_breaks_genes($dbh, $sp1, $sp2);
	
	my $stats = {
		'ranked' => 0,
		'deleted' => 0,
	};
	my $to_delete = [];
	
	$dbh->do('BEGIN');
	foreach my $breakid (sort keys %$breaks) {
		# Retrieve all associated genes
		my $bgenes1 = $breaks_genes->{ $breakid }->{ 1 };
		my $bgenes2 = $breaks_genes->{ $breakid }->{ 2 };
		
		# Scores
		my $score1 = score_genes($bgenes1);
		my $score2 = score_genes($bgenes2);
		
		my %scores = (
			'breakid' => $breakid,
			'tRNA_both' => both_tRNA($score1, $score2),
			'tRNA_both_ext' => both_tRNA_ext($score1, $score2),
			'content1' => content($score1),
			'content2' => content($score2),
			'real_size1' => $score2->{'real_size'},
			'real_size2' => $score1->{'real_size'},
			'paralogs1' => $score1->{'paralogs'},
			'paralogs2' => $score2->{'paralogs'},
			'delta_GC1' => $score1->{'delta_GC'},
			'delta_GC2' => $score2->{'delta_GC'},
		);
		
		my @cols = sort keys( %scores );
		my @vals = ();
		my @valsq = ();
		foreach my $c (@cols) {
			push @vals, $scores{$c};
			push @valsq, '?';
		}
		
		# Clean it?
		if ($clean and bad_break($breaks->{ $breakid }, \%scores)) {
			$stats->{'deleted'}++;
			push @$to_delete, $breakid;
			next;
		}
		
		# Add those scores to the ranking table
		my $query = "INSERT INTO breaks_ranking(" . join(', ', @cols) . ") VALUES(" . join(', ', @valsq) . ")";
		$dbh->do($query, undef, @vals) or die( "Double $breakid\n" );
		$stats->{'ranked'}++;
	}
	$dbh->do('COMMIT');
	
	return $stats, $to_delete;
}

sub delete_breaks
{
	my ($dbh, $breaks) = @_;
	
	$dbh->do('BEGIN');
	my $sth_delete = $dbh->prepare( "DELETE FROM breaks WHERE breakid=?" );
	foreach my $breakid (@$breaks) {
		$sth_delete->execute( $breakid );
	}
	$dbh->do('COMMIT');
	# Rebuild breaks_all, in case some breaks were deleted
	regenerate_breaks_all($dbh);
}

sub bad_break
{
	my ($break, $score) = @_;
	
	if (not defined($score->{ 'real_size1' })) {
		warn Dumper($break);
		die Dumper($score);
	}
	
	# No non-ortholog gene in the break at all: not worth keeping
	if ( $score->{ 'real_size1' } == 0
			and $score->{ 'real_size2' } == 0) {
		return 1;
	}
	# No break with too much orthologs (at least 50% in each)
	if (	(
			$break->{ 'break_size1' } == 0 or
			$score->{ 'real_size1' } / $break->{ 'break_size1' } <= 0.5 or
			$score->{ 'real_size1' } <= 2
		)
		and
		(
			$break->{ 'break_size2' } == 0 or
			$score->{ 'real_size2' } / $break->{ 'break_size2' } <= 0.5 or
			$score->{ 'real_size2' } <= 2
		)
	) {
		return 1;
	}
	
	# One of the breaks is mainly composed of orthologs
	my $min_percent = 1/4;
	if ( ($break->{ 'break_size1' } > 4 and $score->{ 'real_size1' } / $break->{ 'break_size1' } <= $min_percent)
			or ($break->{ 'break_size2' } > 4 and $score->{ 'real_size2' } / $break->{ 'break_size2' } <= $min_percent)) {
			return 1;
	}
	
	return 0;
}

sub both_tRNA
{
	my ($score1, $score2) = @_;
	
	
	if ($score1->{'tRNA'} and $score2->{'tRNA'}) {
		return 2;
	} elsif ($score1->{'tRNA'} or $score2->{'tRNA'}) {
		return 1;
	} else {
		return 0;
	}
}

sub both_tRNA_ext
{
	my ($score1, $score2) = @_;
	
	if ($score1->{'tRNA_ext'} and $score2->{'tRNA_ext'}) {
		return 2;
	} elsif ($score1->{'tRNA_ext'} or $score2->{'tRNA_ext'}) {
		return 1;
	} else {
		return 0;
	}
}

sub content
{
	my ($score) = @_;
	
	my @line = ();
	my @allowed = qw(tRNA SM regulatory resistance transport mobile phage CRISPR);
	foreach my $feat (@allowed) {
		push @line, "$score->{$feat} $feat" if $score->{$feat};
	}
	return join(", ", @line);
}

sub score_genes
{
	my ($genes) = @_;
	
	my %score = (
		'paralogs' => 0,
		'tRNA' => 0,
		'real_size' => 0,
		'delta_GC' => 0,
	);
	my $glengths = 0;
	
	# Mobile elements score
	my $end = $genes ? @$genes - 1 : 0;
	my $i = 0;
	foreach my $g (@$genes) {
		if ($g->{'product'} =~ /(\b|-)(insertion|mobile element|integrase|excisionase|plasmid|DNA ligase|transposase|transfer protein)(\b|-)/i) {
			$score{'mobile'}++;
		}
		elsif ($g->{'product'} =~ /\b(Spd[ABCD])(\b|-)/i) {
			$score{'mobile'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(pro-?)?(phage)(\b|-)/i) {
			$score{'phage'}++;
		}
		elsif ($g->{'product'} =~ /\b(CRISPR)(-.+)?\b/i) {
			$score{'CRISPR'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(regulat|repress)(or|ory|ion)(\b|-)/i) {
			$score{'regulatory'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(transport(|er|ing)|export|permease|efflux)(\b|-)/i) {
			$score{'transport'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(resistance)(\b|-)/i) {
			$score{'resistance'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(PKS|polyketide|beta[ \-]?lactamase|penicillin|antibiotic|acyl[ \-]?carrier)(\b|-)/i) {
			$score{'SM'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(.+[cd]in|.+phenazine|penicillin)(\b|-)/i) {
			$score{'SM'}++;
		}
		elsif ($g->{'product'} =~ /(\b|-)(chitin(|ase))(\b|-)/i) {
			$score{'SM'}++;
		}
		elsif ($g->{'feat'} eq 'tRNA') {
			$score{'tRNA'}++;
			
			# Find tRNA at the beginning or end
			if ($i == 0 or $i == $end) {
				$score{'tRNA_ext'}++;
			}
			# Allow some tolerance if the break is big enough (> 10)
			if ($end > 10 and ( $i <= 2 or $i >= $end - 2 )) {
				$score{'tRNA_ext'}++;
			}
		}
		$i++;
		
		# Number of genes with paralogs
		if ( $g->{'feat'} eq 'CDS' and $g->{ 'paralogs_n' } > 0) {
			$score{ 'paralogs' }++;
			}

		# Number of genes without orthologs
		if ( $g->{'feat'} eq 'CDS' and not $g->{ 'oid' } ) {
			$score{ 'real_size' }++;
		}
		
		# GC (only for CDS, because RNAs have bias)
		if ($g->{'feat'} eq 'CDS') {
			$glengths += $g->{'loc_length'};
			$score{'delta_GC'} += $g->{'loc_length'} * $g->{'delta_GC'};
		}
	}
	$score{'delta_GC'} = $glengths > 0 ? $score{'delta_GC'} / $glengths : 0;
	
	return \%score;
}

sub find_similar_breaks
{
	my ($dbh) = @_;
	my $breaks = get_breaks($dbh);

	# Clusterize the breaks into clusters
	print STDERR "Analyze graphs...\n";
	create_graph_table($dbh);
	my $groups = get_cycles($dbh, $breaks);
	
	# Then foreach break add the number of similar breaks
	$dbh->do('ALTER TABLE breaks_ranking ADD COLUMN cycle INT DEFAULT 0');
	$dbh->do('ALTER TABLE breaks_ranking ADD COLUMN graphid INT DEFAULT 0');
	$dbh->do('UPDATE breaks_ranking SET cycle=0');
	$dbh->do('BEGIN');
	
	print STDERR "Update ranking...\n";
	my $sth = $dbh->prepare("UPDATE breaks_ranking SET cycle=?, graphid=? WHERE breakid=?");
	foreach my $g (@$groups) {
		my @group = map { $g->{'group'}->{$_}->{'breakid'} } keys %{ $g->{'group'} };
		my $cycle = $g->{'cycle'};
		my $graphid = $g->{'graphid'};
		foreach my $brid (@group) {
			$sth->execute($cycle, $graphid, $brid);
		}
	}
	$dbh->do('COMMIT');
}

sub get_cycles
{
	my ($dbh, $breaks) = @_;

	# First create the pool of all breaks
	my $pool = {};
	my $similars = {};
	foreach my $b (@$breaks) {
		$pool->{$b->{'breakid'}} = $b;
		my $key = $b->{ 'left1' } . '|' . $b->{ 'right1' };
		push @{ $similars->{ $key } }, $b;
	}
	
	# Then, retrieve all groups
	my @groups = ();
	foreach my $break (@$breaks) {
		my $breakid = $break->{'breakid'};
		next if not $pool->{ $breakid };
		my $group = {};
		push @groups, get_group($break, $pool, $similars, $group);
	}
	
	# Create nodes and unite them when they are identical (=with the exact same neighbors)
	my @cycles = ();
	my @lines = ();
	my $graphid = 0;
	print STDERR "Check for cycles ";
	foreach my $g (@groups) {
		$graphid++;
		# Unite the nodes with identical neighbors
		my $graph_nodes = unite_nodes($g);
		my $graph_size = keys %$graph_nodes;
		
		push @lines, prepare_graph_lines($graph_nodes, $graphid);
		
		# Only get nodes within circles
		my $cycle_nodes = cycle_nodes($graph_nodes);
		my $cycle_size = keys(%$cycle_nodes);
		
		# Save this cycle
		#save_cycle($dbh, $cycle_size);
		
		# Mark the basic graph informations for each break
		my %cycle = (
			'group' => $g,
			'cycle' => $cycle_size,
			'graphid' => $graphid
		);
		push @cycles, \%cycle;
	}
	print STDERR " $graphid graphs analyzed. Saving... ";
	graph_table($dbh, \@lines);
	print STDERR " done\n";
	
	return \@cycles;
}

# Recursively retrieve each member of a group around a break
sub get_group
{
	my ($break, $pool, $similars, $group) = @_;
	
	# Visit every neighbor
	my $key = $break->{ 'left1' } . '|' . $break->{'right1'};
	foreach my $sim (@{ $similars->{ $key } }) {
		my $sid = $sim->{'breakid'};
		if ($pool->{$sid}) {
			# Add to the group
			$group->{$sid} = $sim;
			# Remove from the pool
			delete $pool->{$sid};
			# Check reciprocal break
			my $recip = $sim->{'opposite'};
			next if not $recip;
			# Check the reciprocal if it is not already in the group
			if (not $group->{$recip}) {
				$group = get_group($pool->{$recip}, $pool, $similars, $group);
			}
		}
	}
	return $group;
}

sub unite_nodes
{
	my ($group) = @_;
	
	my $nodes = {};
	# First, assemble the breaks in nodes
	foreach my $bid (keys %$group) {
		my $b = $group->{$bid};
		$nodes->{$b->{'sp1'}}->{'neighbors'}->{ $b->{'sp2'} } = $b;
	}
	# Second, compare the neighbors of every node
	my %unique_nodes = ();
	foreach my $sp (keys %$nodes) {
		my $neighbors = $nodes->{$sp}->{'neighbors'};
		my $fingerprint = join(' ', sort keys %$neighbors);
		push @{ $unique_nodes{$fingerprint} }, $sp;
	}
	
	my %united_nodes = ();
	my %unames = ();
	foreach my $finger (keys %unique_nodes) {
		my @list = @{ $unique_nodes{$finger} };
		
		# Unify all names
		my $name = join(' ', sort { $a cmp $b } @list);
		my $sp1 = $list[0];
		my $neighbors = $nodes->{$sp1}->{'neighbors'};
		@{ $united_nodes{$name} } = sort keys %$neighbors;
		
		# Match the nodes to the unified name
		foreach my $sp (@list) {
			$unames{ $sp } = $name;
		}
	}
	
	# Use all the unified names for the neighbors
	foreach my $uname (keys %united_nodes) {
		my @neigh = @{ $united_nodes{$uname} };
		my %new_neigh = ();
		foreach my $sp (@neigh) {
			my $name = $unames{$sp} ? $unames{$sp} : $sp;
			$new_neigh{$uname}{ $name } = 1;
		}
		@{ $united_nodes{$uname} } = sort { $a cmp $b } keys %{ $new_neigh{$uname} };
	}
	return \%united_nodes;
}

sub cycle_nodes
{
	my ($nodes) = @_;
	
	$nodes = clone($nodes);
	
	# Delete all leaves one after the other until no leaf remains
	my $deleted = 1;
	while($deleted) {
		$deleted = 0;
		foreach my $n (keys %$nodes) {
			my @neigh = @{ $nodes->{$n} };
			@neigh = filter_deleted($nodes, @neigh);
			if (@neigh < 2) {
				delete $nodes->{$n};
				$deleted = 1;
			}
		}
	}
	return $nodes;
}

sub clone
{
	my ($ref) = @_;
	
	my %ref2 = ();
	foreach my $k (keys %$ref) {
		$ref2{$k} = $ref->{$k};
	}
	return \%ref2;
}

sub filter_deleted
{
	my ($nodes, @neigh) = @_;
	
	
	my @filtered = ();
	foreach my $n (@neigh) {
		my $sp = $n;
		if ($nodes->{$sp}) {
			push @filtered, $n;
		}
	}
	
	return @filtered;
}

sub prepare_graph_lines
{
	my ($nodes, $graphid) = @_;
	
	my @lines = ();
	# From...
	foreach my $uname (keys %$nodes) {
		# To...
		foreach my $neigh (@{$nodes->{ $uname }}) {
			my @line = ($graphid, $uname, $neigh);
			push @lines, \@line;
		}
	}
	return @lines;
}

sub graph_table
{
	my ($dbh, $lines) = @_;
	
	$dbh->do('BEGIN');
	my $sth = $dbh->prepare("INSERT INTO breaks_graph(graphid, from_name, to_name) VALUES(?,?,?)");
	foreach my $line (@$lines) {
		$sth->execute(@$line);
	}
	$dbh->do('COMMIT');
}

sub regenerate_breaks_all
{
	my ($dbh) = @_;
	$dbh->do('DROP TABLE IF EXISTS breaks_all');
	$dbh->do("CREATE TABLE breaks_all AS SELECT * FROM breaks_all_view");
	$dbh->do("CREATE INDEX idx_breaks_all_sp ON breaks_all(sp1, sp2)");
	$dbh->do("CREATE INDEX idx_breaks_all_break ON breaks_all(breakid)");
}

##################################
# MAIN
init();

# Connect
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt{d}","","");
$dbh->do('PRAGMA foreign_keys = ON');

# Regenerate breaks table in case it wasn't finished
regenerate_breaks_all($dbh);

print STDERR "Create ranking table... \n";
create_ranking_table($dbh);

print STDERR "Ranking breaks:\n";
rank_genomes_breaks($dbh, $opt{C});

# Get all breaks data (except for genes)
print STDERR "Find similar breaks:\n";
find_similar_breaks($dbh);

__END__

