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
	
	This script creates a syntebase-like database.
	
	-h        : this (help) message
	
	Input:
	-i <path> : ortholog pairs
	-p <path> : paralog pairs (optional)
	-g <path> : gene data
	-G <path> : genomes data path (optional)
	-d <path> : path to the sqlite3 database to create
	
	Additional info:
	-N <str>  : a name for this data set
	-A <str>  : an author
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hi:p:g:G:d:N:A:', \%opt) or usage();
	usage() if $opt{h};
	usage("Ortholog pairs needed (-i)") unless $opt{i};
	usage("Gene data file needed (-g)") unless $opt{g};
	usage("Sqlite3 database needed (-d)") unless $opt{d};
}

sub add_paralogs
{
	my ($dbh, $para_path) = @_;
	return if not $para_path;

	open(PARS, $para_path) or die("Couldn't read $para_path: $!");
	$dbh->do('BEGIN');
	my $sth = $dbh->prepare('UPDATE genes SET paralogs_n=?, paralogs=? WHERE pid=?');
	while(my $line = <PARS>) {
		chomp($line);
		my ($pid, $npars, $pars) = split(/\t/, $line);
		$sth->execute($npars, $pars, $pid);
	}
	$dbh->do('COMMIT');
	close(PARS);
}

sub add_genome_infos
{
	my ($dbh, $info_path) = @_;
	return if not $info_path;

	open(INFO, $info_path) or die("Couldn't read $info_path: $!");
	my $sth = $dbh->prepare('UPDATE genomes SET name=?, GC=? WHERE sp=?');
	my $head = <INFO>;
	chomp($head);
	my @heads = split(/\t/, $head);
	
	while(my $line = <INFO>) {
		chomp($line);
		my @elts = split(/\t/, $line);
		my %data = ();
		for (my $i = 0; $i < @heads; $i++) {
			$data{ $heads[$i] } = $elts[$i];
		}
		$sth->execute($data{'species'}, $data{'GC'}, $data{'abbr'});
	}
	close(INFO);
}

sub get_genomes
{
	my ($dbh) = @_;
	
	my $query = 'SELECT DISTINCT sp FROM orthos O, genes P WHERE O.pid1=P.pid';
	my $genomes = $dbh->selectcol_arrayref($query) or die("'$query': $!");
	
	return $genomes;
}

sub get_orthos
{
	my ($dbh, $sp1, $sp2) = @_;
	
	my $npairs = 0;
	
	my $select = "SELECT oid, pid1, pid2, G1.pnum_CDS AS pnum_CDS1, G2.pnum_CDS AS pnum_CDS2 FROM orthos O, genes G1, genes G2 WHERE O.pid1=G1.pid AND O.pid2=G2.pid AND G1.sp=? AND G2.sp=?";
	my @values = ($sp1, $sp2);
	my $orthos = $dbh->selectall_arrayref($select, { Slice => {} }, @values);
	my $n = @$orthos;
	print STDERR "\t$n";
	
	return $orthos;
}

sub order_orthos_group
{
	my ($orthos) = @_;

	# Need to partition in the different gparts??

	# Order along the first genome
	my $n1 = 0;
	foreach my $o (sort { $a->{'pnum_CDS1'} <=> $b->{'pnum_CDS1'} } @$orthos) {
		$n1++;
		$o->{'pnum_order1'} = $n1;
	}
	# Order along the second genome
	my $n2 = 0;
	foreach my $o (sort { $a->{'pnum_CDS2'} <=> $b->{'pnum_CDS2'} } @$orthos) {
		$n2++;
		$o->{'pnum_order2'} = $n2;
	}
	return $orthos;
}

sub update_orthos_order
{
	my ($dbh, $all_orthos) = @_;
	
	# First, create the table
	print STDERR "Add order to orthos table... \n";
	$dbh->do("ALTER TABLE orthos ADD COLUMN pnum_order1 int");
	$dbh->do("ALTER TABLE orthos ADD COLUMN pnum_order2 int");

	# Then add all the data
	$dbh->do('BEGIN');
	my $sth = $dbh->prepare('UPDATE orthos SET pnum_order1=?, pnum_order2=? WHERE oid=?');
	foreach my $o (@$all_orthos) {
		my @values = ($o->{'pnum_order1'}, $o->{'pnum_order2'}, $o->{'oid'});
		$sth->execute(@values);
	}
	$dbh->do('COMMIT');
}

sub order_orthos
{
	my ($dbh) = @_;
	
	my $genomes = get_genomes($dbh);
	
	my @all_orthos = ();
	foreach my $sp1 (@$genomes) {
		foreach my $sp2 (@$genomes) {
			next if $sp1 eq $sp2;
			print STDERR sprintf("Order %-12s vs %-12s:", $sp1, $sp2);
			my $orthos = get_orthos($dbh, $sp1, $sp2);
			$orthos = order_orthos_group($orthos);
			push @all_orthos, @$orthos;
			print STDERR "\n";
		}
	}
	print STDERR "Number of ordered orthos: " . (@all_orthos+0) . "\n";
	# Now that all genes have their order inside, time to update!
	update_orthos_order($dbh, \@all_orthos);
}

##################################
# MAIN
init();

# Connect
my $dbh = DBI->connect("dbi:SQLite:dbname=$opt{d}","","");

my $genes_path = $opt{g};
my $orthos_path = $opt{i};

# First creates the genes data table
my $genes_table = 'genes';
print STDERR "Create $genes_table table... \n";
$dbh->do("DROP TABLE IF EXISTS $genes_table");
my $create_genes_table = <<"REQ";
CREATE TABLE $genes_table (
	sp text,
	gpart text,
	pid text PRIMARY KEY,
	pnum_CDS int,
	pnum_all int,
	feat text,
	loc_start int,
	loc_end int,
	strand int,
	loc_length int,
	sequence text,
	product text,
	GC real,
	delta_GC real
);
REQ
$dbh->do($create_genes_table);

# Import genes data
print STDERR "Import data into $genes_table table... \n";
system("sqlite3 $opt{d} -separator '	' \".import $genes_path $genes_table\"") and die("sqlite3 failed. Abort.");
$dbh->do("DELETE FROM $genes_table WHERE sp='sp'");

# Add the display pnum column
$dbh->do("ALTER TABLE $genes_table ADD COLUMN pnum_display int");
$dbh->do("UPDATE $genes_table SET pnum_display=pnum_all");

# Add the paralog columns
print STDERR "Import paralog data into $genes_table table... \n";
$dbh->do("ALTER TABLE $genes_table ADD COLUMN paralogs_n int DEFAULT 0");
$dbh->do("ALTER TABLE $genes_table ADD COLUMN paralogs text");
add_paralogs($dbh, $opt{p});

# Also create a simple table for the list of genomes
print STDERR "Create genomes table... \n";
$dbh->do("DROP TABLE IF EXISTS genomes");
my $create_genomes_table = <<"REQ";
CREATE TABLE genomes (
	sp text,
	name text,
	max_pnum_display int,
	max_loc_display int,
	GC int
);
REQ
$dbh->do($create_genomes_table);

my $populate_genomes_table = <<"REQ";
INSERT INTO genomes (sp, name, max_pnum_display)
SELECT sp, sp, max(pnum_display) FROM genes GROUP BY sp ORDER BY sp
REQ
$dbh->do($populate_genomes_table);

# Add more detailed genome informations if provided
add_genome_infos($dbh, $opt{G});

# Also create a simple table for the list of genome parts
print STDERR "Create genome_parts table... \n";
$dbh->do("DROP TABLE IF EXISTS genome_parts");
my $create_genome_parts_table = <<"REQ";
CREATE TABLE genome_parts (
	sp text,
	gpart text,
	min int,
	max int
);
REQ
$dbh->do($create_genome_parts_table);

my $populate_genome_parts_table = <<"REQ";
INSERT INTO genome_parts (sp, gpart, min, max)
SELECT sp, gpart, min(pnum_display), max(pnum_display) FROM genes GROUP BY gpart ORDER BY sp, min(pnum_display)
REQ
$dbh->do($populate_genome_parts_table);

# Finally, import the orthology data table
print STDERR "Create orthos table... \n";
$dbh->do("DROP TABLE IF EXISTS orthos");
my $create_orthos_table = <<"REQ";
CREATE TABLE orthos (
	oid int PRIMARY KEY,
	pid1 text REFERENCES genes(pid),
	pid2 text REFERENCES genes(pid),
	o_ident	INTEGER,
	o_alen INTEGER
)
REQ
$dbh->do($create_orthos_table);

# Import orthos data
print STDERR "Import data into orthos table... \n";
system("sqlite3 $opt{d} -separator '	' \".import $orthos_path orthos\"") and die("sqlite3 is required. Abort.");
$dbh->do("DELETE FROM orthos WHERE oid='oid'");
$dbh->do("ALTER TABLE orthos ADD COLUMN noblock int DEFAULT 0");

# Order orthos along each genome
order_orthos($dbh);

# Create detailed orthos_all view
$dbh->do("DROP VIEW IF EXISTS orthos_all");
my $create_orthos_view = <<"REQ";
CREATE VIEW orthos_all AS
SELECT
oid,
noblock,
pnum_order1,
pnum_order2,
G1.pid AS pid1,
G2.pid AS pid2,
o_ident,
o_alen,
G1.sp AS sp1,
G2.sp AS sp2,
G1.gpart AS gpart1,
G2.gpart AS gpart2,
G1.product AS product1,
G2.product AS product2,
G1.pnum_CDS AS pnum_CDS1,
G2.pnum_CDS AS pnum_CDS2,
G1.pnum_all AS pnum_all1,
G2.pnum_all AS pnum_all2,
G1.pnum_display AS pnum_display1,
G2.pnum_display AS pnum_display2
FROM orthos O, genes G1, genes G2
WHERE
    O.pid1 = G1.pid
AND O.pid2 = G2.pid
REQ
$dbh->do($create_orthos_view);

# Create info table
create_info($dbh, $opt{N}, $opt{A});

sub create_info
{
	my ($dbh, $description, $author) = @_;
	my $genomes = get_genomes($dbh);
	my $num = @$genomes;
	
	# Prepare data to input
	my %data = (
		"description" => $description ? $description : "",
		"author" => $author ? $author : "",
		"date" => get_date(),
		"num" => $num
	);
	# Create the database
	my $create_info_table = <<"REQ";
CREATE TABLE info (
	description text,
	author text,
	date text,
	num int
)
REQ
	$dbh->do($create_info_table);
	
	# Populate the database
	my $populate_info_table = <<"REQ";
INSERT INTO info (description, author, date, num) VALUES (?, ?, ?, ?);
REQ
	my @values = ($data{description}, $data{author}, $data{date}, $data{num});
	$dbh->do($populate_info_table, undef, @values);
}

sub get_date
{
	use Time::localtime;
	my $tm = localtime;
	my $date = sprintf("%04d-%02d-%02d", $tm->year+1900, ($tm->mon)+1, $tm->mday);
	return $date;
}

print STDERR "Initial database creation done\n";

__END__

