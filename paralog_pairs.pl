#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;

our %opt;	# Getopt options
my $default_id = 40;

#################################################
# Message about this program and how to use it
sub usage {
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script creates a list of paralogs for each gene in a blast.
	
	-h        : this (help) message
	
	Input
	-i <path> : blasthits input file
	-g <path> : genes data file
	Output
	-o <path> : paralogs output file

	-s <num>  : identity threshold for paralogs (default: $default_id)
	
	-v		  : verbose
EOF
	exit 1;
}

##################################
# Command line options processing
sub init {
	getopts('hi:g:o:s:v', \%opt) or usage();
	usage() if $opt{h};
	usage("Blasthits file needed (-i)") unless $opt{i};
	usage("Genes data file needed (-g)") unless $opt{g};
	usage("Pairing file needed (-o)") unless $opt{o};
	$default_id = $opt{s} if defined( $opt{s} );
}

sub get_genes_data {
	my ($path, $protdata) = @_;
	my %protdata = ();
	my %protorder = ();
	
	open(DATA, $path) or die("$path: $!");
	
	# Get fields names
	my $head = <DATA>;
	chomp($head);
	$head =~ s/^--//;
	my @prothead = split(/\t/, $head);
	
	while(my $line = <DATA>) {
		next if $line =~ /^--/;
		chomp($line);
		my @fields = split(/\t/, $line);
		
		my %data = ();
		for (my $i = 0; $i < @prothead; $i++) {
			$data{ $prothead[$i] } = $fields[$i];
		}
		$protdata{ $data{'pid'} } = \%data;
	}
	close(DATA);
	
	return \%protdata;
}

sub para_search {
	my ($inpath, $protdata) = @_;
	
	my %hits = ();
	my $nhits = 0;
	my $nhits_all = 0;
	my $condemned = 0;
	
	# Read the blasthits
	open(BH, $inpath) or die("$inpath: $!");
	while(my $line = <BH>) {
		next if $line =~ /^#/;
		chomp($line);
		$nhits_all++;
		# query 	subject 	%id 	alignmentlength	mismatches 	gapopenings 	querystart	queryend	subjectstart	subjectend 	Evalue	bitscore
		my ($query, $subj, $id, $allen, $mis, $gaps, $qstart, $qend, $sstart, $send, $eval, $bits) = split(/\t/, $line);
		
		# Filter
		if (not defined($protdata->{ $query })) {
			warn("No protein data for '$query'") if $condemned < 10;
			warn("...\n") if $condemned == 10;
			$condemned++;
			next;
		}
		if (not defined($protdata->{ $subj })) {
			warn("No protein data for '$subj'") if $condemned < 10;
			warn("...\n") if $condemned == 10;
			$condemned++;
			next;
		}
		
		# Don't care about different genome comparisons
		if ( $protdata->{ $query }->{ 'sp' } ne $protdata->{ $subj }->{ 'sp' } or
			       	$query eq $subj ) {
			next;
		}
		my $sp = $protdata->{ $query }->{ 'sp' };
		
		my $qlen = $protdata->{ $query }->{'length'};	# /3: nucleic length!
		$qlen = $qlen/3;
		my $slen = $protdata->{ $subj  }->{'length'};
		$slen = $slen/3;
		my $smaller = $qlen < $slen ? $qlen : $slen;
		
		# Filter by length, identity and evalue
		if ($allen < $smaller * 0.5) {
			next;
		}
		elsif ($id < $default_id) {
			next;
		}
		elsif ($eval > 1E-20) {
			next;
		# KEEP THIS HIT!
		} else {
			if (not $hits{ $sp }{ $query }{ $subj } or $hits{ $sp }{ $query }{ $subj } < $id) {
				$hits{ $sp }{ $query }{ $subj } = $id;
				$nhits++;
			}
		}
	}
	close(BH);
	
	warn("$nhits / $nhits_all filtered hits\n") if $opt{v};
	
	if ($condemned) {
		die("Some sequences have no data in the proteins table: can't continue. ($condemned errors)");
	}
	return \%hits;
}

sub print_para
{
	my ($outpath, $paras) = @_;

	open(OUT, ">$outpath") or die("$outpath: $!");
	my $n = 0;
	foreach my $sp (keys %$paras) {
		foreach my $query (sort { $a cmp $b } keys %{ $paras->{ $sp } }) {
			my @line = ();
			push @line, $query;
			my @para_list = ();
			foreach my $subj (sort { $a cmp $b } keys %{ $paras->{ $sp }->{ $query } }) {
				my $id = int($paras->{ $sp }->{ $query }->{ $subj });
				my $text = "$subj ($id\%)";
				push @para_list, $text;
			}
			push @line, $#para_list+1;
			push @line, join(", ", @para_list);
			print OUT join("\t", @line) . "\n";
		}
	}
	close(OUT);
}


##################################
# MAIN
init();

my $protdata = get_genes_data( $opt{g} );
my $para = para_search( $opt{i}, $protdata );
print_para( $opt{o}, $para );

__END__

