#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;

our %opt;	# Getopt options
#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script crudely annotates a set of sequences based on an annotated fasta file (uniprot-like).
	
	-h        : this (help) message
	
	Input
	-i <path> : blasthit file
	-a <path> : annotation fasta file
	Output
	-o <path> : protein id - annotation file
EOF
	exit 1;
}

##################################
# Command line options processing
sub init
{
	getopts('hi:a:o:', \%opt) or usage();
	usage() if $opt{h};
	usage("Blasthits file needed (-i)") unless $opt{i};
	usage("Annotation fasta file needed (-a)") unless $opt{a};
	usage("Output file needed (-o)") unless $opt{o};
}

sub get_annotations
{
	my ($inpath) = @_;
	my %annot = ();
	open(FASTA, $inpath) or die("$inpath: $!");
	while(my $line = <FASTA>) {
		if ($line =~ /^>(.+)$/) {
			my $title = $1;
			if ($title =~ /^([^ ]+) ([^=]+) [A-Z]{2}=/) {
				my $id = $1;
				my $annot_text = $2;
				$annot{ $id } = $annot_text;
			}
		}
	}
	return \%annot;
}

sub get_blast
{
	my ($inpath) = @_;
	
	my %hits = ();
	
	# Read the blasthits
	open(BH, $inpath) or die("$inpath: $!");
	while(my $line = <BH>) {
		next if $line =~ /^#/;
		chomp($line);
		# query 	subject 	%id 	alignmentlength	mismatches 	gapopenings 	querystart	queryend	subjectstart	subjectend 	Evalue	bitscore
		my ($query, $subject, $id, $allen, $mis, $gaps, $qstart, $qend, $sstart, $send, $eval, $bits) = split(/\t/, $line);
		my %bhit = ('query' => $query, 'subject' => $subject, 'id' => $id, 'length' => $allen, 'mis' => $mis, 'gaps' => $gaps, 'eval' => $eval, 'bits' => $bits);
		push @{ $hits{ $query } }, \%bhit;
	}
	close(BH);
	return \%hits;
}

sub best_matches
{
	my ($blast) = @_;
	my @matches = ();

	ID : foreach my $id (keys %$blast) {
		my $hit = undef;
		my @hits = @{ $blast->{ $id } };

		if (@hits == 0) {
			next ID;
		} elsif (@hits == 1) {
			$hit = $hits[0];
		} else {
			# Filter
			my $min_e = 9999;
			my $best_hit = undef;
			foreach my $h (@hits) {
				if ($h->{ 'eval' } < $min_e) {
					$min_e = $h->{ 'eval' };
					$best_hit = $h;
				}
			}
			$hit = $best_hit;
		}
		# Minimum check
		if (defined($hit) and
			#$hit->{ 'id' } > 40 and
			$hit->{ 'eval' } < 1E-5
		) {
			push @matches, $hit;
		}
	}
	return \@matches;
}

sub annotate
{
	my ($annot, $matches) = @_;
	my @annotations = ();
	foreach my $hit (@$matches) {
		my $annotation_text = $annot->{ $hit->{ 'subject' } };
		if ($annotation_text) {
			if ($hit->{ 'subject' } =~ /^..\|.+?\|(.+?)$/) {
				my $uniprot_id = $1;
				$annotation_text .= " ($hit->{ 'id' }% $uniprot_id)";
			}
			my %annotation = (
				'pid' => $hit->{ 'query' },
				'annotation' => $annotation_text,
				'match' => $hit->{ 'subject' },
				'id' => $hit->{ 'id' },
			);
			push @annotations, \%annotation;
		}
	}
	return \@annotations;
}

sub print_annotated
{
	my ($outpath, $pairs) = @_;
	
	my @ortho_headers = qw(pid annotation id match);
	open(OUT, ">$outpath") or die("$outpath: $!");
	print OUT join("\t", @ortho_headers) . "\n";
	my $n = 0;
	foreach my $pair (sort { $a->{'pid'} cmp $b->{'pid'} } @$pairs) {
		$n++;
		$pair->{'oid'} = $n;
		my @line = ();
		foreach my $h (@ortho_headers) {
			push @line, $pair->{$h};
		}
		print OUT join("\t", @line) . "\n";
	}
	close(OUT);
}

##################################
# MAIN
init();

my $annot = get_annotations($opt{a});
my $nannot = keys %$annot;
my ($first) = keys %$annot;
print STDERR "$nannot annotations ($first...)\n";
my $blast = get_blast($opt{i});
my $nblast = keys %$blast; print STDERR "$nblast hits blasts\n";
my $matches = best_matches($blast);
my $nmatches = @$matches; print STDERR "$nmatches good matches\n";
my $annotated = annotate($annot, $matches);
my $nannotated = @$annotated; print STDERR "$nannotated annotations\n";
print_annotated($opt{o}, $annotated);

__END__

