#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Data::Dumper;
our %opt;	# Getopt options

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script creates a gbk file from fasta, glimmer3, trna predictions.
	
	usage: $0 [-h] -f infile.fna -i infile.predict -o outfile.gbk
	
	-h        : this (help) message
	
	Input:
	-f <path> : path to the nucleotidic fasta sequence
	-i <path> : path to the Glimmer3 prediction file
	-t <path> : path to the trnascan-SE prediction file
	-a <path> : path to the annotation file (pair pid\tannotation)
	
	Output :
	-o <path> : gbk genome file
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hi:f:t:o:a:', \%opt) or usage();
	usage() if $opt{h};
	usage("Glimmer3 prediction file needed (-i)") unless $opt{i};
	usage("Nucleotide sequence fasta file needed (-f)") unless $opt{f};
	usage("Output genbank file needed (-o") unless $opt{o};
}

sub get_annotations
{
	my ($path) = @_;
	my %annot = ();
	return \%annot if not $path;

	open(ANNOT, $path) or die("Couldn't read $path: $!\n");
	while(my $line = <ANNOT>) {
		my ($pid, $annotation_text) = split(/\t/, $line);
		$annot{ $pid } = $annotation_text;
	}
	close(ANNOT);
	return \%annot;
}

sub add_glimmer_features
{
	my ($seq, $glpath, $annot) = @_;

	my $rep = '';
	my %feats = ();
	my $n = 0;
	my $first = 0;
	my $gbok = 0;
	
	open(GL, $glpath) or die("Couldn't read $glpath: $!");
	LINE : while(my $line = <GL>) {
		chomp($line);
		$line =~ s/\r//g;
		if ($line =~ /^>(.+?)($|\s)/) {
			$rep = $1;
			$first = 1;
			if ($rep eq $seq->id) {
				$gbok = 1;
			} else {
				$gbok = 0;
			}
		} elsif ($gbok) {
			chomp($line);
			my ($orf, $start, $end, $phase, $score) = split(/\s+/, $line);
			if ($phase > 0 and ($start > $end)
				or
			    $phase < 0 and ($start < $end)
				) {
				print STDERR "No circular genes allowed (after $rep - $n)\n";
				next LINE;
			}
			$n++;
			my $locus_tag = $rep . '_' . sprintf("%0.5d", $n);
			my $product = $annot->{ $locus_tag } ? $annot->{ $locus_tag } : undef;

			my $feat = new Bio::SeqFeature::Generic(
				-start => $start,
				-end => $end,
				-primary_tag => 'CDS',
				-tag => {
					'locus_tag' => $locus_tag
				}
			);
			$feat->add_tag_value("product", $product) if $product;
			$seq->add_SeqFeature($feat);
		}
	}
	print STDERR "$n glimmer annotations found\n";
	close(GL);
	return $seq;
}

sub add_trnascan_features
{
	my ($seq, $trpath) = @_;
	return $seq if not $trpath;

	my $n = 0;
	my $start = 0;
	open(TR, $trpath) or die("Couldn't read $trpath: $!");
	LINE : while(my $line = <TR>) {
		if (not $start and $line =~ /^--/) {
			$start = 1;
		} elsif ($start) {
			my ($rep, $orf, $start, $end, $type, $codon) = split(/\s+/, $line);
			if ($rep eq $seq->id) {
				$n++;
				my $feat = new Bio::SeqFeature::Generic(
					-start => $start,
					-end => $end,
					-primary_tag => 'tRNA',
					-tag => {
						'locus_tag' => $rep . '_t' . sprintf("%0.2d", $n),
						'product' => 'tRNA-' . $type . '-' . $codon
					}
				);
				$seq->add_SeqFeature($feat);
			}
		}
	}
	print STDERR "$n tRNAs found\n";
	close(TR);
	return $seq;
}

##################################
# MAIN
init();

# Get the fasta file
my $seqfasta = Bio::SeqIO->new(
	-file   => "<$opt{f}",
	-format => 'fasta',
);

my $gbk = Bio::SeqIO->new(
		-file	=> ">$opt{o}",
		-format	=> 'genbank',
);

my $annot = get_annotations($opt{a});

while (my $inseq = $seqfasta->next_seq) {
	print STDERR $inseq->id . "\n";	
	my $seq = add_glimmer_features($inseq, $opt{i}, $annot);
	$seq = add_trnascan_features($seq, $opt{t});
	$gbk->write_seq($seq);
}

__END__

