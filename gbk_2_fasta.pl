#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
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
	
	This script extracts fasta file(s) from a gbk file.
	
	usage: $0 [-h] -i file.gbk -n file.fna
	
	-h        : this (help) message
	
	Input:
	-i <path> : gbk file

	Output:
	-n <path> : nucleotid sequence fasta file
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hi:n:', \%opt) or usage();
	usage() if $opt{h};
	usage("Genbank file needed (-i)") unless $opt{i};
	usage("Output fasta file needed (-i)") unless $opt{n};
}

sub gbk_to_fastan
{
	my ($gbk_path, $fasta_path) = @_;

	# Open the stream
	my $gbk_file;
	
	if ($gbk_path =~ /\.(dat|txt|embl)$/) {
		eval {
			$gbk_file = Bio::SeqIO->new(
					-file	=> "<$gbk_path",
					-format	=> 'embl',
			);
		};
	} elsif ($gbk_path =~ /\.(gb.*|genbank)$/) {
		eval {
			$gbk_file = Bio::SeqIO->new(
					-file	=> "<$gbk_path",
					-format	=> 'genbank',
			);
		};
	} else {
		print STDERR "Invalid file extension\n";
	}
	if ($@) {
		die("Error: unsupported format: $@");
	}
		
	# Open the fasta
	my $fasta_file = Bio::SeqIO->new(
			-file	=> ">$fasta_path",
			-format	=> 'fasta',
	);
	
	SEQ : while(my $seq = $gbk_file->next_seq) {
		$fasta_file->write_seq($seq);
	}
}

####################################
# MAIN
init();
gbk_to_fastan($opt{i}, $opt{n});

__END__
