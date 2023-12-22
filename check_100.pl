#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;

our %opt;	# Getopt options

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script quickly compares two genbank genomes and finds the exactly matching sequences.
	
	usage: $0 [-h] -f file
	
	-h        : this (help) message
	
	-i <path> : path to the first gbk file
	-j <path> : path to the second gbk file
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts( 'hi:j:', \%opt ) or usage();
	usage() if $opt{h};
	usage("First gbk needed (-i)") unless $opt{i};
	usage("First gbk needed (-j)") unless $opt{j};
}

sub get_seqs
{
	my ($gbk) = @_;
	
	my $range = 200;
	open(GB, "$gbk") or die("$gbk: $!");
	my ($id, $seq) = ('', '');
	my %sequences = ();
	while(my $line = <GB>) {
		if ($line =~ /^\s+CDS\s+(.+?)$/) {
			$id = $1;
		}
		elsif ($line =~ /^\s+\/translation="(.+?)"/) {
			$seq = $1;
			$seq = substr($seq, -$range, $range);
			push @{ $sequences{$seq} }, $seq if $id and $seq;
			($id, $seq) = ('', '');
		}
		elsif ($line =~ /^\s+\/translation="(.+?)$/) {
			$seq = $1;
			SEQ : while(my $seqline = <GB>) {
				if ($seqline =~ /^\s+(.+?)"$/) {
					last SEQ;
				}
				elsif ($seqline =~ /^\s+(.+?)$/) {
					$seq .= $1;
				}
				else {
					die "Invalid sequence line: $seqline\n";
				}
			}
			$seq = substr($seq, -$range, $range);
			push @{ $sequences{$seq} }, $seq if $id and $seq;
			($id, $seq) = ('', '');
		}
	}
	close(GB);
	
	return \%sequences;
}

sub matching
{
	my ($seqs1, $seqs2) = @_;
	
	my @matches = ();
	
	foreach my $s1 (keys %$seqs1) {
		if (defined($seqs2->{$s1})) {
			my @match = ($seqs1->{$s1}, $seqs2->{$s1});
			push @matches, \@match;
		}
	}
	return \@matches;
}
	
##################################
# MAIN
init();

my $seqs1 = get_seqs($opt{i});
my $n1 = keys %$seqs1;
print STDERR "$n1 sequences in genome 1 ($opt{i})\n";

my $seqs2 = get_seqs($opt{j});
my $n2 = keys %$seqs2;
print STDERR "$n2 sequences in genome 2 ($opt{j})\n";

my $matches = matching($seqs1, $seqs2);
my $num = @$matches;
print STDERR "$num exact matches\n";

__END__


