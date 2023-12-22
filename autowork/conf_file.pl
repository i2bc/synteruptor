#!/usr/bin/perl
use strict;
use JSON;

my $file = $ARGV[0];
my $key = $ARGV[1];
my $val = $ARGV[2];

my $content;
if ($file and $key) {
	# Read content of file first if it exists
	my $content = {};
	if (-e $file) {
		open(FILE, "$file") or die("Couldn't read $file: $!");
		my $JScontent = <FILE>;
		close(FILE);
		#chomp($JScontent);
		#print STDERR "$file:\n";
		#print STDERR "$JScontent\n";
		$content = decode_json($JScontent);
	}
	# Change value if provided
	if ($val) {
		$content->{$key} = $val;
		open(OUT, ">$file") or die ("Couldn't write $file: $!");
		print OUT encode_json($content);
	# Otherwise, print value
	} else {
		print STDOUT $content->{$key};
	}
}
