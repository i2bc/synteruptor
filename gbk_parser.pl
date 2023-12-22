#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under GPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
our %opt;	# Getopt options

# Saved fields
my @genes_fields = qw( sp gpart pid pnum_CDS pnum_all feat loc_start loc_end strand length sequence product GC delta_GC );
my @genomes_fields = qw( abbr species strain taxonomy GC );
my @gpart_fields = qw( abbr gpart definition date accession version nbbasepair nbprotein );

#################################################
# Message about this program and how to use it
sub usage
{
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script quickly parses a Genbank genome file.
	The gbk filename will be used as the abbreviation for this species.
	
	usage: $0 [-h] -i file.gbk -f file.faa -o file.genes -g file.genome
	
	-h        : this (help) message
	
	Input:
	-i <path> : path to the gbk file to read.
				IF this parameter is a file descriptor, then all files will be concatenated in one output each,
				except for fasta files
	
	Output :
	-f <path> : fasta file with all the CDS (amino-acid sequences).
	-o <path> : tabulated file with all genes data (one gene per line)
	-g <path> : tabulated file with the genome data (one genome per line)
	-G <path> : tabulated file with the genome parts data (one gpart per line)
	
	Additionnal:
	-p <str>  : prefix for locus_tags (default locus_tags used: LOCUS display name in the gbk sequence headers).
	-v <str>  : verbose
EOF
	exit 1;
}

##################################
# Command line options processing
sub init()
{
	getopts('hi:f:o:g:G:p:v', \%opt) or usage();
	usage() if $opt{h};
	usage("Genbank file needed (-i)") unless $opt{i};
}

# GBK parser, using Bioperl
sub parse_gbk
{
	my ($gbk_paths, $prefix) = @_;
	
	my %genes = ();
	my @genomes = ();
	my @gparts = ();
	my @sequences = ();
	
	foreach my $gbk_path (glob(`ls $gbk_paths`)) {
		print STDERR "Parsing $gbk_path...\n";
		
		# Open the stream
		my $gbk_file;
		
		if ($gbk_path =~ /\.(dat|embl|txt)$/) {
			eval {
				$gbk_file = Bio::SeqIO->new(
						-file	=> "<$gbk_path",
						-format	=> 'embl',
				);
			};
		# Default to gbk
		} else {
			eval {
				$gbk_file = Bio::SeqIO->new(
						-file	=> "<$gbk_path",
						-format	=> 'genbank',
				);
			};
		}
		if ($@) {
			die("Error: unsupported format: $@");
		}
		
		# Deduce genome species abbreviation from the gbk filename (easier to customize without modifying the gbk file)
		my $abbr = '';
		if ($gbk_path =~ /([^\/]+)\.[a-z]{2,3}/) {
			$abbr = $1;
		} else {
			warn("Can't guess the genome abbreviation from the file name $gbk_path");
		}
		
		# Read each sequence entry (there may be several entries, for several gparts or several contigs)
		my %gene = ();
		my %genome = (
			'abbr' => $abbr,
			'species' => $abbr,
			'strain' => '',
			'taxonomy' => '',
			'GC' => 0,
			'prot_lengths' => 0
		);
		
		SEQ : while(my $seq = $gbk_file->next_seq) {
			# First, read the genome and gpart data
			my %gpart = ();
			my $rname = $seq->display_name;
			$gpart{'gpart'} = $rname;
			$gpart{'abbr'} = $abbr;
			
			if (defined($seq->species)) {
				if (defined($seq->species->node_name)) {
					$genome{'species'} = $seq->species->node_name;
				} if (defined($seq->species->classification)) {
					$genome{'taxonomy'} = join(',', $seq->species->classification);
				}
			}
			
			$gpart{'definition'} = $seq->desc;
			$gpart{'date'} = join(',', $seq->get_dates);
			$gpart{'accession'} = $seq->accession;
			$gpart{'version'} = $seq->seq_version;
			$gpart{'nbbasepair'} = $seq->length;
			
			# Read every gene/CDS etc.
			my $num = 0;
			my $ncds = 0;
			FEAT : foreach my $feat ($seq->get_SeqFeatures) {
				# Other data on the gpart and genome
				if ($feat->primary_tag eq 'source') {
					if ($feat->has_tag('serovar')) {
						$gpart{'strain'} = join(',', $feat->get_tag_values('serovar'));
						$genome{'strain'} = $gpart{'strain'};
					}
					next FEAT;
					#if ($feat->has_tag('mol_type')) {
					#	$gpart{'gpart'} = join(',', $feat->get_tag_values('mol_type'));
					#}
				} else {
					# Non source = feature: give it a unique number
					$num++;
				}
				
				# Explicit locus_tag for this feature?
				my $locus = '';
				if ($feat->has_tag('locus_tag')) {
					($locus) = $feat->get_tag_values('locus_tag');
					
				} else {
					# No locus_tag? Create one from the gpart name (and a prefix if provided)
					my $pnum = sprintf("%05d", $num);
					if (not defined($prefix)) {
						$locus = $rname.'_'.$pnum;
					} else {
						$locus = $prefix.$rname.'_'.$pnum;
					}
				}
				
				my %rgene = ();
				
				# Save main genome data
				$rgene{'sp'} = $abbr;
				$rgene{'gpart'} = $gpart{'gpart'};
				
				# Save data for this locus
				$rgene{'pid'} = $locus;
				$rgene{'feat'} = $feat->primary_tag;
				$rgene{'product'} = $feat->has_tag('product') ? join(',', $feat->get_tag_values('product')) : '';
				$rgene{'product'} =~ s/"//g;
				
				# Position data
				$rgene{'loc_start'} = $feat->location->start;
				$rgene{'loc_end'} = $feat->location->end;
				$rgene{'length'} = abs($rgene{'loc_start'} - $rgene{'loc_end'});
				$rgene{'strand'} = $feat->location->strand;
				$rgene{'sequence'} = '';
				# Assumed strand 1 if no strand info
				if (not defined($rgene{'strand'})) {
					print STDERR "Warning: Assuming strand 1 for '$locus' (start at $rgene{'loc_start'})\n";
					$rgene{'strand'} = 1;
				}
				
				# If CDS, also save amino-acid sequence to put in the fasta file
				if ($feat->primary_tag eq 'CDS') {
					$ncds++;
					my $sequence = '';
					# Get the translation, or translate by myself
					if ($feat->has_tag('translation')) {
						$sequence = join(',', $feat->get_tag_values('translation'));
					}
					# Even if we get the sequence, we need to get the GC from the nucleotids
					my $nucl = '';
					if ($seq->length() > 0) {
						# Get nucl sequence
						$nucl = $feat->spliced_seq->seq;
						$rgene{'GC'} = nucleotid_GC($nucl);
						# Add this proportion to the total GC mean
						$genome{'GC'} += $rgene{'GC'} * $rgene{'length'};
						$genome{'prot_lengths'} += $rgene{'length'};
						
						# Translate to AA if needed
						if ($sequence eq '') {
							my $seqobj = Bio::Seq->new(-seq => $nucl, -alphabet => 'dna');
							$sequence = $seqobj->translate->seq;
							$sequence =~ s/\*//;
						}
					}
					
					
					if ($sequence ne '') {
						my $seqobj = Bio::Seq->new( -display_id => $locus, -seq =>  $sequence );
						push @sequences, $seqobj;
						$rgene{'sequence'} = $sequence;
					}
				# Keep "pseudogenes"
				} elsif ($feat->primary_tag eq 'gene' and $feat->has_tag('pseudo')) {
					$rgene{'feat'} = 'pseudo';
					$rgene{'sequence'} = $feat->spliced_seq->seq if $seq->length() > 0;
					# Compute GC
					$rgene{'GC'} = nucleotid_GC($rgene{'sequence'});
				# Take tRNA, rRNA etc.
				} elsif ($feat->primary_tag =~ /RNA/) {
					$rgene{'sequence'} = $feat->spliced_seq->seq if $seq->length() > 0;
					# Compute GC
					$rgene{'GC'} = nucleotid_GC($rgene{'sequence'});
				# But don't take misc_feature and operon
				} else {
					print STDERR "Skip ".$feat->primary_tag." (potential $locus)\n" if $opt{v};
					$num--;
					next FEAT;
				}
				$gene{$rname}{$locus} = \%rgene;
			}
			$gpart{'nbprotein'} = $ncds;
			push @gparts, \%gpart;
		}
		# Finish computing GC for the genome
		$genome{'GC'} = $genome{'prot_lengths'} > 0 ? $genome{'GC'} / $genome{'prot_lengths'} : 0;
		push @genomes, \%genome;
		# And for every gene, compute delta GC
		my $gene_final = delta_GC(\%gene, \%genome);
		
		# Save all CDS for this genome
		$genes{$abbr} = $gene_final;
		
	}
	return \@sequences, \%genes, \@genomes, \@gparts;
}

sub delta_GC
{
	my ($genes, $genome) = @_;
	
	my $base_GC = $genome->{'GC'};
	
	foreach my $gpart (keys %$genes) {
		foreach my $locus (keys %{$genes->{$gpart}}) {
			my $g = $genes->{$gpart}->{$locus};
			$g->{'delta_GC'} = $g->{'GC'} - $base_GC;
			$genes->{$gpart}->{$locus} = $g;
		}
	}
	return $genes;
}

sub nucleotid_GC
{
	my ($seq) = @_;
	my @nucls = split(//, $seq);
	my $GC = 0;
	my $total = 0;
	foreach my $n (@nucls) {
		if ($n =~ /^[GC]$/i) {
			$GC++;
		}
		$total++;
	}
	return $total > 0 ? $GC/$total : 0;
}

# To sort the CDS
sub direc
{
	my ($cds) = @_;
	
	# We use the end of the ORFs, as it is more accurate than the start
	# This is useful for ORFs that are too long and overlap with other ORFs
	# also, we substract X pb to avoid overlapping of orfs in different strands (it happens)
	# This diff can't be shorter the half the length (to avoid getting out)
	my $diff = 300;
	my $midlen = $cds->{'length'}/2;
	$diff = $midlen if $midlen < $diff;
	if ($cds->{'strand'} == 1) {
		return $cds->{'loc_end'} - $diff;
	} else {
		return $cds->{'loc_start'} + $diff;
	}
}

# Order all genes in the genome (assuming linear sequence) and give them numbers
# 3 numbers:
# pnum_all	= number in the whole genome
# pnum_CDS		= number of the CDS
sub order_genes
{
	my ($genes) = @_;
	
	foreach my $sp (sort keys %$genes) {
		my $nglob = 0;
		my $ncds = 0;
		foreach my $gpart (sort keys %{$genes->{$sp}}) {
			my $repdata = $genes->{$sp}->{$gpart};
			foreach my $n (keys %$repdata) {
				$repdata->{$n}->{'ord'} = direc($repdata->{$n});
			}
			
			# Order all genes in this gpart
			foreach my $name (sort { $repdata->{$a}->{'ord'} <=> $repdata->{$b}->{'ord'} } keys %$repdata) {
				$nglob++;
				$repdata->{$name}->{'pnum_CDS'} = -1;
				$repdata->{$name}->{'pnum_all'} = $nglob;
				
				# CDS?
				if ($repdata->{$name}->{'feat'} eq 'CDS') {
					$ncds++;
					$repdata->{$name}->{'pnum_CDS'} = $ncds;
				}
			}
			$genes->{$sp}->{$gpart} = $repdata;
		}
	}
	return $genes;
}

sub print_fasta
{
	my ($seqs, $outpath) = @_;
	return if not $outpath;
	
	# Open the stream
	my $fasta = Bio::SeqIO->new(
			-file	=> ">$outpath",
			-format	=> 'fasta',
	);
	
	foreach my $seq (@$seqs) {
		$fasta->write_seq($seq);
	}
}

# Print the genes data in tabulated form (easy for database import)
sub print_gene_data
{
	my ($data, $outpath, $fields) = @_;
	return if not $outpath;
	
	open(DATA, ">$outpath") or die("$outpath: $!");
	
	print DATA join("\t", @$fields) . "\n";
	
	my $warned = 0;
	
	foreach my $sp (sort keys %$data) {
		foreach my $gpart (sort keys %{$data->{$sp}}) {
			my $repdata = $data->{$sp}->{$gpart};
			
			# Print the sequences in order of location
			foreach my $name (sort { $repdata->{$a}->{'pnum_all'} <=> $repdata->{$b}->{'pnum_all'} } keys %$repdata) {
				my @line = ();
				foreach my $field (@$fields) {
					# Print value if defined
					if (defined($repdata->{$name}->{$field})) {
						push @line, $repdata->{$name}->{$field};
					# Warn if a value is not defined (default: empty string)
					} else {
						if ($warned < 10) {
							warn("No value for field '$field' for data '$name'\n");
						}
						warn("...\n") if $warned == 10;
						$warned++;
						push @line, '';
					}
				}
				print DATA join("\t", @line) . "\n";
			}
		}
	}
	close(DATA);
	print "$warned warnings!\n" if $warned >= 10;
}

# Same as the print_gene_data but more generic: used for genome table
sub print_data
{
	my ($data, $outpath, $fields) = @_;
	return if not $outpath;
	
	open(DATA, ">$outpath") or die("$outpath: $!");
	
	print DATA join("\t", @$fields) . "\n";
	
	my $warned = 0;
	
	foreach my $d (@$data) {
		my @line = ();
		foreach my $field (@$fields) {
			if (defined($d->{$field})) {
				push @line, $d->{$field};
			} else {
				if ($warned < 10) {
					warn("No value for field '$field' for data '$d'\n");
				}
				warn("...\n") if $warned == 10;
				$warned++;
				push @line, '';
			}
		}
		print DATA join("\t", @line) . "\n";
	}
	close(DATA);
	print "$warned warnings!\n" if $warned >= 10;
}

##################################
# MAIN
init();

# One file or several?
my ($sequences, $genes_data, $genomes_data, $gparts_data) = parse_gbk($opt{i}, $opt{p});
$genes_data = order_genes($genes_data);
print_fasta($sequences, $opt{f});
print_gene_data($genes_data, $opt{o}, \@genes_fields);
print_data($genomes_data, $opt{g}, \@genomes_fields);
print_data($gparts_data, $opt{G}, \@gpart_fields);

__END__
