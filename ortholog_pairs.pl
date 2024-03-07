#!/usr/bin/perl -w

# Copyright 2023 Matthieu Barba
# This program is free software under AGPLv3 license
# License terms are in the LICENSE file, or at <http://www.gnu.org/licenses/>.

use strict;
use warnings;
use Getopt::Std;
use Bio::SeqIO;
use Data::Dumper;

our %opt;	# Getopt options

my %methods = (
	brh => {
		e_value => 1E-10,
		evalue_tolerance => 1,
		length => 0.40,
		identity => 0.40,
		synteny => 1,
	}
);
my @methods_order = sort keys(%methods);
my $methods_list = '"' . join('", "', @methods_order) . '"';

#################################################
# Message about this program and how to use it
sub usage {
	print STDERR "[ $_[0] ]\n" if $_[0];
	print STDERR << "EOF";
	
	This script creates ortholog pairs between two genomes based on a blasthit result.
	
	-h        : this (help) message
	
	Input
	-i <path> : blasthits input file
	-g <path> : genes data file
	Output
	-o <path> : pairs output file
	
	-m <str>  : select the method (available: $methods_list) (default: "$methods_order[0]")
	-p <str>  : method parameters in the form name1=val,name2=val (default parameter for those missing)
	-s        : don't solve multiple hits with synteny (on by default)
	
	-v        : verbose
EOF
	print STDERR "\nDefault parameters list:\n";
	foreach my $m (@methods_order) {
		print_parameters($m, $methods{$m});
	}
	exit 1;
}

##################################
# Command line options processing
sub init {
	getopts('hi:g:o:p:m:v', \%opt) or usage();
	usage() if $opt{h};
	usage("Blasthits file needed (-i)") unless $opt{i};
	usage("Genes data file needed (-g)") unless $opt{g};
	usage("Pairing file needed (-o)") unless $opt{o};
	# What method will be used?
	if (defined($opt{m}) and not $methods{$opt{m}}) {
		die("Method '$opt{m}' is not defined");
	} elsif (not defined($opt{m})) {
		$opt{m} = $methods_order[0];
		print STDERR ("Default method: $opt{m}\n") if $opt{v};
	}
}

# Read the genes data file
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
		chomp($line);
		my @fields = split(/\t/, $line);
		
		my %data = ();
		for (my $i = 0; $i < @prothead; $i++) {
			$data{ $prothead[$i] } = $fields[$i];
		}
		$protdata{ $data{'pid'} } = \%data;
		$protorder{ $data{'sp'} }{ $data{'pnum_CDS'} } = $data{'pid'};
	}
	close(DATA);
	
	return \%protdata, \%protorder;
}

###############################################################################################
# HITS IMPORT FUNCTIONS

# Keep the blast hit if it is the best
sub add_match {
	my ($matches, $hit, $protdata) = @_;
	
	my $from = $hit->{'query'};
	my $to = $hit->{'subject'};
	
	# First check if both sequences are in different genomes
	if (not defined($protdata->{$from}) or not defined($protdata->{$to})) {
		die "No protein data for one or both of these: $from, $to";
	}
	my $sp1 = $protdata->{$from}->{'sp'} or die("No sp for $from");
	my $sp2 = $protdata->{$to}->{'sp'} or die("No sp for $to");
	
	if ($sp1 eq $sp2) {
		return;
	}
	
	# Is there already a hit for this sequence between the same two genomes?
	my $matches_sps = $matches->{$sp1}->{$sp2};
	if (defined($matches_sps->{$from})) {
		# Compare the current best with the new hit
		my $eval = $matches_sps->{$from}->{'e_value'};
		my $id = $matches_sps->{$from}->{'identity'};
		my $length = $matches_sps->{$from}->{'length'};
		
		# Same evalue (with some tolerance)?
		if ($hit->{'e_value'} == $eval or
			($eval > 0 and $hit->{'e_value'} / $eval <= $methods{brh}{'evalue_tolerance'})) {
			
			# If so, maybe same identity too?
			if ($hit->{'identity'} == $id) {
				push @{ $matches->{$sp1}->{$sp2}->{$from}->{'matches'} }, $to;
			}
			# Or a better identity?
			elsif ($hit->{'identity'} > $id) {
				$matches->{$sp1}->{$sp2}->{$from} = init_match($hit, $to);
			}
		}
		# Or better evalue?
		elsif ($eval > $hit->{'e_value'}) {
				$matches->{$sp1}->{$sp2}->{$from} = init_match($hit, $to);
		}
		
		# Skip otherwise
	}
	# First time for this: init
	else {
		$matches->{$sp1}->{$sp2}->{$from} = init_match($hit, $to);
	}
}

# Create a hit "object" with its first match
# Return this object
sub init_match {
	my ($bhit, $to) = @_;
	
	my $match = {
		'e_value'	=> $bhit->{'e_value'},
		'identity'	=> $bhit->{'identity'},
		'length' => $bhit->{'length'}
	};
	@{ $match->{'matches'} } = ($to);
	
	return $match;
}

# Filter the hits before further analysis
sub filter_hit {
	my ($hit, $protdata, $pars) = @_;
	
	my $query = $hit->{'query'};
	my $subj = $hit->{'subject'};
	
	# First, check that every sequence has data
	if (not defined($protdata->{ $query })) {
		die("No protein data for '$query'");
	}
	if (not defined($protdata->{ $subj })) {
		die("No protein data for '$subj'");
	}
	
	# Don't care about same genome comparisons
	if ( $protdata->{ $query }->{ 'sp' } eq $protdata->{ $subj }->{ 'sp' } ) {
		return 0;
	}
	
	my $qlen = $protdata->{ $query }->{'length'} / 3;	# /3: nucleic length!
	my $slen = $protdata->{ $subj  }->{'length'} / 3;
	my $smaller = $qlen < $slen ? $qlen : $slen;
	
	# Filter by length, identity and evalue
	if ($hit->{'length'} < $smaller * $pars->{'length'} or
			$hit->{'identity'} < $pars->{'identity'} or
			$hit->{'e_value'}+0 > $pars->{'e_value'}) {
		return 0;
	} else {
		return 1;
	}
}

# Create ortholog pairs based on BRH
# print directly the result to the file, returns nothing
sub method_brh {
	my ($inpath, $outpath, $protdata, $protorder, $pars) = @_;
	
	my %matches = ();
	my $nhits = 0;
	my $nhits_all = 0;
	
	# Read the blasthits
	open(BH, $inpath) or die("$inpath: $!");
	while(my $line = <BH>) {
		next if $line =~ /^#/;
		chomp($line);
		my %hit = ();
		(
			$hit{'query'},
			$hit{'subject'},
			$hit{'identity'},
			$hit{'length'},
			$hit{'mismatches'},
			$hit{'gap_openings'},
			$hit{'query_start'},
			$hit{'query_end'},
			$hit{'subject_start'},
			$hit{'subject_end'},
			$hit{'e_value'},
			$hit{'bit_score'},
		) = split(/\t/, $line);
		
		$nhits_all++;
		if (filter_hit(\%hit, $protdata, $pars)) {
			add_match(\%matches, \%hit, $protdata);
			$nhits++;
		}
	}
	close(BH);
	warn("$nhits / $nhits_all filtered hits\n") if $opt{v};
	
	# From the matches, compute all the best pairs
	my @pairs = make_all_brh_pairs(\%matches, $protdata, $protorder, $pars);
	print_pairs($outpath, \@pairs);
}

###############################################################################################
# BRH FUNCTIONS

sub add_to_group {
	my ($best, $recip, $grp_counter, $from_grp, $to_grp) = @_;
	
	# Add them to a group (create the group if no group found)
	my $grp = -1;
	# Search group from the query
	foreach my $q (@$recip) {
		if (defined($from_grp->{$q})) {
			$grp = $from_grp->{$q};
			last;
		}
	}
	# Not found? Then search group from the pairs
	if ($grp == -1) {
		foreach my $p (@$best) {
			if (defined($to_grp->{$p})) {
				$grp = $to_grp->{$p};
				last;
			}
		}
	}
	
	# No group found? Create a new one
	if ($grp == -1) {
		$grp_counter++;
		$grp = $grp_counter;
	}
	
	# Add all elements in this group
	foreach my $q (@$recip) {
		$from_grp->{$q} = $grp;
	}
	foreach my $p (@$best) {
		$to_grp->{$p} = $grp;
	}
	
	return($grp_counter, $from_grp, $to_grp);
}

sub make_all_brh_pairs {
	my ($matches, $protdata, $protorder, $pars) = @_;
	
	my @pairs = ();
	foreach my $sp1 (sort keys %$matches) {
		foreach my $sp2 (sort keys %{ $matches->{$sp1} }) {
			push @pairs, make_pairs($matches, $protdata, $protorder, $sp1, $sp2, $pars);
		}
	}
	
	my $npairs = @pairs;
	warn("SUMMARY\tFound pairs:\t$npairs\n") if $opt{v};
	
	return @pairs;
}

sub make_pairs {
	my ($matches, $protdata, $protorder, $sp1, $sp2, $pars) = @_;
	
	my %multis = ();
	my @solved_pairs = ();
	
	my ($from_grp, $to_grp, $grp_counter, $orthos_from, $orthos_to, $solved) = make_groups($matches, $sp1, $sp2);
	push @solved_pairs, @$solved;
	
	# Now solve the multiple groups and try to get pairs out of them
	if ($pars->{'synteny'}) {
		push @solved_pairs, solve_multiples($from_grp, $to_grp, $grp_counter, $protdata, $protorder, $orthos_from, $orthos_to, $sp1, $sp2);
		# Stats
		my $nfrom = keys %$from_grp;
		my $nto = keys %$to_grp;
		my $nsolved = @solved_pairs;
		my @vals = (
			sprintf("%-12s", $sp1),
			sprintf("%-12s", $sp2),
			sprintf("group_from: %-5s", $nfrom),
			sprintf("group_to: %-5s", $nto),
			sprintf("matches: %-5s", $nsolved),
		);
		print STDERR join("\t", @vals) . "\n" if $opt{v};
	}
	return @solved_pairs;
}

sub make_groups {
	my ($matches, $sp1, $sp2) = @_;	
	my $hits = $matches->{$sp1}->{$sp2};
	my $rhits = $matches->{$sp2}->{$sp1};
	
	my $from_grp = {};
	my $to_grp = {};
	my $grp_counter = -1;
	my @solved_pairs = ();
	my $npairs = 0;
	my %orthos_from = ();
	my %orthos_to = ();
	
	# For this pair of genomes, find all best reciprocal matches
	foreach my $query (keys %$hits) {
		# One or several best hits?
		my @best = @{ $hits->{$query}->{'matches'} };
		
		if (@best == 0) {
			print STDERR "NOMATCH_NOPAIR\t$query is not paired to anything\n" if $opt{v};
		}
		
		# Only one best? That's easier
		elsif (@best == 1) {
			my $pair = $best[0];
			
			# Get all reciprocal hits if any
			if ($rhits->{$pair}) {
				my @recip = @{ $rhits->{$pair}->{'matches'} };
				
				# Is there one and only one best reciprocal hit?
				if (@recip == 1) {
					# But is it the one we want (the query)?
					if ($recip[0] eq $query) {
						# FOUND! Keep this pair
						print STDERR "PAIR $query -> $pair\n" if $opt{v};
						my %pairdata = (
							'pid1' => $query,
							'pid2' => $pair,
							'o_ident' => $rhits->{$pair}->{'identity'},
							'o_alen' => $rhits->{$pair}->{'length'},
						);
						$npairs++;
						push @solved_pairs, \%pairdata;
						# Keep the info for the more controversial choices
						$orthos_from{$query} = $pair;
						$orthos_to{$pair} = $query;
					# Not the best: give up on the query
					} else {
						if ($opt{v}) {
							print STDERR "NOPAIR\t$query is not paired with $pair (the best was $recip[0])\n" if $opt{v};
						}
					}
				# There are several matches back!
				} else {
					# Let's check that our query is at least among them
					my $rematch = 0;
					foreach my $back_query (@recip) {
						if ($back_query eq $query) {
							$rematch = 1;
							last;
						}
					}
					
					# We did find a match back to our query, damn
					if ($rematch == 1) {
						# This is a multiple group
						print STDERR "MULTI_GROUP_BACK\t$query - $pair <- is multiply paired back (".join(',', @recip).")\n" if $opt{v};
						my @recipq = (@recip, $query);
						($grp_counter, $from_grp, $to_grp) = add_to_group(\@best, \@recipq, $grp_counter, $from_grp, $to_grp);
					}
				}
			} else {
				# No match back...
				print STDERR "NOMATCH\tPair $pair does not best match back $query.\n" if $opt{v};
			}
		} else {
			# Check if there is any match back at all
			my @match_back = ();
			foreach my $pair (@best) {
				if ($rhits->{$pair}) {
					my @recip = @{ $rhits->{$pair}->{'matches'} };
					foreach my $rec (@recip) {
						push @match_back, $pair if $rec eq $query;
					}
				}
			}
			
			# If there is no match back to the query, no need to go further
			if (@match_back == 0) {
				print STDERR "MULTI_NOPAIR\t$query has multiple pairs (".join(', ', @best).") but none match back\n" if $opt{v};
			}
			elsif (@match_back == 1) {
				print STDERR "MULTI_PAIR\t$query - $match_back[0] (the only one that matched back)\n" if $opt{v};
				my @recipq = ($query);
				($grp_counter, $from_grp, $to_grp) = add_to_group(\@best, \@recipq, $grp_counter, $from_grp, $to_grp);
			}
			else {
				print STDERR "MULTI_GROUP\t$query has multiple pairs matching back\n" if $opt{v};
				my @recipq = ($query);
				($grp_counter, $from_grp, $to_grp) = add_to_group(\@best, \@recipq, $grp_counter, $from_grp, $to_grp);
			}
		}
	}
	return $from_grp, $to_grp, $grp_counter, \%orthos_from, \%orthos_to, \@solved_pairs;
}

sub solve_multiples {
	my ($from_grp, $to_grp, $grp_counter, $protdata, $protorder, $orthos_from, $orthos_to, $sp1, $sp2) = @_;
	
	my @solved_pairs = ();
	my $solving = 1;
	while($solving) {
		$solving = 0;
		my @new_solved_pairs = ();
	
		# Restore the groups
		my @groups = ();
		my ($nfrom, $nto) = (0,0);
		
		foreach my $f (keys %$from_grp) {
			my $grp = $from_grp->{$f};
			# Check if it has not been put in a pair already, just to be sure
			if (not $orthos_from->{$f}) {
				push @{$groups[$grp]->{'from'}}, $f;
				$nfrom++;
			} else {
				print STDERR "FROM\t$f is already in a pair!\n" if $opt{v};
			}
		}
		foreach my $t (keys %$to_grp) {
			my $grp = $to_grp->{$t};
			# Check if it has not been put in a pair already, just to be sure
			if (not $orthos_to->{$t}) {
				push @{$groups[$grp]->{'to'}}, $t;
				$nto++;
			} else {
				print STDERR "TO\t$t is already in a pair!\n" if $opt{v};
			}
		}
		print STDERR "Groups: from $nfrom\tto $nto\n" if $opt{v};
		
		# Solve each group
		GROUP : for (my $g = 0; $g <= $grp_counter; $g++) {
			# Check that this group has both From and To elements
			if (not (defined($groups[$g]->{'from'}) and defined($groups[$g]->{'to'}))) {
				next GROUP;
			}
			
			my @froms = @{ $groups[$g]->{from} };
			my @tos = @{ $groups[$g]->{to} };
			
			my $fromt = join(', ', sort @froms);
			my $tot = join(', ', sort @tos);
			print STDERR "$g\t$fromt   ->   $tot\n" if $opt{v};
			
			# Is there anything left?
			if ($#froms+1 == 0) {
				print STDERR "No 'from' left for ".join(', ', @tos)."\n" if $opt{v};
			}
			elsif ($#tos+1 == 0) {
				print STDERR "No 'to' left for ".join(', ', @froms)."\n" if $opt{v};
			}
			# Is there maybe just one of each left?
			elsif (@froms == 1 and @tos == 1) {
				my $q = $froms[0];
				my $p = $tos[0];
				my %pair = (
					'pid1' => $q,
					'pid2' => $p,
					'o_ident' => 0,
					'o_alen' => 0,
				);
				push @solved_pairs, \%pair;
				$orthos_from->{$q} = $p;
				$orthos_to->{$p} = $q;
				print STDERR "LEFTOVER_MATCH\t$q - $p\n" if $opt{v};
			}
			# Otherwise: lots of work to find the matches with synteny
			else {
				my @synt_solved = solve_synteny(\@froms, \@tos, $protdata, $protorder, $orthos_from, $orthos_to, $sp1, $sp2);
				my $nsynt = @synt_solved;
				print STDERR "$nsynt Synteny pairs found\n" if $opt{v};
				push @new_solved_pairs, @synt_solved;
			}
		}
		my $nsolved = $#solved_pairs+1;
		push @solved_pairs, @new_solved_pairs;
		my $nchanges = $#new_solved_pairs+1;
		if ($nchanges > 1) {
			print STDERR "SYNTENY_SEARCH: $sp1 vs $sp2\tanother round\n" if $opt{v};
			$solving = 1;
		} else {
			$solving = 0;
		}
	}
	
	return @solved_pairs;
}

sub solve_synteny {
	my ($froms, $tos, $protdata, $protorder, $orthos_from, $orthos_to, $sp1, $sp2) = @_;
	
	my @potential_pairs = ();
	
	FROM : foreach my $f (@$froms) {
		my $fnum = $protdata->{$f}->{'pnum_CDS'};
		my $fbefore = $protorder->{$sp1}->{$fnum-1};
		my $fafter  = $protorder->{$sp1}->{$fnum+1};
		my $fbefore_ortho = ($fbefore and $orthos_from->{$fbefore}) ? $orthos_from->{$fbefore} : '';
		my $fafter_ortho  = ($fafter  and $orthos_from->{$fafter})  ? $orthos_from->{$fafter}  : '';
		
		if ($fbefore_ortho or $fafter_ortho) {
			foreach my $t (@$tos) {
				my $tnum = $protdata->{$t}->{'pnum_CDS'};
				my $tbefore = $protorder->{$sp2}->{$tnum-1};
				$tbefore = '' if not $tbefore;
				my $tafter  = $protorder->{$sp2}->{$tnum+1};
				$tafter = '' if not $tafter;
				
				if
				   ($tbefore eq $fbefore_ortho
				or $tbefore eq $fafter_ortho
				or $tafter eq $fbefore_ortho
				or $tafter eq $fafter_ortho) {
					my %pair = (
						'pid1' => $f,
						'pid2' => $t,
						'o_ident' => 0,
						'o_alen' => 0,
					);
					push @potential_pairs, \%pair;
				}
			}
		}
	}
	
	# Check potential pairs and only keep those that only use each CDS one and only one time
	my %from_done = ();
	my %to_done   = ();
	my @synt_solved = ();
	
	my $n_pot = @potential_pairs;
	print STDERR "$n_pot potential matches\n" if $opt{v};
	foreach my $pair (@potential_pairs) {
		push @{ $from_done{ $pair->{'pid1'} } }, $pair;
		push @{   $to_done{ $pair->{'pid2'} } }, $pair;
	}
	
	# Check FROM
	foreach my $pidfrom (keys %from_done) {
		my @list = @{ $from_done{$pidfrom} };
		if (@list > 1) {
			print STDERR "Several pairs possible from $pidfrom: can't keep them\n" if $opt{v};
			delete $from_done{$pidfrom};
			foreach my $pair (@list) {
				delete $to_done{ $pair->{'pid2'} };
			}
		}
	}
	
	# Check TO
	foreach my $pidto (keys %to_done) {
		my @list = @{ $to_done{$pidto} };
		if (@list > 1) {
			print STDERR "Several pairs possible to $pidto: can't keep them\n" if $opt{v};
			delete $to_done{$pidto};
			foreach my $pair (@list) {
				delete $from_done{ $pair->{'pid1'} };
			}
		}
	}
	
	# Finally, add the correct pairs
	foreach my $pid (keys %from_done) {
		my $pair = ${ from_done{$pid} }->[0];
		push @synt_solved, $pair;
		my $q = $pair->{'pid1'};
		my $p = $pair->{'pid2'};
		next if not $to_done{$p};
		$orthos_from->{$q} = $p;
		$orthos_to->{$p} = $q;
		print STDERR "SYNTENY_MATCH\t$q - $p\n" if $opt{v};
	}
	
	return @synt_solved;
}

sub print_pairs {
	my ($outpath, $pairs) = @_;

	my @ortho_headers = qw(oid	pid1	pid2	o_ident	o_alen);
	open(OUT, ">$outpath") or die("$outpath: $!");
	print OUT join("\t", @ortho_headers) . "\n";
	my $n = 0;
	foreach my $pair (sort { $a->{'pid1'} cmp $b->{'pid1'} } @$pairs) {
		$n++;
		$pair->{'oid'} = $n;
		my @line = ();
		foreach my $h (@ortho_headers) {
			push @line, $pair->{$h} ? $pair->{$h} : '';
		}
		print OUT join("\t", @line) . "\n";
	}
	close(OUT);
}

sub change_method_parameters {
	my ($pars, $changes) = @_;
	
	# Copy default values
	my %new_pars = %$pars;
	
	# Changes custom values
	my @assigns = split(/,/, $changes);
	foreach my $as (@assigns) {
		if ($as =~ /^([^=]+)=([^=]+)$/) {
			$new_pars{$1} = $2;
		} else {
			warn("Wrong parameter assignation (missing '=' ?)");
		}
	}
	
	return \%new_pars;
}

sub print_parameters {
	my ($name, $pars) = @_;
	return if not defined($name) or not defined($methods{$name});
	return if not $opt{v};
	print STDERR "Method: $name\n";
	foreach my $k (sort keys %$pars) {
		print STDERR "\t$k\t$pars->{$k}\n";
	}
}

##################################
# MAIN
init();

# Prepare the genes data before analyzing the blast
my ($protdata, $protorder) = get_genes_data($opt{g});
my $method_name = $opt{m};
my $method_pars = $methods{ $method_name };
if (defined($method_pars)) {
	my $pars = $opt{p} ? change_method_parameters($method_pars, $opt{p}) : $method_pars;
	print_parameters($method_name, $pars);
	# BRH
	if ($opt{m} eq 'brh') {
		method_brh($opt{i}, $opt{o}, $protdata, $protorder, $pars);
	}
}

__END__
