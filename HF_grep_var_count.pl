#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Sam;		#Install
use Bio::Cigar;			#Install
use Try::Tiny;			#Install
use Data::Dumper;
use List::Util qw(sum max min shuffle pairs);
use Switch;			#Install
use Storable 'dclone';
use Cwd qw(getcwd abs_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage;
use File::Basename;
use Time::HiRes qw(gettimeofday);
use List::MoreUtils qw(uniq);
use Storable 'dclone';

use Dir::Self;
use lib __DIR__ . '/lib';
use Score;
use Cigar;
use Sample;


$| = 1;

#	{
#	"contig" : 
#	"position" : 
#	"ref" : 
#	"alt" :
#	"name" : 
#	"reads" : 
#		[
#			{
#			"name" : 
#			"BQ" : 
#			"strand" : 
#			"vote" : "ref/alt"
#			"amplicon" : 
#			}
#		]
#	}

#-------------------------------------------------------------------------------------------
#---    CONSTANTS
#-------------------------------------------------------------------------------------------

my $qscore_averaging_range      = 1; # Phred quality score is average in this window/ This value defines half of the window length.
my $qscore_min  = 16; # Ignore positions with base quality score lower then this value;
my $minimum_coverage = 2; # Positions with coverage lower this value will be ignored (defined as non-detectable)
$qscore_min = Score->new($qscore_min);

head();

#-------------------------------------------------------------------------------------------
#---    CORE SUBROUTINES              
#-------------------------------------------------------------------------------------------

sub get_stat { # see pipeline function
	my $Mutation = shift;
	my $alignment = shift;
	#aps,ape,rps,rpe,qps,qpe - 0-based coordinates
	#aps/ape - start/end position of variant site in alignment
	#rps/rpe - start/end positions of variation site in reference
	my $qps; my $rps; my $aps; my $qpe; my $rpe; my $ape;
	my $opqs; my $opqe;
	$rps = ($Mutation->{position} - $alignment->start);
	$rpe = ($Mutation->{position} - $alignment->start) + length($Mutation->{ref});
	return undef if $rps < 1;
	return undef if ($rps) > ($rpe);

	my ($ref, $match, $query) = $alignment->padded_alignment;

	my $ts = $ref; $ts =~ s/-//g;
	return undef if ($rpe) > length($ts);

	my $cigar = Cigar->new($alignment->cigar_str);
	try {($qps, $opqs) = $cigar->rpos_to_qpos($rps);};
	return undef unless defined $qps;
	$qpe = $qps + length($Mutation->{alt});
	return undef if $qpe < 0;
	return undef if $qps < 0;
	$aps = $rps + $cigar->get_shift($rps);
	$ape = $aps + max(length($Mutation->{alt}),length($Mutation->{ref}));
	my $stat;
	$stat->{oref} = substr($ref, $aps, $ape-$aps);
	$stat->{oalt} = substr($query, $aps, $ape-$aps);
	$stat->{oref} =~ s/-//g;
	$stat->{oalt} =~ s/-//g;

	$stat->{match} = $match;
	$stat->{aps} = $aps; $stat->{ape} = $ape;
	$stat->{qps} = $qps; $stat->{qpe} = $qpe;

	return $stat;
	}

sub get_qscore { # see pipeline function
	my $alignment	= shift;
	my $stat	= shift;
	my @scores = @{$alignment->qscore};
	my $qscore_start = max(0, $stat->{qps} - $qscore_averaging_range);
	my $qscore_end   = min($stat->{qpe} - 1 + $qscore_averaging_range, (scalar @scores) - 1);
	my $qscore = 0;
	map {$qscore = $qscore + Score->new($_)->prob} @scores[($qscore_start)..($qscore_end)];
	$qscore = $qscore / ($qscore_end - $qscore_start + 1);
	return $qscore;
	}

sub pipeline {
	my $sam			= shift;
	my $segment		= shift;
	my $ampliconsHash	= shift;
	
	my $sam_segment = $sam->segment($segment->{contig}, $segment->{start}, $segment->{end});
	return undef unless defined $sam_segment;
	my @all_alignments = $sam_segment->features;
	foreach my $alignment (@all_alignments) {
		foreach my $Mutation (@{$segment->{mutations}}) {
			if (defined($alignment->get_tag_values("SUPPLEMENTARY"))) {
				next if $alignment->get_tag_values("SUPPLEMENTARY") eq '1';
				}
			next if $alignment->get_tag_values("UNMAPPED") eq '1';
			next if $alignment->get_tag_values("NOT_PRIMARY") eq '1';
			
			my $stat = get_stat($Mutation, $alignment);
			next unless defined($stat);
			my $qscore = get_qscore($alignment, $stat);
			#print "!",$alignment->qname,"\t",$stat->{oref},"\t",$stat->{oalt},"\t",$qscore,"\n";
			next if $qscore >= $qscore_min->prob;
			my $read;
			$read->{name}		= $alignment->qname;
			$read->{BQ}		= $qscore;
			$read->{strand}		= $alignment->strand;
			$read->{amplicon}	= select_amplicon($Mutation, $alignment, $ampliconsHash);
			if (($stat->{oref} eq ($Mutation->{ref})) and ($stat->{oalt} eq ($Mutation->{alt}))) {
					$read->{vote} = 'alt';
					push(@{$Mutation->{reads}}, $read);
				} elsif (((length($stat->{oref}) eq length($Mutation->{ref}))and
					(length($stat->{oalt}) eq length($Mutation->{ref})))or
					(substr($stat->{match}, $stat->{aps}, $stat->{ape}-$stat->{aps}) =~ /^\|*$/)) {
					$read->{vote} = 'ref';
					push(@{$Mutation->{reads}}, $read);
					} else {
						#print STDERR "WHAT?\n";
					}
			}
		}
	#return $mutation;
	}

sub select_amplicon {
	my $Mutation		= shift;
	my $alignment		= shift;
	my $ampliconsHash	= shift;
	my $amplicons = $ampliconsHash->{$Mutation->{contig}.":".$Mutation->{position}};
	if (scalar @{$amplicons} eq 0) {
		return undef;
		} elsif ((scalar @{$amplicons}) eq 1) {
		return $amplicons->[0];
		}
	my $min = 999999999999;
	my $selection;
	foreach my $Amplicon (@{$amplicons}) {
		my $start;
		my $end;
		if ($Amplicon =~ /(\d+)-(\d+)/) {
			$start = $1;
			$end = $2;
			}
		$start	= abs($alignment->start - $start);
		$end	= abs($alignment->end - $end);
		if ($start + $end < $min) {
			$min = $start + $end;
			$selection = $Amplicon;
			}
		}
	return $selection;
	}

sub define_amplicons {
	my $inputFile	= shift;
	my $header	= shift;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";
	
	my $amplicon;
	while (<$panel_fh>) {
		my @mas = split/\t/;
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in BED file. Non-concordant with input BAM file\n" unless defined $header->{$mas[0]};
		die "Segment end position ($mas[2]) out of contig length for chromosome $mas[0]\n" if $mas[2] > $header->{$mas[0]};
		for (my $i = $mas[1] + 1; $i <= $mas[2]; $i++) {
			my $name = "$mas[0]:$i";
			$amplicon->{$name} = [] unless defined($amplicon->{$name});
			push @{$amplicon->{$name}}, "$mas[0]:$mas[1]-$mas[2]";
			}
		}

	close $panel_fh;
	return $amplicon;
	}

sub define_segments { # simple merge BED file
	my $inputFile	= shift;
	my $header	= shift;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";

	my $segments;
	my $current_segment;
	while (<$panel_fh>) {
		my @mas = split/\t/;
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in BED file. Non-concordant with input BAM file\n" unless defined $header->{$mas[0]};
		die "Segment end position ($mas[2]) out of contig length for chromosome $mas[0]\n" if $mas[2] > $header->{$mas[0]};
		unless (defined $current_segment) {
			$current_segment = {contig => $mas[0], start => $mas[1], end => $mas[2], mutations => []};
			next;
			}
		if (($mas[0] eq $current_segment->{contig}) and ($mas[1] < $current_segment->{start})) {
			die "Input BED file is not sorted at position $mas[0]:$mas[1] <-; Use \"sort -k1,1 -k2,2n <input_bed>\" to sort bed file\n";
			}
		if (($mas[0] eq $current_segment->{contig})
				and ($mas[1] > $current_segment->{start})
				and ($mas[1] < $current_segment->{end})) {
			if ($mas[2] > $current_segment->{end}) {
				$current_segment->{end} = $mas[2];
				} else {
				next;
				}
			} else {
			push (@{$segments}, $current_segment);
			$current_segment = {contig => $mas[0], start => $mas[1], end => $mas[2], mutations => []};
			}
		}
	return $segments;
	}

sub make_mutation_connection {
	my $mutation;
	
	my $mutation_hash_start;
	foreach my $Mutation (@{$mutation}) {
		if (not(defined($mutation_hash_start->{$Mutation->{contig}}))) {
			$mutation_hash_start->{$Mutation->{contig}} = \$Mutation;
			}
		my $min = 9999999999999;
		my $next;
		foreach my $MutationNext (grep {$_->{contig} eq $Mutation->{contig}} @{$mutation}) {
			my $equal = 1;
			foreach my $arg (qw(position ref alt)) {
				if ($MutationNext->{$arg} ne $Mutation->{$arg}) {$equal = 0}
				}
			next if $equal eq 1;
			next if defined($MutationNext->{prev});
			next if $MutationNext->{position} - $Mutation->{position} < 0;
			if ($MutationNext->{position} - $Mutation->{position} < $min) {
				$min = $MutationNext->{position} - $Mutation->{position};
				$next = \$MutationNext;
				}
			}
		unless (defined $next) {
			$Mutation->{next} = undef;
			} else {
			$Mutation->{next} = $next;
			$$next->{prev} = \$Mutation;
			}
		}
	return $mutation;
	}

sub load_mutations {
	my $vcf = shift;
	my $header = shift;

	open (my $vcf_fh, "<$vcf");

	my $mutations = [];
	while (<$vcf_fh>) {
		chomp;
		next if m!#!;
		my @mas = split/\t/;
		next unless $mas[1] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in VCF file. Non-concordant with input BAM file\n" unless defined $header->{$mas[0]};
		die "Mutation position ($mas[1]) out of contig length for chromosome $mas[0]\n" if $mas[1] > $header->{$mas[0]};

		my $current = {};
		$current->{contig} = $mas[0];
		$current->{position} = $mas[1];
		$current->{ref} = $mas[3];
		$current->{alt} = $mas[4];
		$current->{reads} = [];
		$current->{name} = '.';
		push (@{$mutations}, $current);
		}
	close $vcf_fh;

	if (defined($mutations)) {
		$mutations = [sort {$a->{position} <=> $b->{position}} @{$mutations}];
		}

	return $mutations;
	}	

sub map_mutations_to_segments { # For each segment define which mutations fall within this segment and fill $segment->{mutations} arrays with the corresponding mutation links
	my $segments	= shift;
	my $mutations	= shift;
	foreach my $seg (@{$segments}) {
		foreach my $Mutation (grep {$_->{contig} eq $seg->{contig}} @{$mutations}) {
			if (($Mutation->{position} >= $seg->{start})and($Mutation->{position} <= $seg->{end})) {
				push(@{$seg->{mutations}}, $Mutation);
				}
			}
		}
	return $segments;
	}

sub readCountCalc {
	my $reads = shift;
	my $type  = shift; # Either 'positive' or 'negative'
	$type = 'positive' unless defined $type;
	my $sum;
	foreach my $Read (@{$reads}) {
		if ($type eq 'positive') {
			$sum += 1 - ($Read->{BQ});
			}
		if ($type eq 'negative') {
			$sum += $Read->{BQ}/3;
			}
		}
	return $sum;
	}

sub readCount {
	my $Mutation	= shift;
	my $info	= shift;
	my $count = 0;
	my $reads = dclone $Mutation->{reads};
	my $opposite_reads = [];
	if (defined($info->{strand})) {
		$reads = [grep {$_->{strand} eq $info->{strand}} @{$reads}];
		}
	if (defined($info->{amplicon})) {
		$reads = [grep {$_->{amplicon} eq $info->{amplicon}} @{$reads}];
		}
	if (defined($info->{vote})) {
		$opposite_reads = [grep {$_->{vote} ne $info->{vote}} @{$reads}];
		$reads = [grep {$_->{vote} eq $info->{vote}} @{$reads}];
		}
	if (scalar @{$reads} > 0) {
		$count += readCountCalc($reads, 'positive');
		#print STDERR "->",readCountCalc($reads, 'positive'),"\n";
		}
	if (scalar(@{$opposite_reads}) > 0) {
		$count += readCountCalc($opposite_reads, 'negative');
		#print STDERR "->",readCountCalc($opposite_reads, 'negative'),"\n";
		}
	return $count;
	}

#head($inputBam, $refFile, $refVCF, $outFile, $outVCF, $version, $cmdline);
sub head {
	my $inputBam	= $ARGV[0];
	my $panelFile	= $ARGV[1];
	my $vcfFile	= $ARGV[2];
	my $Sample	= Sample->new($inputBam);
	

	my $segments		= define_segments($panelFile, $Sample->header);
	my $ampliconsHash	= define_amplicons($panelFile, $Sample->header);
	my $mutations		= load_mutations($vcfFile, $Sample->header);
	die "Can not identify genome regions in input BED file ($panelFile)\n" unless defined $segments;
	die "Can not identify genome regions in input BED file ($panelFile)\n" unless defined $ampliconsHash;
	die "Can not identify mutations in input VCF file ($vcfFile)\n" unless defined $mutations;

	map_mutations_to_segments($segments, $mutations);
	

	foreach my $seg (@{$segments}) {
		next if scalar (@{$seg->{mutations}}) eq 0;
		print $seg->{contig},"\t",$seg->{start},"\t",$seg->{end},"\t",scalar (@{$seg->{mutations}}),"\n";
		pipeline($Sample->sam, $seg, $ampliconsHash);
		last;
		}
	
	foreach my $Mutation (@{$mutations}) {
		next if $Mutation->{position} ne '6529203';
		foreach my $amplicon (uniq (map {$_->{amplicon}} @{$Mutation->{reads}})) {
			#print STDERR (scalar @{$Mutation->{reads}}),"\n";
			#print STDERR (scalar (grep {$_->{strand} eq '1'} @{$Mutation->{reads}})),"\n";
			#print STDERR (scalar (grep {$_->{strand} eq '1' && $_->{vote} eq 'ref'} @{$Mutation->{reads}})),"\n";
			#print STDERR (scalar (grep {$_->{strand} eq '1' && $_->{vote} eq 'alt'} @{$Mutation->{reads}})),"\n";
			print STDERR readCount($Mutation, {vote => 'ref', strand => '1', amplicon => $amplicon}),"\n";
			print STDERR readCount($Mutation, {vote => 'alt', strand => '1', amplicon => $amplicon}),"\n";
			}
		}

	exit;
	
	foreach my $Mutation (@{$mutations}) {
		next if $Mutation->{contig} ne 'chr22';
		print STDERR Dumper $Mutation;
		foreach my $Read (@{$Mutation->{reads}}) {
			#next
			#next if $Read->{vote} eq 'ref';
			#print Dumper $Read;
			}
		}
	exit;
	}






















