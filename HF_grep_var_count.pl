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
use threads;
use Storable qw ( freeze thaw );
use Thread::Queue;

use Dir::Self;
use lib __DIR__ . '/lib';
use Score;
use Cigar;
use Sample;
use Design;


$| = 1;

#-------------------------------------------------------------------------------------------
#---    CONSTANTS
#-------------------------------------------------------------------------------------------

my $qscore_averaging_range      = 1; # Phred quality score is average in this window/ window length is 1+2*$qscore_averaging_range.
my $qscore_min  = 16; # Ignore positions with base quality score lower then this value;
my $minimum_coverage = 2; # Positions with coverage lower this value will be ignored (defined as non-detectable)
$qscore_min = Score->new($qscore_min);

head();
sub head {
	my $inputBam		= $ARGV[0];
	my $panelFile		= $ARGV[1];
	my $vcfFile		= $ARGV[2];
	my $Design		= Design->new();
	my $Sample		= $Design->newSample($inputBam);
	$Sample->init();
	$Design->init({seqdic => $Sample->header, VCF => $vcfFile, BED => $panelFile});
	$Design->{config}->{qscore_averaging_range} = $qscore_averaging_range;
	$Sample->normalizeBQ();
	#print STDERR Dumper $Sample->{BQMatrix};
	foreach my $seg (@{$Design->segments}) {
		#print STDERR Dumper $seg;
		next if scalar (@{$seg->{variations}}) eq 0;
		#print $seg->{contig},"\t",$seg->{start},"\t",$seg->{end},"\t",scalar (@{$seg->{variations}}),"\n";
		$Sample->pipeline($seg);
		foreach my $CandidateVariation (@{$seg->{variations}}) {
			my $index = $CandidateVariation->{index};
			#next if $CandidateVariation->{position} ne '6529203';
			#print STDERR (scalar @{$Sample->allele($index)->{reads}}),"\n";
			foreach my $amplicon (uniq (map {$_->{amplicon}} @{$Sample->allele($index)->{reads}})) {
				foreach my $strand (qw(-1 1)) {
					#my $refCount = $Sample->allele($index)->readCount({vote => 'ref', strand => $strand, amplicon => $amplicon});
					my $altCount = $Sample->allele($index)->readCount({vote => 'alt', strand => $strand, amplicon => $amplicon});
					my $DP = scalar (grep {($_->{amplicon} eq $amplicon)and($_->{strand} eq $strand)} @{$Sample->allele($index)->{reads}});
					#next if $refCount + $altCount <= 0;
					#print STDERR Dumper $Sample->allele($index);
					print "$index\t$amplicon\t$strand\t$altCount\t$DP\n";
					}
				}
			#my $freq = $altCountSum/($refCountSum + $altCountSum);
			#print "$index\t$freq\n";
			}
		}
	exit;
	}











