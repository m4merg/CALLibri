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
use Design;


$| = 1;

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

sub head {
	my $inputBam		= $ARGV[0];
	my $panelFile		= $ARGV[1];
	my $vcfFile		= $ARGV[2];
	my $controlBamList	= $ARGV[3];
	my $Design		= Design->new();
	my $Sample		= $Design->Sample($inputBam);

	$Design->init({seqdic => $Sample->header, VCF => $vcfFile, BED => $panelFile});
	$Design->{config}->{qscore_averaging_range} = $qscore_averaging_range;
	
	die "Can not identify genome regions in input BED file ($panelFile)\n" unless defined $Design->segments;
	
	open (READ, "<$controlBamList");

	while (<READ>) {
		chomp;
		next if m!^#!;
		$Design->add_control($_);
		}

	close READ;

	foreach my $seg (@{$Design->segments}) {
		next if scalar (@{$seg->{variations}}) eq 0;
		print $seg->{contig},"\t",$seg->{start},"\t",$seg->{end},"\t",scalar (@{$seg->{variations}}),"\n";
		$Sample->pipeline($seg);
		my $n = 0;
		next;
		foreach my $Control (@{$Design->controls}) {
			print STDERR "$n\n";
			$Control->pipeline($seg);
			++$n;
			}
		foreach my $CandidateVariation (@{$seg->{variations}}) {
			my $index = $CandidateVariation->{index};
			next if $CandidateVariation->{position} ne '6529203';
			print STDERR (scalar @{$Sample->allele($index)->{reads}}),"\n";
			foreach my $amplicon (uniq (map {$_->{amplicon}} @{$Sample->allele($index)->{reads}})) {
				foreach my $strand (qw(-1 1)) {
					#print STDERR "$amplicon\t$strand\tTOTAL\t",(scalar (grep {$_->{strand} eq $strand && $_->{amplicon} eq $amplicon} @{$Sample->allele($index)->{reads}})),"\n";
					#print STDERR "$amplicon\t$strand\tREF\t",(scalar (grep {$_->{strand} eq $strand && $_->{amplicon} eq $amplicon && $_->{vote} eq 'ref'} @{$Sample->allele($index)->{reads}})),"\n";
					#print STDERR "$amplicon\t$strand\tALT\t",(scalar (grep {$_->{strand} eq $strand && $_->{amplicon} eq $amplicon && $_->{vote} eq 'alt'} @{$Sample->allele($index)->{reads}})),"\n";
					print STDERR "$amplicon\t$strand\t",$Sample->allele($index)->readCount({vote => 'ref', strand => $strand, amplicon => $amplicon, BQ => Score->new(16)->prob}),"\n";
					print STDERR "$amplicon\t$strand\t",$Sample->allele($index)->readCount({vote => 'alt', strand => $strand, amplicon => $amplicon, BQ => Score->new(16)->prob}),"\n";
					}
				}
			}
		last;
		}
	exit;
	}






















