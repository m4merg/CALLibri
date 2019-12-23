#!/usr/bin/env perl

use strict;
use warnings;
use Dir::Self;
use lib __DIR__ . '/lib';
use Sample;
use Design;
use Data::Dumper;
use threads;
use Storable qw ( freeze thaw );
use Thread::Queue;

my $qscore_averaging_range      = 1; # Phred quality score is average in this window/ This value defines half of the window length.
my $minimum_coverage = 2; # Positions with coverage lower this value will be ignored (defined as non-detectable)

my $work_q   = Thread::Queue->new;

sub worker {
	while ( my $passed = $work_q->dequeue ) {
		my $bam		= $passed->[0];
		my $panel	= $passed->[1];
		my $vcf		= $passed->[2];
		my $seed	= $passed->[3];
		print STDERR "Started $bam\n";
		my $cmd = "perl HF_grep_var_count.pl $bam $panel $vcf > temp/$seed";
		print STDERR "$cmd\n";
		`$cmd`;
		}
	}

threads->create( \&worker ) for 1 .. 13;

open (READ, "<bamListRef");

my $n = 0;
while (<READ>) {
	chomp;
	next if m!^#!;
	my $bam = $_;
	my $panel = 'CCP.bed';
	my $vcf = '81485-01-01.vcf';
	my $seed = int(rand(999999999999999));
	$seed = "control_N$seed";
	$work_q->enqueue( [$bam, $panel, $vcf, $seed] );
	}

close READ;

$work_q->end;
$_->join for threads->list;




