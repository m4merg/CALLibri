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
		#$bam = "/home/onco-admin/RnD/UEBAcall/81485-01-01.bam";
		print STDERR "Started $bam\n";
		my $cmd = "perl HF_grep_var_count.pl $bam $panel $vcf > test/$seed";
		#print STDERR "$cmd\n";
		`$cmd`;
		}
	}

threads->create( \&worker ) for 1 .. 13;

open (READ, "<bamListRef");

my $n = 0;
my @control;
while (<READ>) {
	chomp;
	next if m!^#!;
	my $bam = $_;
	my $panel = 'CCP.bed';
	my $vcf = 'test.vcf';
	my $seed = int(rand(999999999999999999));
	$seed = "control_N$seed";
	push(@control, $seed);
	$work_q->enqueue( [$bam, $panel, $vcf, $seed] );
	}

close READ;

my $panel = 'CCP.bed';
my $vcf = 'test.vcf';
my $sampleSeed = int(rand(999999999999999999));
$sampleSeed = "sample_N$sampleSeed";
my $bam = "/home/onco-admin/RnD/UEBAcall/81485-01-01.bam";
$work_q->enqueue( [$bam, $panel, $vcf, $sampleSeed] );

$work_q->end;
$_->join for threads->list;

#exit;

# Reading  output files with read counts and forming inner data structure

my $controlData = [];
foreach my $cFile (@control) {
	open (CFILEINPUT, "<test/$cFile");
	while (<CFILEINPUT>) {
		chomp;
		my @mas = split/\t/;
		my $data;
		$data->{seed} = $cFile;
		$data->{index} = $mas[0];
		$data->{amplicon} = $mas[1];
		$data->{strand} = $mas[2];
		$data->{altCnt} = $mas[3];
		$data->{depth} = $mas[4];
		my $weight = `R --slave -f lib/get_weight.r --args $mas[3] $mas[4] 0.05`;
		chomp $weight;
		$data->{weight} = $weight;
		push @{$controlData}, $data;
		}
	close CFILEINPUT;
	}

open (CFILEINPUT, "<test/$sampleSeed");

my $sampleData = [];
while (<CFILEINPUT>) {
	chomp;
	my @mas = split/\t/;
	my $data;
	$data->{seed} = $sampleSeed;
	$data->{index} = $mas[0];
	$data->{amplicon} = $mas[1];
	$data->{strand} = $mas[2];
	$data->{altCnt} = $mas[3];
	$data->{depth} = $mas[4];
	my $weight = `R --slave -f lib/get_weight.r --args $mas[3] $mas[4] 0.05`;
	chomp $weight;
	$data->{weight} = $weight;
	push @{$sampleData}, $data;
	}

close CFILEINPUT;

# Creating input for R scripts - generating noize functions

open (SAMPLEINPUT, "<test/$sampleSeed");

while (<SAMPLEINPUT>) {
	chomp;
	my @mas = split/\t/;
	my $amplicon = $mas[1];
	my $index = $mas[0];
	my $strand = $mas[2];
	my @data = grep{($_->{amplicon} eq $amplicon)and($_->{index} eq $index)and($_->{strand} eq $strand)} @{$controlData};
	my $seed = int(rand(999999999999999999));
	$seed = "indexdata_N$seed";
	open (SAMPLEOUTPUT, ">test/$seed");
	
	foreach my $arg (@data) {
		print SAMPLEOUTPUT "",$arg->{altCnt},"\t",$arg->{depth},"\t",$arg->{weight},"\n";
		}
	
	close SAMPLEOUTPUT;
	}

close SAMPLEINPUT;

















