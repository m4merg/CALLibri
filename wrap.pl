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
use List::Util qw(min);
use Thread::Queue;

my $current_dir = __DIR__;
my $test_folder = "$current_dir/test";
my $log_file = "$test_folder/log";
my $n_threads = 13;
my $list_control = $ARGV[0]; # bamListHRDRef
my $sample_bam = $ARGV[1]; # /home/onco-admin/RnD/UEBAcall/5099.bam
my $input_panel = $ARGV[2]; # /home/onco-admin/ATLAS_software/aod-pipe/panel_info/AODHRD15/AODHRD15.designed.bed
my $input_vcf = $ARGV[3]; # test.vcf

open (my $log_fh, ">$log_file");

my $qscore_averaging_range      = 1; # Phred quality score is average in this window/ This value defines half of the window length.
my $minimum_coverage = 2; # Positions with coverage lower this value will be ignored (defined as non-detectable)

my $work_HF   = Thread::Queue->new;
my $work_YU   = Thread::Queue->new;

sub generate_seed {
        my @set = ('0' ..'9', 'A' .. 'Z', 'a' .. 'z');
        my $str = join '' => map $set[rand @set], 1 .. 20;
        return $str
        }

sub worker_HF {
	while ( my $passed = $work_HF->dequeue ) {
		my $bam		= $passed->[0];
		my $panel	= $passed->[1];
		my $vcf		= $passed->[2];
		my $seed	= $passed->[3];
		#$bam = "/home/onco-admin/RnD/UEBAcall/81485-01-01.bam";
		print $log_fh "Started HF $bam\n";
		my $cmd = "perl $current_dir/HF_grep_var_count.pl $bam $panel $vcf > test/$seed";
		`$cmd`;
		}
	}

sub worker_YU {
	while ( my $passed = $work_YU->dequeue ) {
		my $seed	= $passed->[0];
		my $ADobs	= $passed->[1];
		my $DPobs	= $passed->[2];
		print $log_fh "Started YU $seed\n";
		my $cmd = "R --slave -f $current_dir/YU_beta_error_approx.R --args $current_dir/fitdistr/ $test_folder/YU_data_$seed 1 $ADobs $DPobs > $test_folder/YU_result_$seed";
		`$cmd`;
		}
	}

threads->create( \&worker_HF ) for 1 .. $n_threads;

open (READ, "<$list_control");

my $n = 0;
my @control;
while (<READ>) {
	chomp;
	next if m!^#!;
	my $bam = $_;
	my $panel = $input_panel;
	my $vcf = $input_vcf;
	my $seed = generate_seed();
	$seed = "control_N$seed";
	push(@control, $seed);
	$work_HF->enqueue( [$bam, $panel, $vcf, $seed] );
	}

close READ;

my $panel = $input_panel;
my $vcf = $input_vcf;
my $sampleSeed = generate_seed();
$sampleSeed = "sample_N$sampleSeed";
$work_HF->enqueue( [$sample_bam, $panel, $vcf, $sampleSeed] );

$work_HF->end;
$_->join for threads->list;

#exit;

# Reading  output files with read counts and forming inner data structure

my $controlData = [];
foreach my $cFile (@control) {
	open (CFILEINPUT, "<$test_folder/$cFile");
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
		#my $weight = `R --slave -f lib/get_weight.r --args $mas[3] $mas[4] 0.05`;
		my $weight = log($mas[3] + 2.7183)*log($mas[4] + 2.7183)*log($mas[4] + 2.7183);
		chomp $weight;
		$data->{weight} = $weight;
		push @{$controlData}, $data;
		}
	close CFILEINPUT;
	}

open (CFILEINPUT, "<$test_folder/$sampleSeed");

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
	#my $weight = `R --slave -f lib/get_weight.r --args $mas[3] $mas[4] 0.05`;
	my $weight = log($mas[3] + 2.7183)*log($mas[4] + 2.7183)*log($mas[4] + 2.7183);
	chomp $weight;
	$data->{weight} = $weight;
	push @{$sampleData}, $data;
	}

close CFILEINPUT;

# Creating input for R scripts - generating noize functions

#print Dumper $sampleData;
open (SAMPLEINPUT, "<$test_folder/$sampleSeed");

my $pval = {};
threads->create( \&worker_YU ) for 1 .. $n_threads;

my %YU_result;
while (<SAMPLEINPUT>) {
	chomp;
	my @mas = split/\t/;
	my $amplicon = $mas[1];
	my $index = $mas[0];
	my $strand = $mas[2];
	print $log_fh "$index - $amplicon - $strand\n";
	my $ADobs = [grep{($_->{amplicon} eq $amplicon)and($_->{index} eq $index)and($_->{strand} eq $strand)} @{$sampleData}]->[0]->{'altCnt'};
	my $DPobs = [grep{($_->{amplicon} eq $amplicon)and($_->{index} eq $index)and($_->{strand} eq $strand)} @{$sampleData}]->[0]->{'depth'};
	#print $log_file "$amplicon\t$index\t$strand\t$ADobs\n";
	my @data = grep{($_->{amplicon} eq $amplicon)and($_->{index} eq $index)and($_->{strand} eq $strand)} @{$controlData};
	my $seed = generate_seed();
	print $log_fh "$seed\n";
	open (SAMPLEOUTPUT, ">$test_folder/YU_data_$seed");
	
	foreach my $arg (@data) {
		print SAMPLEOUTPUT "",$arg->{altCnt},"\t",$arg->{depth},"\t",$arg->{weight},"\n";
		}
	
	close SAMPLEOUTPUT;
	
	$work_YU->enqueue( [$seed, $ADobs, $DPobs] );
	$YU_result{$index} = [] unless defined $YU_result{$index};
	push @{$YU_result{$index}}, {'amplicon' => $amplicon, 
					'index' => $index,
					'strand' => $strand,
					'seed' => $seed,
					'ADobs' => $ADobs,
					'DPobs' => $DPobs};
	#print $log_fh `R --slave -f $current_dir/YU_beta_error_approx.R --args $current_dir/fitdistr/ $test_folder/$seed 1 $ADobs $DPobs`;
	#my $p = `tail -n1 $log_file`;
	#chomp $p;
	#$pval->{$index} = [] unless defined $pval->{$index};
	#push @{($pval->{$index})}, $p;
	#print $log_fh "$index - $amplicon - $strand - $p\n";
	
	}

close SAMPLEINPUT;

$work_YU->end;
$_->join for threads->list;
print STDERR Dumper \%YU_result;

exit();

foreach my $index (keys %{$pval}) {
	my $p = join(' ', (sort {$a <=> $b} @{$pval->{$index}}));
	print STDERR "$index\t$p\n";
	$p = `R --slave -f $current_dir/lib/Fisher.r --args $p`;
	chomp $p;
	$p = min(1, ($p * 1200000));
	print "$index\t$p\n";
	}















