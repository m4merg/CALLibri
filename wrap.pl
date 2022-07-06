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
use List::MoreUtils qw(uniq);
use Thread::Queue;

my $current_dir = __DIR__;
my $test_folder = "$current_dir/test";
my $log_file = "$test_folder/log";
my $n_threads = 13;
my $list_control = $ARGV[0]; # bamListHRDRef
#my $sample_bam = $ARGV[1]; # /home/onco-admin/RnD/UEBAcall/5099.bam
my $input_panel = $ARGV[1]; # /home/onco-admin/ATLAS_software/aod-pipe/panel_info/AODHRD15/AODHRD15.designed.bed
my $input_vcf = $ARGV[2]; # test.vcf
my $panel_size = $ARGV[3];
$panel_size = 1 unless defined $panel_size;

open (my $log_fh, ">$log_file");

my $qscore_averaging_range      = 1; # Phred quality score is average in this window/ This value defines half of the window length.
my $minimum_coverage = 2; # Positions with coverage lower this value will be ignored (defined as non-detectable)

my $work_HF   = Thread::Queue->new;
my $work_YU   = Thread::Queue->new;
my $work_PPB  = Thread::Queue->new;

sub generate_seed {
        my @set = ('0' ..'9', 'A' .. 'Z', 'a' .. 'z');
        my $str = join '' => map $set[rand @set], 1 .. 15;
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
		my $cmd = "perl $current_dir/HF_grep_var_count.pl $bam $panel $vcf > $test_folder/$seed";
		`$cmd`;
		}
	}

sub worker_YU {
	while ( my $passed = $work_YU->dequeue ) {
		my $seed	= $passed->[0];
		print $log_fh "Started YU $seed\n";
		my $cmd = "R --slave -f $current_dir/YU_beta_error_approx.R --args $current_dir/fitdistr/ $test_folder/YU_data_$seed 1 > $test_folder/YU_result_$seed";
		`$cmd`;
		}
	}

sub worker_PPB {
	while ( my $passed = $work_PPB->dequeue ) {
		my $seed	= $passed->[0];
		my $counter	= $passed->[1];
		print $log_fh "Started PPB $seed $counter\n";
		my $cmd = `R --slave -f $current_dir/lib/ppb.r --args $test_folder/ppb_$seed.in.$counter $test_folder/ppb_$seed.out.total.p$counter $test_folder/ppb_$seed.out.detailed.p$counter 2> $test_folder/ppb_$seed.log.p$counter`;
		`$cmd`;
		}
	}

threads->create( \&worker_HF ) for 1 .. $n_threads;

# Generating count data

open (READ, "<$list_control");

my $n = 0;
my $sample_data;
my %sample_file;
while (<READ>) {
	chomp;
	next if m!^#!;
	my @mas = split/\t/;
	if (defined($sample_file{$mas[0]})) {
		die "Sample name '$mas[0]' met twice in input\n";
		}
	my $bam = $mas[1];
	my $panel = $input_panel;
	my $vcf = $input_vcf;
	my $seed = generate_seed();
	$seed = "bam_data_N$seed";
	$sample_file{$mas[0]} = $seed;
	$work_HF->enqueue( [$bam, $panel, $vcf, $seed] );
	}

close READ;

$work_HF->end;
$_->join for threads->list;

#exit;

# Reading  output files with count data with read counts and forming inner data structure

my %job_list;
foreach my $seed (values %sample_file) {
	open (CFILEINPUT, "<$test_folder/$seed");
	$sample_data->{({reverse %sample_file}->{$seed})} = {};
	while (<CFILEINPUT>) {
		chomp;
		my @mas = split/\t/;
		my $weight = log($mas[3] + 2.7183)*log($mas[4] + 2.7183)*log($mas[4] + 2.7183);
		chomp $weight;
		$job_list{"$mas[0]\@$mas[1]\@$mas[2]"} = 1;
		$sample_data->{({reverse %sample_file}->{$seed})}->{"$mas[0]\@$mas[1]\@$mas[2]"} = {"altCnt" => $mas[3], "depth" => $mas[4], "weight" => $weight};
		}
	close CFILEINPUT;
	}


# fitting distributions

threads->create( \&worker_YU ) for 1 .. $n_threads;

my %indexData;
my %beta;
foreach my $index (uniq(map {$_ = substr($_, 0, index($_, '@'))} keys %job_list)) {
	foreach my $job_element (grep(/$index/, (keys %job_list))) {
		my $amplicon;
		my $strand;
		if ($job_element =~ /(\S+)@(\S+)@(\S+)/) {
			$amplicon = $2;
			$strand = $3;
			}
		my $seed = generate_seed();
		$beta{$seed} = {};
		$beta{$seed}->{amplicon} = $amplicon;
		$beta{$seed}->{index} = $index;
		$beta{$seed}->{strand} = $strand;
		open (SAMPLEOUTPUT, ">$test_folder/YU_data_$seed");
		
		foreach my $sample (keys %{$sample_data}) {
			next unless defined $sample_data->{$sample}->{$job_element};
			print SAMPLEOUTPUT "",$sample_data->{$sample}->{$job_element}->{altCnt},"\t",$sample_data->{$sample}->{$job_element}->{depth},"\t",$sample_data->{$sample}->{$job_element}->{weight},"\n";
			}
		close SAMPLEOUTPUT;

		$work_YU->enqueue( [$seed] );
		}
	}

$work_YU->end;
$_->join for threads->list;

foreach my $seed (keys %beta) {
	print $log_fh "",$beta{$seed}->{index}," / AMP: ",$beta{$seed}->{amplicon}," / STRAND: ",$beta{$seed}->{strand},"$seed\n";

	open (YURES, "<$test_folder/YU_result_$seed");
	my $line;
	while (<YURES>) {
		$line = $_;
		chomp $line;
		print $log_fh "$line\n";
		}
	close YURES;
	
	my $key = $beta{$seed}->{index}.'@'.$beta{$seed}->{amplicon}.'@'.$beta{$seed}->{strand};
	$beta{$key} = {};
	$beta{$key}->{alpha} = ((split/\t/,$line)[0] || "NA");
	$beta{$key}->{beta}  = ((split/\t/,$line)[1] || "NA");
	$beta{$key}->{mean}  = ((split/\t/,$line)[2] || "NA");
	$beta{$key}->{alpha} =~ s/^\s+|\s+$//g;
	$beta{$key}->{beta}  =~ s/^\s+|\s+$//g;
	$beta{$key}->{mean}  =~ s/^\s+|\s+$//g;
	$beta{$key}->{seed} = $seed;
	}

my $pval_calc_seed = generate_seed();

my %group_seeds;
my $group_seed_count = 0;
my $group_seed_counter = 1;
open (my $pval_calc_fh, ">$test_folder/ppb_$pval_calc_seed.in.$group_seed_counter");
threads->create( \&worker_PPB ) for 1 .. $n_threads;

foreach my $sample (keys %{$sample_data}) {
	foreach my $index (uniq(map {$_ = substr($_, 0, index($_, '@'))} keys %job_list)) {
		my $pval = [];
		my $ad = [];
		my $dp = [];
		my $seed = generate_seed();
		$group_seeds{$seed} = {'index' => $index, 'sample' => $sample};
		foreach my $job_element (grep(/$index/, (keys %job_list))) {
			my $altCnt;
			my $depth;
			unless (defined $sample_data->{$sample}->{$job_element}) {
				$altCnt = "NA";
				$depth = "NA";
				} else {
				$altCnt = $sample_data->{$sample}->{$job_element}->{altCnt};
				$depth  = $sample_data->{$sample}->{$job_element}->{depth};
				}
			my $alpha_val  = $beta{$job_element}->{alpha};
			my $beta_val   = $beta{$job_element}->{beta};
			my $mean_val   = $beta{$job_element}->{mean};
			#my $cmd = "R --slave -f $current_dir/lib/ppb.r --args $altCnt $depth $alpha_val $beta_val $mean_val $panel_size";
			print $pval_calc_fh "$seed\t$altCnt\t$depth\t$alpha_val\t$beta_val\t$mean_val\t1\n";
			print $log_fh "$index\t$sample\t$seed\t$altCnt\t$depth\t$alpha_val\t$beta_val\t$mean_val\t1\n";
			++$group_seed_count;
			}
		if ($group_seed_count > 3000) {
			close $pval_calc_fh;
			$work_PPB->enqueue( [$pval_calc_seed, $group_seed_counter] );
			$group_seed_count = 0;
			$group_seed_counter += 1;
			open ($pval_calc_fh, ">$test_folder/ppb_$pval_calc_seed.in.$group_seed_counter");

			}
		}
	}
close $pval_calc_fh;
if ($group_seed_count > 0) {
	$work_PPB->enqueue( [$pval_calc_seed, $group_seed_counter] );
	}

$work_PPB->end;
$_->join for threads->list;

`cat $test_folder/ppb_$pval_calc_seed.out.total.p* > $test_folder/ppb_$pval_calc_seed.out.total`;
`cat $test_folder/ppb_$pval_calc_seed.out.detailed.p* > $test_folder/ppb_$pval_calc_seed.out.detailed`;

my @pval_by_group;
open (PVALBYGROUP, "<$test_folder/ppb_$pval_calc_seed.out.detailed");

while (<PVALBYGROUP>) {
	chomp;
	my @mas = split/\t/;
	my $pval = $mas[3];
	$pval = ((-1)*int(10*log($pval)/log(10))/1) unless $pval eq 'NA';
	my $ad = $mas[1];
	$ad = int($ad) unless $ad eq 'NA';
	my $alpha = $mas[4];
	my $beta = $mas[5];
	$alpha = int(1000*$alpha)/1000 unless $alpha eq 'NA';
	$beta = int(1000*$beta)/1000 unless $beta eq 'NA';
	push @pval_by_group, {'seed' => $mas[0], "AD" => $ad, "DP" => $mas[2], "P" => $pval, "A" => $alpha, "B" => $beta};
	}

close PVALBYGROUP;

open (PVALTOTAL, "<$test_folder/ppb_$pval_calc_seed.out.total");

while (<PVALTOTAL>) {
	chomp;
	my @mas = split/\t/;
	my $pval = $mas[1];
	$pval = ((-1) * int(10*log(min(1, ($pval * $panel_size)))/log(10))/1) unless $pval eq 'NA';
	my $i = 0;
	my $index  = $group_seeds{$mas[0]}->{index};
	my $sample = $group_seeds{$mas[0]}->{sample};
	print "$sample\t$index\t$pval\t",join(';', map {++$i; "AODAD$i=".$_->{AD}.",".$_->{DP}.";AODP$i=".$_->{P}.";AODA$i=".$_->{A}.";AODB$i=".$_->{B}} (grep {$_->{seed} eq $mas[0]} @pval_by_group)),"\n";
        }

close PVALTOTAL;


exit();




