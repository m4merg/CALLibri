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
use Getopt::Long 'GetOptions';
use Pod::Usage;


my $work_PPB  = Thread::Queue->new;
my @knownTags = qw(ShortOverlap);

sub generate_seed {
        my @set = ('0' ..'9', 'A' .. 'Z', 'a' .. 'z');
        my $str = join '' => map $set[rand @set], 1 .. 15;
        return $str
        }

sub worker_PPB {
	while ( my $passed = $work_PPB->dequeue ) {
		my $seed	= $passed->[0];
		my $counter	= $passed->[1];
		my $current_dir	= $passed->[2];
		my $test_folder	= $passed->[3];
		#print STDERR "STARTED THERAD\n";
		my $cmd = "R --slave -f $current_dir/lib/ppb.r --args $test_folder/ppb_$seed.in.$counter $test_folder/ppb_$seed.out.total.p$counter $test_folder/ppb_$seed.out.detailed.p$counter 2> $test_folder/ppb_$seed.log.p$counter";
		#print STDERR "$cmd\n";
		`$cmd`;
		}
	}


sub make_call {
	my $options = shift;
	my $sample_data = shift;
	my $beta = shift;
	my $job_list = shift;

	my $pval_calc_seed = generate_seed();
	$options->{pval_calc_seed} = $pval_calc_seed;
	
	my $test_folder = $options->{test_folder};
	my $current_dir = $options->{current_dir};

	my $group_seeds = {};
	my $group_seed_count = 0;
	my $group_seed_counter = 1;
	open (my $pval_calc_fh, ">".$options->{test_folder}."/ppb_$pval_calc_seed.in.$group_seed_counter");
	threads->create( \&worker_PPB ) for 1 .. ($options->{threads});
	#print Dumper $job_list;
	foreach my $index (uniq(map {$_ = substr($_, 0, index($_, '@'))} keys %{$job_list})) {
		#print "NEW index\n";
		foreach my $tag (@knownTags, "CLEAR", "ALL") {
		#foreach my $tag ("ALL") {
			#print "NEW tag\n";
			my $pval = [];
			my $ad = [];
			my $dp = [];
			my $seed = generate_seed();
	
			$group_seeds->{$seed} = {'index' => $index, "tag" => $tag};
			#$group_seeds->{$seed} = {'index' => $index};
			foreach my $job_element (grep(/$index\@/, (keys %{$job_list}))) {
				#print "NEW ELEMENT\n";
				my $job_element_tagged = "$job_element@" . $tag;
				my ($altCnt, $depth, $alpha_val, $beta_val, $mean_val);
				unless (defined $sample_data->{$job_element_tagged}) {
					$altCnt = "NA";$depth = "NA";
					} else {
					$altCnt = $sample_data->{$job_element_tagged}->{altCnt};
					$depth  = $sample_data->{$job_element_tagged}->{depth};
					}
				if (defined($beta->{$job_element})) {
					$alpha_val  = $beta->{$job_element}->{alpha};
					$beta_val   = $beta->{$job_element}->{beta};
					$mean_val   = $beta->{$job_element}->{mean};
					} else {
					$alpha_val  = 'NA';
					$beta_val   = 'NA';
					$mean_val   = 'NA';
					}
				#print "$seed\t$job_element\t$job_element_tagged\n";
				print $pval_calc_fh "$seed\t$altCnt\t$depth\t$alpha_val\t$beta_val\t$mean_val\t1\n";
				#print $log_fh "$index\t$seed\t$altCnt\t$depth\t$alpha_val\t$beta_val\t$mean_val\t1\n";
				++$group_seed_count;
				}
			if ($group_seed_count > 100) {
				close $pval_calc_fh;
				$work_PPB->enqueue( [$pval_calc_seed, $group_seed_counter, $current_dir, $test_folder] );
				$group_seed_count = 0;
				$group_seed_counter += 1;
				open ($pval_calc_fh, ">".$options->{test_folder}."/ppb_$pval_calc_seed.in.$group_seed_counter");
	
				}
			}
		}
	close $pval_calc_fh;
	if ($group_seed_count > 0) {
		$work_PPB->enqueue( [$pval_calc_seed, $group_seed_counter, $current_dir, $test_folder] );
		}
	$work_PPB->end;
	$_->join for threads->list;
	#print STDERR "WHAT?\n";
	my @cmd;
	push @cmd, "cat $test_folder/ppb_$pval_calc_seed.out.total.p* > $test_folder/ppb_$pval_calc_seed.out.total";
	push @cmd, "cat $test_folder/ppb_$pval_calc_seed.out.detailed.p* > $test_folder/ppb_$pval_calc_seed.out.detailed";
	foreach my $arg (@cmd) {
		#print STDERR "$arg\n";
		`$arg`;
		}
	#`cat $test_folder/ppb_$pval_calc_seed.out.total.p* > $test_folder/ppb_$pval_calc_seed.out.total`;
	#`cat $test_folder/ppb_$pval_calc_seed.out.detailed.p* > $test_folder/ppb_$pval_calc_seed.out.detailed`;
	return $group_seeds;
	}

sub print_results {
	my $options = shift;
	my $group_seeds = shift;
	my $pval_calc_seed = $options->{pval_calc_seed};
	#print Dumper $group_seeds;
	
	my @pval_by_group;
	open (PVALBYGROUP, "<".$options->{test_folder}."/ppb_$pval_calc_seed.out.detailed");
	
	while (<PVALBYGROUP>) {
		chomp;
		my @mas = split/\t/;
		my $pval = $mas[3];
		#$pval = ((-1)*int(10*log($pval)/log(10))/1) unless $pval eq 'NA';
		my $ad = $mas[1];
		$ad = int($ad) unless $ad eq 'NA';
		my $alpha = $mas[4];
		my $beta = $mas[5];
		$alpha = int(1000*$alpha)/1000 unless $alpha eq 'NA';
		$beta = int(1000*$beta)/1000 unless $beta eq 'NA';
		push @pval_by_group, {'seed' => $mas[0], "AD" => $ad, "DP" => $mas[2], "P" => $pval, "A" => $alpha, "B" => $beta};
		}
	
	close PVALBYGROUP;
	
	open (PVALTOTAL, "<".$options->{test_folder}."/ppb_$pval_calc_seed.out.total");
	
	my %uniq_indexes;
	my %printed_indexes;
	foreach my $key (keys %{$group_seeds}) {
		$uniq_indexes{$group_seeds->{$key}->{index}} = 1;
		}
	while (<PVALTOTAL>) {
		chomp;
		my @mas = split/\t/;
		next if ($group_seeds->{$mas[0]}->{tag} ne 'ALL');
		my $pval = $mas[1];
		#$pval = ((-1) * int(10*log(min(1, ($pval * ($options->{panel_size}))))/log(10))/1) unless $pval eq 'NA';
		my $index  = $group_seeds->{$mas[0]}->{index};
		next if (defined($printed_indexes{$index}));
		my @info;
		foreach my $key (keys %{$group_seeds}) {
			my $i = 0;
			next if $index ne $group_seeds->{$key}->{index};
			my $tag = $group_seeds->{$key}->{tag};
			my @info_local;
			if ($tag eq 'ALL') {
				$tag = '';
				@info_local = map {++$i;"AODAD${i}$tag=".$_->{AD}.",".$_->{DP}.";AODP${i}$tag=".$_->{P}.";AODA${i}=".$_->{A}.";AODB$i=".$_->{B}} (grep {$_->{seed} eq $key} @pval_by_group);
				} else {
				$tag = "_$tag" unless $tag eq 'ALL';
				@info_local = map {++$i;"AODAD${i}$tag=".$_->{AD}.",".$_->{DP}.";AODP${i}$tag=".$_->{P}} (grep {$_->{seed} eq $key} @pval_by_group);
				}

			#$tag = "_$tag" unless $tag eq 'ALL';
			#$tag = "" if $tag eq 'ALL';
			#my @info_local = map {++$i;"AODAD${i}$tag=".$_->{AD}.",".$_->{DP}.";AODP${i}$tag=".$_->{P}} (grep {$_->{seed} eq $key} @pval_by_group);
			push @info, @info_local;
			#foreach my $element (@pval_by_group) {
			#	next if ($element->{seed})
			#	}
			#print "$index\t",join(";", @info_local),"\n";
			#print "$index\t$pval\t",join(';', map {++$i; "AODAD${i}_$tag=".$_->{AD}.",".$_->{DP}.";AODP${i}_tag=".$_->{P}.";AODA${i}=".$_->{A}.";AODB$i=".$_->{B}} (grep {$_->{seed} eq $mas[0]} @pval_by_group)),"\n";
			}
		print "$index\t$pval\t",join(";", @info),"\n";
		$printed_indexes{$index} = 1;
	        }
	
	close PVALTOTAL;
	}

sub get_sample_data {
	my $options = shift;

	open (CDATA, "<".$options->{cdata});
	
	my $sample_data;
	while (<CDATA>) {
		chomp;
		my @mas = split/\t/;
		my $job_element = "$mas[0]\@$mas[1]\@$mas[2]\@$mas[5]";
		$sample_data->{$job_element} = {};
		$sample_data->{$job_element}->{altCnt} = $mas[3];
		$sample_data->{$job_element}->{depth} = $mas[4];
		}
	
	close CDATA;
	return $sample_data;
	}

sub get_beta {
	my $options = shift;

	open (BDATA, "<".$options->{bdata});
	
	my $beta;
	while (<BDATA>) {
		chomp;
		my @mas = split/\t/;
		$beta->{$mas[0]} = {};
		$beta->{$mas[0]}->{alpha} = $mas[1];
		$beta->{$mas[0]}->{beta}  = $mas[2];
		$beta->{$mas[0]}->{mean}  = $mas[3];
		}
	
	close BDATA;

	return $beta;
	}

sub get_job_list {
	my $options = shift;
	my $sample_data = shift;
	my $beta = shift;

	open (VCF, "<".$options->{vcf});
	
	my $job_list = {};
	my $cnt = 0;
	while (<VCF>) {
		chomp;
		next if m!^#!;
		my @mas = split/\t/;
		my @alt_array = split/,/,$mas[4];
		foreach my $alt (@alt_array) {
			++$cnt;
			my $index = "$mas[0]:$mas[1]".uc($mas[3]).">".uc($alt);
			foreach my $job_element (keys %{$sample_data}) {
				if ($job_element =~ /$index\@/) {
					my $job_element_untagged = $job_element;
					foreach my $tag (@knownTags, "ALL", 'CLEAR') {
						$job_element_untagged =~ s/\@$tag//;
						}
					$job_list->{$job_element_untagged} = 1;
					}
				}
			foreach my $job_element (keys %{$beta}) {
				if ($job_element =~ /$index\@/) {
					$job_list->{$job_element} = 1;
					}
				}
			}
		}
	
	close VCF;
	
	return $job_list;
	}

sub get_panel_size {
	my $options = shift;
	$options->{panel_size} = 1;
	}

sub run {
	my $options = shift;
	#print STDERR "HERE1\n";
	my $sample_data = get_sample_data($options);
	#print STDERR "HERE2\n";
	my $beta = get_beta($options);
	#print STDERR "HERE3\n";
	my $job_list = get_job_list($options, $sample_data, $beta);
	
	#print STDERR "HERE4\n";
	get_panel_size($options);
	#print STDERR "HERE5\n";
	my $group_seeds = make_call($options, $sample_data, $beta, $job_list);
	#print STDERR "HERE6\n";
	print_results($options, $group_seeds);
	}

sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'        => \$opts{'h'},
		'v|vcf=s'   => \$opts{'vcf'},
		'n|threads=s' => \$opts{'threads'},
            	'bdata|output=s'  => \$opts{'bdata'},
		'cdata|cdata=s' => \$opts{'cdata'},
		'seed|seed=s'   => \$opts{'seed'}
	);
	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'vcf'});
	pod2usage(1) if(!$opts{'bdata'});
	pod2usage(1) if(!$opts{'cdata'});
	unless (defined($opts{'threads'})) {
		$opts{'threads'} = 1;
		}
	return \%opts;
	}

{
	my $options = option_builder();
	my $current_dir = __DIR__;
	$options->{current_dir} = $current_dir;
	$options->{test_folder} = "$current_dir/test";
	$options->{log_file} = $options->{test_folder}."/log";
	open (my $log_fh, ">".$options->{log_file});
	$options->{log_fh} = $log_fh;
	$options->{seed} = generate_seed() unless defined($options->{seed});

	eval {run($options)};
	if ($@) {
		print STDERR "$@\n";
		} else {
		}

	close $log_fh;
}


__END__

=head1 NAME

AOD_AMPL_CALL SINGLE SAMPLE CALLER

=head1 SYNOPSIS

CALL VARIANTS FROM SINGLE SAMPLE

Options:

    -v  [REQUIRED] - list of target variants in .vcf format
    -n  [OPTIONAL] - number of threads to be used
    -bdata  [REQUIRED] - input/output file with site-specific distribution parameters
    -cdata  [REQUIRED] - input data with counts generated by 'get_counts_for_sample.pl'




















