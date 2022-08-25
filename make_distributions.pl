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


my $work_HF   = Thread::Queue->new;
my $work_YU   = Thread::Queue->new;

sub generate_seed {
        my @set = ('0' ..'9', 'A' .. 'Z', 'a' .. 'z');
        my $str = join '' => map $set[rand @set], 1 .. 15;
        return $str
        }

sub worker_HF {
	while ( my $passed = $work_HF->dequeue ) {
		my $bam		= $passed->[0];
		my $seed	= $passed->[1];
		my $panel	= $passed->[2];
		my $vcf		= $passed->[3];
		my $test_folder	= $passed->[4];
		my $current_dir	= $passed->[5];
		
		my $cmd = "perl $current_dir/HF_grep_var_count.pl $bam $panel $vcf | grep 'ALL' > $test_folder/$seed";
		`$cmd`;
		}
	}

sub worker_YU {
	while ( my $passed = $work_YU->dequeue ) {
		my $seed	= $passed->[0];
		my $test_folder = $passed->[1];
		my $current_dir	= $passed->[2];
		
		my $cmd = "R --slave -f $current_dir/YU_beta_error_approx.R --args $current_dir/fitdistr/ $test_folder/YU_data_$seed 1 > $test_folder/YU_result_$seed";
		`$cmd`;
		}
	}

sub generate_count_data {
	my $options = shift;

	threads->create( \&worker_HF ) for 1 .. ($options->{threads});
	open (READ, "<".$options->{list});
	
	my $n = 0;
	my $sample_file;
	while (<READ>) {
		chomp;
		next if m!^#!;
		my @mas = split/\t/;
		if (defined($sample_file->{$mas[0]})) {
			die "Sample name '$mas[0]' met twice in input\n";
			}
		my $bam = $mas[1];
		my $seed = generate_seed();
		$seed = "bam_data_N$seed";
		$sample_file->{$mas[0]} = $seed;
		my $panel = $options->{panel};
		my $vcf = $options->{target_vcf};
		my $test_folder = $options->{test_folder};
		my $current_dir = $options->{current_dir};
		$work_HF->enqueue( [$bam, $seed, $panel, $vcf, $test_folder, $current_dir] );
		}
	
	close READ;
	
	$work_HF->end;
	$_->join for threads->list;

	return $sample_file;
	}

sub read_count_data {
	my $options = shift;
	my $sample_file = shift;

	my $job_list;
	my $sample_data;
	foreach my $seed (values %{$sample_file}) {
		open (CFILEINPUT, "<".$options->{test_folder}."/$seed");
		$sample_data->{({reverse %{$sample_file}}->{$seed})} = {};
		while (<CFILEINPUT>) {
			chomp;
			my @mas = split/\t/;
			my $weight = log($mas[3] + 2.7183)*log($mas[4] + 2.7183)*log($mas[4] + 2.7183);
			chomp $weight;
			$job_list->{"$mas[0]\@$mas[1]\@$mas[2]"} = 1;
			$sample_data->{({reverse %{$sample_file}}->{$seed})}->{"$mas[0]\@$mas[1]\@$mas[2]"} = {"altCnt" => $mas[3], "depth" => $mas[4], "weight" => $weight};
			}
		close CFILEINPUT;
		}

	return ($sample_data, $job_list);
	}

sub fitting_distributions {
	my $options = shift;
	my $sample_data = shift;
	my $job_list = shift;
	
	my $beta;
	foreach my $index (uniq(map {$_ = substr($_, 0, index($_, '@'))} keys %{$job_list})) {
		foreach my $job_element (grep(/$index/, (keys %{$job_list}))) {
			my $amplicon;
			my $strand;
			if ($job_element =~ /(\S+)@(\S+)@(\S+)/) {
				$amplicon = $2;
				$strand = $3;
				}
			my $seed = generate_seed();
			$beta->{$seed} = {};
			$beta->{$seed}->{amplicon} = $amplicon;
			$beta->{$seed}->{index} = $index;
			$beta->{$seed}->{strand} = $strand;
			open (SAMPLEOUTPUT, ">".$options->{test_folder}."/YU_data_$seed");
			
			foreach my $sample (keys %{$sample_data}) {
				next unless defined $sample_data->{$sample}->{$job_element};
				print SAMPLEOUTPUT "",$sample_data->{$sample}->{$job_element}->{altCnt},"\t",$sample_data->{$sample}->{$job_element}->{depth},"\t",$sample_data->{$sample}->{$job_element}->{weight},"\n";
				}
			close SAMPLEOUTPUT;
			}
		}
	
	$sample_data = undef;
	$job_list = undef;

	threads->create( \&worker_YU ) for 1 .. ($options->{threads});

	foreach my $seed (keys %{$beta}) {
		my $test_folder = $options->{test_folder};
		my $current_dir = $options->{current_dir};
		$work_YU->enqueue( [$seed, $test_folder, $current_dir] );
		}

	$work_YU->end;
	$_->join for threads->list;

	return $beta;
	}

sub read_distributions {
	my $options = shift;
	my $beta = shift;

	my $log_fh = $options->{log_fh};
	foreach my $seed (keys %{$beta}) {
		print $log_fh "",$beta->{$seed}->{index}," / AMP: ",$beta->{$seed}->{amplicon}," / STRAND: ",$beta->{$seed}->{strand},"$seed\n";
	
		open (YURES, "<".$options->{test_folder}."/YU_result_$seed");
		my $line;
		while (<YURES>) {
			$line = $_;
			chomp $line;
			print $log_fh "$line\n";
			}
		close YURES;
		
		my $key = $beta->{$seed}->{index}.'@'.$beta->{$seed}->{amplicon}.'@'.$beta->{$seed}->{strand};
		$beta->{$key} = {};
		$beta->{$key}->{alpha} = ((split/\t/,$line)[0] || "NA");
		$beta->{$key}->{beta}  = ((split/\t/,$line)[1] || "NA");
		$beta->{$key}->{mean}  = ((split/\t/,$line)[2] || "NA");
		$beta->{$key}->{alpha} =~ s/^\s+|\s+$//g;
		$beta->{$key}->{beta}  =~ s/^\s+|\s+$//g;
		$beta->{$key}->{mean}  =~ s/^\s+|\s+$//g;
		$beta->{$key}->{seed} = $seed;
		}
	return $beta;
	}

sub get_target_vcf {
	my $options = shift;
	
	if ($options->{'mode'} eq 'create') {
		$options->{target_vcf} = $options->{vcf};
		return 0;
		}

	my $target_vcf = $options->{test_folder}."/".generate_seed().".vcf";
	if (open(BDATA_FH, "<".$options->{bdata})) {
		if ((-s ($options->{bdata})) > 0) {
			my %known;
			while (<BDATA_FH>) {
				chomp;
				my @mas = split/\t/;
				my @target = split/\@/, $mas[0];
				my $var = uc($target[0]);
				if ($var =~ /(\S+):(\d+)(\S+)>(\S+)/) {
					$known{$var} = 1;
					} else {
					die "Unknown Format of input .bdata file"
					}
				}
			open (READ_VCF, "<".$options->{vcf});
			open (WRITE_VCF, ">$target_vcf");
			while (<READ_VCF>) {
				chomp;
				next if m!^#!;
				my @mas = split/\t/;
				my @alt_array = split/,/, $mas[4];
				foreach my $alt (@alt_array) {
					my $var = uc("$mas[0]:$mas[1]$mas[3]>$alt");
					next if defined $known{$var};
					print WRITE_VCF "$mas[0]\t$mas[1]\t.\t$mas[3]\t$alt\t.\t.\t.\n";
					}
				}
			close READ_VCF;
			close WRITE_VCF;
			$options->{target_vcf} = $target_vcf;
			} else {
			$options->{target_vcf} = $options->{vcf};
			}
		close BDATA_FH;	
		} else {
		$options->{target_vcf} = $options->{vcf};
		}
	}

sub bdata_write {
	my $bdata_write = shift;
	my $beta = shift;
	foreach my $group (sort {$a cmp $b} keys %{$beta}) {
		next if defined $beta->{$group}->{index};
		my $alpha_p = $beta->{$group}->{alpha};
		my $beta_p = $beta->{$group}->{beta};
		my $mean_p = $beta->{$group}->{mean};
		print $bdata_write "$group\t$alpha_p\t$beta_p\t$mean_p\n";
		}
	}

sub print_distributions {
	my $options = shift;
	my $beta = shift;
	
	if ($options->{mode} eq 'create') {
		open (my $bdata_write, ">".$options->{bdata});
		bdata_write($bdata_write, $beta);
		close $bdata_write;
		return 0;
		}

	if (open(BDATA_FH, "<".$options->{bdata})) {
		if ((-s ($options->{bdata})) > 0) {
			open (my $bdata_write, ">>".$options->{bdata});
			bdata_write($bdata_write, $beta);
			close $bdata_write;
			} else {
			open (my $bdata_write, ">".$options->{bdata});
			bdata_write($bdata_write, $beta);
			close $bdata_write;
			}
		} else {
		open (my $bdata_write, ">".$options->{bdata});
		bdata_write($bdata_write, $beta);
		close $bdata_write;
		}
	}

sub run {
	my $options = shift;
	get_target_vcf($options);
	if ((-s ($options->{target_vcf})) eq 0) {
		return 0;
		}
	my $sample_file = generate_count_data($options);
	my ($sample_data, $job_list) = read_count_data($options, $sample_file);
	my $beta = fitting_distributions($options, $sample_data, $job_list);
	$beta = read_distributions($options, $beta);
	print_distributions($options, $beta);
	}

sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'        => \$opts{'h'},
		'v|vcf=s'   => \$opts{'vcf'},
		'l|list=s'   => \$opts{'list'},
		'n|threads=s' => \$opts{'threads'},
		'bdata|output=s'  => \$opts{'bdata'},
		'p|panel=s' => \$opts{'panel'},
		'mode|mode=s' => \$opts{'mode'},
		'seed|seed=s'   => \$opts{'seed'}
	);
	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'vcf'});
	pod2usage(1) if(!$opts{'list'});
	pod2usage(1) if(!$opts{'bdata'});
	pod2usage(1) if(!$opts{'panel'});
	unless (defined($opts{mode})) {
		$opts{mode} = 'append';
		}
	$opts{mode} = lc($opts{mode});
	if (($opts{mode} eq 'append')or($opts{'mode'} eq 'create')) {
		} else {
		pod2usage(1);
		}
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

AOD_AMPL_CALL DISTRIBUTION MAKER

=head1 SYNOPSIS

CREATE DISTRIBUTION PARAMETERS

Options:

    -l  [REQUIRED] - list of input bam files (either first or second collumn is path to .bam file)
    -v  [REQUIRED] - list of target variants in .vcf format
    -p  [REQUIRED] - path to bed file (amplicon panel)
    -n  [OPTIONAL] - number of threads to be used
    -bdata  [REQUIRED] - input/output file with site-specific distribution parameters
    -mode  [OPTIONAL] - CREATE/APPEND[DEFAULT]. Create will overwrite existing data in .bdata output file. Append will only append new sites from input .vcf file into output file















