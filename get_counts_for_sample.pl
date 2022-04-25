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
		my $panel	= $passed->[1];
		my $vcf		= $passed->[2];
		my $current_dir	= $passed->[3];
		my $output	= $passed->[4];
		
		my $cmd = "perl $current_dir/HF_grep_var_count.pl $bam $panel $vcf > $output";
		`$cmd`;
		}
	}

sub generate_count_data {
	my $options = shift;
	my $vcf_array = shift;

	threads->create( \&worker_HF ) for 1 .. ($options->{threads});
	
	my $n = 0;
	foreach my $seed (keys %{$vcf_array}) {
		my $bam = $options->{sample};
		my $panel = $options->{panel};
		my $vcf = $vcf_array->{$seed};
		my $current_dir = $options->{current_dir};
		my $output = $options->{test_folder}."/$seed.data";
		$work_HF->enqueue( [$bam, $panel, $vcf, $current_dir, $output] );
		}
	$work_HF->end;
	$_->join for threads->list;
	}

sub write_output {
	my $options = shift;
	my $vcf_array = shift;
	
	open (WRITE, ">".$options->{output});
	close WRITE;

	foreach my $seed (keys %{$vcf_array}) {
		my $data = $options->{test_folder}."/$seed.data";
		my $cmd = "cat $data >> ".$options->{output};
		`$cmd`;
		}
	}

sub get_vcf_array {
	my $options = shift;
	
	my $vcf_array = {};
	open (VCF, "<".$options->{vcf});
	
	my $total_count = 0;
	my %var_list;
	while (<VCF>) {
		chomp;
		next if m!^#!;
		my @mas = split/\t/;
		my @alt_array = split/,/,$mas[4];
		foreach my $alt (@alt_array) {
			my $var = "$mas[0]:$mas[1]".uc($mas[3]).">".uc($alt);
			$var_list{$var} = 1;
			$total_count += 1;
			}
		}

	close VCF;

	if ((scalar(keys %var_list)) > 100) {
		my $portion = int((scalar(keys %var_list))/($options->{threads}));
		
		my $seed = generate_seed();
		my $vcf = "".$options->{test_folder}."/$seed.vcf";
		$vcf_array->{$seed} = $vcf;
		my $count = 0;
		my $added = 0;
		open (VCF, ">$vcf");
		foreach my $var (sort {$a cmp $b} keys %var_list) {
			if ($var =~ /(\S+):(\d+)(\S+)>(\S+)/) {
				print VCF "$1\t$2\t.\t$3\t$4\t.\t.\t.\n";
				$count += 1;
				$added += 1;
				if ($count > $portion) {
					if ($total_count - $added < 120) {
						} else {
						close VCF;
						$count = 0;
						my $seed = generate_seed();
						$vcf = "".$options->{test_folder}."/$seed.vcf";
						$vcf_array->{$seed} = $vcf;
						open (VCF, ">$vcf");
						}
					}
				}
			}
		close VCF;
		} else {
		my $seed = generate_seed();
		$vcf_array->{$seed} = ($options->{vcf});
		}
	return $vcf_array;
	}

sub run {
	my $options = shift;
	my $vcf_array = get_vcf_array($options);
	generate_count_data($options, $vcf_array);
	write_output($options, $vcf_array);
	}

sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'        => \$opts{'h'},
		'v|vcf=s'   => \$opts{'vcf'},
		's|sample=s'   => \$opts{'sample'},
		'o|output=s' => \$opts{'output'},
		'p|panel=s' => \$opts{'panel'},
		'n|output=s' => \$opts{'threads'},
		'seed|seed=s'   => \$opts{'seed'}
	);
	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'vcf'});
	pod2usage(1) if(!$opts{'sample'});
	pod2usage(1) if(!$opts{'output'});
	pod2usage(1) if(!$opts{'panel'});
	$opts{threads} = 1 unless defined $opts{threads};
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

AOD_AMPL_CALL COUNTER

=head1 SYNOPSIS

CREATE DISTRIBUTION PARAMETERS

Options:

    -s  [REQUIRED] - input .bam file
    -v  [REQUIRED] - list of target variants in .vcf format
    -o  [REQUIRED] - output file
    -p  [REQUIRED] - path to bed file (amplicon panel)
    -n  {REQUIRED] - number of threads















