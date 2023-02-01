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

my @knownTags = qw(ShortOverlap);

sub generate_seed {
        my @set = ('0' ..'9', 'A' .. 'Z', 'a' .. 'z');
        my $str = join '' => map $set[rand @set], 1 .. 15;
        return $str
        }

sub pass_by_tag {
	my $options = shift;
	my $info = shift;
	my $tag = shift;
	my @info_mas = split/;/, $info;
	my $filter = 1;
	if ($tag eq 'ALL') {$tag = ''} else {$tag = "_$tag"}
	
	my $is_there_not_NA = 0;
	my @p_mas;
	foreach my $arg (@info_mas) {
		if ($arg =~ /AODP(\d+)$tag=(\d+|NA)/) {
			my $number = 1;
			my $local = $2;
			if ($local eq 'NA') {
				$local = 0;
				} 
			push @p_mas, $local;
			if ($local ne 'NA') {$is_there_not_NA = 1} else {next}
			if ($local < $options->{local}) {
				$filter = 0;
				}
			}
		}
	#my $p_string = join(" ", map {10**((-1)*($_/10))} @p_mas);
	my $p_string = join(" ", @p_mas);
	my $current_dir = $options->{current_dir};
	if ($p_string eq '') {return 0}
	my $p_global = `R --slave -f $current_dir/lib/Fisher.r --args $p_string`; chomp $p_global;
	if ($p_global > 0) {
		#$p_global = (-10)*log($p_global)/log(10)
       		} else {
		$p_global = $options->{global} + 10;
		}
	if (($filter eq 1)or($p_global > $options->{global})) {
		return $p_global;
		} else {
		return 0;
		}
	$filter = 0 if $is_there_not_NA eq 0;
	return $filter;
	}

sub get_AF_by_tag {
	my $options = shift;
	my $info = shift;
	my $tag = shift;
	my @info_mas = split/;/, $info;
	my $closest = 2;
	my $result = '0';
	if ($tag eq 'ALL') {$tag = ''} else {$tag = "_$tag"}
	
	foreach my $arg (@info_mas) {
		if ($arg =~ /AODAD\d+$tag=(\d+),(\d+)/) {
			my $ad = $1;
			my $dp = $2;
			my $af;
			if ($dp eq 0) {$af = 0} else {$af = $ad/$dp}
			if (abs($af - 0.5) < $closest) {
				$result = $af;
				$closest = abs($af - 0.5);
				}
			}
		}
	return int(10000*$result)/10000;
	}

sub get_element_from_info {
	my $info = shift;
	my $element = shift;
	my @info_mas = split/;/, $info;
	foreach my $arg (@info_mas) {
		if ($arg =~ /$element=/) {
			return $arg;
			}
		}
	return undef;
	}

sub filter_vcf {
	my $options = shift;

	my $output;
	my $input = $options->{input};
	if (defined($options->{output})) {
		$output = $options->{output}
		} else {
		$output = $options->{test_folder} . "/" . $options->{seed} . ".vcf";
		}

	open (VCF, "<$input");
	open (WRITE, ">$output");
	
	while (<VCF>) {
		chomp;
		if (m!^#!) {
			print WRITE "$_\n";
			next;
			}
		my @mas = split/\t/;
		#next unless $mas[1] eq '108114661';
		#next unless $mas[3] eq 'AT';
		my $global = $mas[5];
		my $filter = 'FAIL';
		my @info = split/;/, $mas[7];
		if (pass_by_tag($options, $mas[7], 'CLEAR')) {
			$filter = 'PASS'
			} else {
			my @passed_tags;
			foreach my $tag (@knownTags) {
				if (pass_by_tag($options, $mas[7], $tag)) {
					push @passed_tags, $tag;
					}
				}
			if ((scalar @passed_tags) > 0) {
				$filter = join(';', @passed_tags);
				}
			}
		for (my $i = 1; $i < 10; $i++) {
			if (defined(get_element_from_info($mas[7], "AODA$i"))) {
				#push @info, get_element_from_info($mas[7], "AODA$i");
				}
			if (defined(get_element_from_info($mas[7], "AODB$i"))) {
				#push @info, get_element_from_info($mas[7], "AODB$i");
				}
			}
		if ($mas[7] =~ /AODAF=/) {
			$mas[7] = join(";", @info);
			} else {
			$mas[7] = join(";", (@info, 'AODAF='.get_AF_by_tag($options, $mas[7], 'ALL')));
			}

		$mas[6] = $filter;
		my $out_string = join("\t", @mas);
		print WRITE "$out_string\n";
		}
	
	close VCF;
	close WRITE;
	unless (defined($options->{output})) {
		if (-s $output) {
			my $line_count_input = `wc -l < $input`;chomp $line_count_input;
			my $line_count_output = `wc -l < $output`;chomp $line_count_output;
			if ($line_count_input eq $line_count_output) {
				`mv $output $input`;
				} else {print STDERR "FILTER FAILED\n";}
			} else {print STDERR "FILTER FAILED\n"}
		}
	}

sub run {
	my $options = shift;
	filter_vcf($options);
	}

sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'        => \$opts{'h'},
		'input|input=s'   => \$opts{'input'},
		'output|output=s' => \$opts{'output'},
		'global|global=s' =>\$opts{'global'},
		'local|local=s' =>\$opts{'local'},
		'seed|seed=s'   => \$opts{'seed'}
	);
	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'input'});
	if ((not(defined($opts{global})))or(not($opts{global} =~ /^-?\d+$/))) {
		$opts{global} = 70;
		}
	if ((not(defined($opts{local})))or(not($opts{local} =~ /^-?\d+$/))) {
		$opts{local} = 30;
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

VCF FILTER

=head1 SYNOPSIS

FILTER VCF FILE GENERATED BY makeVCF.pl

Options:

    -input   [REQUIRED] - input .vcf file
    -output  [OPTIONAL] - output .vcf file (input will be rewrited if not defined)
    -global  [OPTIONAL] - threshold for global p-value (DEFAULT:70)
    -local   [OPTIONAL] - threshold for p-value for each read group (DEFAULT:30)




















