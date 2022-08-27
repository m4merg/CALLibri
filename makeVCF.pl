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


sub generate_seed {
        my @set = ('0' ..'9', 'A' .. 'Z', 'a' .. 'z');
        my $str = join '' => map $set[rand @set], 1 .. 15;
        return $str
        }

sub create_vcf {
	my $options = shift;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

	open (CALL_DATA, "<".$options->{input});
	open (VCF, ">".$options->{output});
	
	my $sample = $options->{sample};
	my $sam_header = `samtools view -H $sample | grep "^\@SQ"`;
	my @contig;
	foreach my $arg (split/\n/, $sam_header) {
		my @mas = split/\t|\s+/, $arg;
		$mas[1] =~ s/^SN://;
		$mas[2] =~ s/^LN://;
		push @contig, "##contig=<ID=$mas[1],length=$mas[2],assembly=unknown>"
		}
	print VCF "##fileformat=VCFv4.1\n";
	$year = 1900 + $year;
	$mon  = (length($mon+1) eq 1 ? "0".($mon+1) : "".($mon+1)."");
	$mday = (length($mday) eq 1 ? "0$mday" : "$mday");
	print VCF "##fileDate=$year$mon$mday\n";
	print VCF "##source=AODcall\n";
	print VCF "##content=AODcall non-polymorphism variant calls\n";
	print VCF "",join("\n", @contig),"\n";
	print VCF "##INFO=<ID=AODAD,Number=.,Type=Integer,Description=\"Total allelic depths for the ref and alt alleles in the order listed across all read groups\">\n";

	#print VCF "##INFO=<ID=AODAD1,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed in the 1st group of reads\">\n";
	#print VCF "##INFO=<ID=AODAD2,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed in the 2d group of reads\">\n";
	#print VCF "##INFO=<ID=AODAD3,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed in the 3d group of reads\">\n";
	#print VCF "##INFO=<ID=AODAD4,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed in the 4th group of reads\">\n";
	
	foreach my $group (qw(1 2 3 4)) {
		print VCF "##INFO=<ID=AODAD$group,Number=.,Type=String,Description=\"Allelic depths for the ref and alt alleles in the order listed in the read group $group\">\n";
		print VCF "##INFO=<ID=AODP$group,Number=.,Type=String,Description=\"P value for variant across reads in the read group number $group\">\n";
		print VCF "##INFO=<ID=AODA$group,Number=.,Type=String,Description=\"Alpha value of beta distibution for read group $group\">\n";
		print VCF "##INFO=<ID=AODB$group,Number=.,Type=String,Description=\"Beta value of beta distibution for read group $group\">\n";
		}
	#print VCF "##INFO=<ID=AODP2,Number=.,Type=String,Description=\"P value for variant across reads in the 2d read group\">\n";
	#print VCF "##INFO=<ID=AODP3,Number=.,Type=String,Description=\"P value for variant across reads in the 3d read group\">\n";
	#print VCF "##INFO=<ID=AODP4,Number=.,Type=String,Description=\"P value for variant across reads in the 4th read group\">\n";
	
	#print VCF "##INFO=<ID=AODA1,Number=.,Type=String,Description=\"Alpha value of beta distibution for read group 1\">\n";
	#print VCF "##INFO=<ID=AODA2,Number=.,Type=String,Description=\"Alpha value of beta distibution for read group 2\">\n";
	#print VCF "##INFO=<ID=AODA3,Number=.,Type=String,Description=\"Alpha value of beta distibution for read group 3\">\n";
	#print VCF "##INFO=<ID=AODA4,Number=.,Type=String,Description=\"Alpha value of beta distibution for read group 4\">\n";
	
	#print VCF "##INFO=<ID=AODB1,Number=.,Type=String,Description=\"Beta value of beta distibution for read group 1\">\n";
	#print VCF "##INFO=<ID=AODB2,Number=.,Type=String,Description=\"Beta value of beta distibution for read group 2\">\n";
	#print VCF "##INFO=<ID=AODB3,Number=.,Type=String,Description=\"Beta value of beta distibution for read group 3\">\n";
	#print VCF "##INFO=<ID=AODB4,Number=.,Type=String,Description=\"Beta value of beta distibution for read group 4\">\n";
	
	print VCF "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";

	print VCF "##FILTER=<ID=PASS,Description=\"PASS\">\n";
	print VCF "##FILTER=<ID=FAIL,Description=\"FAIL\">\n";


	print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
	#print VCF "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	
	my %write_data;
	while (<CALL_DATA>) {
		chomp;
		my @mas = split/\t/;
		my $AD;
		my $RD;
		my @info = split/;/, $mas[2];
		foreach my $arg (@info) {
			if ($arg =~ /AODAD.*=(\d+),(\d+)/) {
				if ($arg =~ /_/) { # SKIPPING TAGS

					} else {
					$AD += $1;
					$RD += $2;
					}
				}
			}
		$AD = 0 unless defined $AD;
		$RD = 0 unless defined $RD;
		my ($chr, $pos, $ref, $alt);
		if ($mas[0] =~ /(\S+):(\d+)(\S+)>(\S+)/) {
			$chr = $1;
			$pos = $2;
			$ref = $3;
			$alt = $4;
			}
		$mas[1] = 0 if uc($mas[1]) eq 'NA';
		$write_data{$mas[0]} = {};
		$write_data{$mas[0]}->{chr} = $chr;
		$write_data{$mas[0]}->{pos} = $pos;
		$write_data{$mas[0]}->{ref} = $ref;
		$write_data{$mas[0]}->{alt} = $alt;
		$write_data{$mas[0]}->{qual} = $mas[1];
		$write_data{$mas[0]}->{filter} = '.';
		$write_data{$mas[0]}->{info} = "AODAD=$AD,$RD;$mas[2]";
		$write_data{$mas[0]}->{format} = "GT";
		$write_data{$mas[0]}->{sample} = ".";
		}
	
	foreach my $variant (sort {$write_data{$a}->{chr} cmp $write_data{$b}->{chr} || $write_data{$a}->{pos} <=> $write_data{$b}->{pos}} keys %write_data) {
		print VCF "",$write_data{$variant}->{chr};
		print VCF "\t",$write_data{$variant}->{pos};
		print VCF "\t.";
		print VCF "\t",$write_data{$variant}->{ref};
		print VCF "\t",$write_data{$variant}->{alt};
		print VCF "\t",$write_data{$variant}->{qual};
		print VCF "\t",$write_data{$variant}->{filter};
		print VCF "\t",$write_data{$variant}->{info};
		print VCF "\t",$write_data{$variant}->{format};
		print VCF "\t",$write_data{$variant}->{sample};
		print VCF "\n";
		}

	close CALL_DATA;
	close VCF;
	}

sub run {
	my $options = shift;
	create_vcf($options);
	}

sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'        => \$opts{'h'},
		'input|input=s'   => \$opts{'input'},
		'output|output=s' => \$opts{'output'},
		'sample|sample=s' =>\$opts{'sample'},
		'seed|seed=s'   => \$opts{'seed'}
	);
	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'input'});
	pod2usage(1) if(!$opts{'output'});
	pod2usage(1) if(!$opts{'sample'});
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

AOD_AMPL_CALL VCF CREATETOR

=head1 SYNOPSIS

CREATE VCF FILE FROM AOD_AMPL_CALL DATA

Options:

    -input   [REQUIRED] - input .cdata file
    -output  [REQUIRED] - output .vcf file
    -sample  [REQUIRED] - input .bam file




















