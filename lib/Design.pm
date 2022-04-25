package Design;

use strict;
use warnings;
use Dir::Self;

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(max);
use Mojo::Base -base;
use Sample;

our @ISA = qw(Exporter);
our @EXPORT     = qw//;

sub new {
	my $class = shift;
	my $self = {};
	$self->{controls} = [];
	$self->{config} = {};
	return (bless $self, $class);
	}

sub config {
	my $class = shift;
	return $class->{config};
	}	

sub init {
	my $class	= shift;
	my $info	= shift;
	$class->load_seqDic($info->{seqdic});
	$class->define_segments($info->{BED});
	$class->load_VarDict($info->{VCF});
	$class->load_amplicons($info->{BED});
	}

sub segments {
	my $class = shift;
	return $class->{segments};
	}

sub load_seqDic {
	my $class = shift;
	my $seq_dic = shift;
	$class->{seqdic} = $seq_dic;
	}

sub seqdic {
	my $class = shift;
	return $class->{seqdic};
	}

sub load_VarDict {
	my $class	= shift;
	my $vcf		= shift;

	open (my $vcf_fh, "<$vcf");
	
	while (<$vcf_fh>) {
		chomp;
		next if m!#!;
		my @mas = split/\t/;
		next unless $mas[1] =~ /^\d+$/;
		
		die "Unknown contig name $mas[0] in VCF file. Non-concordant with input BAM file\n" unless defined $class->seqdic->{$mas[0]};
		die "Mutation position ($mas[1]) out of contig length for chromosome $mas[0]\n" if $mas[1] > ($class->seqdic->{$mas[0]});

		my @altAlleles = split/,/, $mas[4];
		foreach my $alt (@altAlleles) {
			my $CandidateVariation = {};
			$CandidateVariation->{contig} = $mas[0];
			$CandidateVariation->{position} = $mas[1];
			$CandidateVariation->{ref} = $mas[3];
			$CandidateVariation->{alt} = $alt;
			$CandidateVariation->{name} = $mas[2];
			$CandidateVariation->{index} = "$mas[0]:$mas[1]$mas[3]>$mas[4]";
			my $count = 0;
			my @added;
			foreach my $seg (@{$class->segments}) {
				next if $CandidateVariation->{contig} ne $seg->{contig};
				if (($CandidateVariation->{position} > $seg->{start})and($CandidateVariation->{position} <= $seg->{end})) {
					push(@{$seg->{variations}}, $CandidateVariation);
					push (@added, $seg);
					++$count;
					}
				}
			die "Multiple maps to segments for mutation $mas[0]:$mas[1]$mas[3]>$alt\n" if $count > 1;
			#warn "Mutation $mas[0]:$mas[1]$mas[3]>$alt falls out from designed amplicons - it will be ignored\n" if $count eq 0;
			}
		}
	close $vcf_fh;
	}

sub define_segments { # simple merge BED file
	my $class = shift;
	my $inputFile   = shift;
	my $mergeDist	= shift;
	$mergeDist = 0 unless defined $mergeDist;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";

	my $current_segment;
	while (<$panel_fh>) {
		my @mas = split/\t/;
		#print STDERR "@mas\n" if $mas[1] eq '29130328';
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in BED file. Non-concordant with input BAM file\n" unless defined $class->seqdic->{$mas[0]};
		die "Segment end position ($mas[2]) out of contig length for chromosome $mas[0]\n" if $mas[2] > ($class->seqdic->{$mas[0]});
		unless (defined $current_segment) {
			$current_segment = {contig => $mas[0], start => $mas[1], end => $mas[2], variations => []};
			next;
			}
		if (($mas[0] eq $current_segment->{contig}) and ($mas[1] < $current_segment->{start})) {
			die "Input BED file is not sorted at position $mas[0]:$mas[1] <-; Use \"sort -k1,1 -k2,2n <input_bed>\" to sort bed file\n";
			}
		if (($mas[0] eq $current_segment->{contig})
				and ($mas[1] > $current_segment->{start})
				and ($mas[1] < $current_segment->{end} + $mergeDist)) {
			if ($mas[2] > $current_segment->{end}) {
				$current_segment->{end} = $mas[2];
				} else {
				next;
				}
			} else {
			#print STDERR Dumper $current_segment if $mas[1] eq '29130328';
			push (@{$class->{segments}}, $current_segment);
			$current_segment = {contig => $mas[0], start => $mas[1], end => $mas[2], variations => []};
			}
		}
	push (@{$class->{segments}}, $current_segment);
	
	}

sub load_amplicons {
	my $class	= shift;
	my $inputFile	= shift;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";

	while (<$panel_fh>) {
		my @mas = split/\t/;
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		foreach my $seg (@{$class->segments}) {
			next if scalar (@{$seg->{variations}}) eq 0;
			next if $seg->{contig} ne $mas[0];
			foreach my $CandidateVariation (@{$seg->{variations}}) {
				if (($CandidateVariation->{position} > $mas[1])and($CandidateVariation->{position} <= $mas[2])) {
					push(@{$CandidateVariation->{amplicons}}, [$mas[1], $mas[2]]);
					}
				}
			}
		}

	close $panel_fh;
	}

sub newSample {
	my $class	= shift;
	my $bam		= shift;
	my $Sample = Sample->new($bam);
	$Sample->{Design} = $class;
	push (@{$class->{samples}}, $Sample);
	return $Sample;
	}

sub newControl {
	my $class	= shift;
	my $bam		= shift;
	my $Sample = $class->newSample($bam);
	push (@{$class->{controls}}, $Sample);
	$Sample->set_value('type', 'control');
	return $Sample;
	}

sub controls {
	my $class = shift;
	return $class->{controls};
	}

sub samples {
	my $class = shift;
	return $class->{samples};
	}

sub Sample {
	my $class = shift;
	my $name = shift;
	foreach my $sample (@{$class->samples}) {
		next unless defined $sample->{name};
		next if $sample->{name} ne $name;
		return $sample;
		}
	return undef;
	}








1;
