package Design;

use strict;
use warnings;
use Dir::Self;
use Segment;
use Sample;
use Allele;

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(max);
use Mojo::Base -base;

our @ISA = qw(Exporter);
our @EXPORT     = qw//;

sub new {
	my $class = shift;
	my $value = shift;
	my $self = {};
	$self->{controls}	= [];
	$self->{segments}	= {};
	$self->{amplicons}	= {};
	$self->{ampliconHash}	= {};
	$self->{seqdic}		= {};
	return (bless $self, $class);
	}

sub getSegmentKey {
	my $segment = shift;
	my $key = $segment->{contig} . ":" . $segment->{start} . "-" . $segment->{end};
	return $key;
	}

sub loadPanel {
	my $class	= shift;
	my $panelFile	= shift;
	$class->{segments}	= $class->define_segments($panelFile);
	my %hash = map { getSegmentKey($_) => $_} @{$class->segments};
	$class->{segments} = \%hash;
	$class->{amplicons}	= $class->define_amplicons($panelFile);
	$class->map_amplicons;
	die "Can not identify genome regions in input BED file ($panelFile)\n" unless defined $class->segments;
	die "Can not identify genome regions in input BED file ($panelFile)\n" unless defined $class->ampliconHash;
	}

sub segments {
	my $class	= shift;
	my $info	= shift;
	if (defined($info)) {
		my @segments = @{$class->{segments}};
		foreach my $key (keys %info) {
			@segments = grep {$_->{$key} eq $info->{$key}} @segments;
			}
		}
	return $class->{segments};
	}

sub amplicons {
	my $class	= shift;
	return $class->{amplicons};
	}

sub ampliconHash {
	my $class	= shift;
	return $class->{ampliconHash};
	}

sub loadSeqDic {
	my $class	= shift;
	my $seqdic	= shift;
	$class->{seqdic} = $seqdic;
	}

sub seqdic {
	my $class = shift;
	return $class->{seqdic};
	}

sub loadVarDic {
	my $class = shift;
	my $vcf = shift;
	my $header = shift;

	open (my $vcf_fh, "<$vcf");

	my $mutations = [];
	while (<$vcf_fh>) {
		chomp;
		next if m!#!;
		my @mas = split/\t/;
		next unless $mas[1] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in VCF file. Non-concordant with input BAM file\n" unless defined $header->{$mas[0]};
		die "Mutation position ($mas[1]) out of contig length for chromosome $mas[0]\n" if $mas[1] > $header->{$mas[0]};

		my $current = {};
		$current->{contig} = $mas[0];
		$current->{position} = $mas[1];
		$current->{ref} = $mas[3];
		$current->{alt} = $mas[4];
		$current->{reads} = [];
		$current->{name} = '.';
		foreach my $Segment ($class->segments)
		push (@{$mutations}, $current);
		}
	close $vcf_fh;

	if (defined($mutations)) {
		$mutations = [sort {$a->{position} <=> $b->{position}} @{$mutations}];
		}

	return $mutations;
	}

sub define_segments { # simple merge BED file
	my $class	= shift;
	my $inputFile   = shift;
	my $header      = $class->seqdic;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";
	
	my $segments;
	my $current_segment;
	while (<$panel_fh>) {
		my @mas = split/\t/;
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in BED file. Non-concordant with input BAM file\n" unless defined $header->{$mas[0]};
		die "Segment end position ($mas[2]) out of contig length for chromosome $mas[0]\n" if $mas[2] > $header->{$mas[0]};
		unless (defined $current_segment) {
			$current_segment = Segment->new({contig => $mas[0], start => $mas[1], end => $mas[2], mutations => []});
			next;
			}
		if (($mas[0] eq $current_segment->contig) and ($mas[1] < $current_segment->start)) {
			die "Input BED file is not sorted at position $mas[0]:$mas[1] <-; Use \"sort -k1,1 -k2,2n <input_bed>\" to sort bed file\n";
			}
		if (($mas[0] eq $current_segment->contig)
				and ($mas[1] > $current_segment->start)
				and ($mas[1] < $current_segment->end)) {
			if ($mas[2] > $current_segment->end) {
				$current_segment->{end} = $mas[2];
				} else {
				next;
				}
			} else {
			push (@{$segments}, $current_segment);
			$current_segment = Segment->new({contig => $mas[0], start => $mas[1], end => $mas[2], mutations => []});
			}
		}
	return $segments;
	}

sub define_amplicons {
	my $class	= shift;
	my $inputFile	= shift;
	my $amplicons;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";
	while (<$panel_fh>) {
		my @mas = split/\t/;
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		my $Amplicon = Segment->new({contig => $mas[0], start => $mas[1], end => $mas[2]});
		my $name = "$mas[0]:$mas[1]-$mas[2]";
		$amplicons->{$name} = $Amplicon;
		}
	close $panel_fh;
	return $amplicons;
	}

sub map_amplicons {
	my $class	= shift;
	my $header      = $class->seqdic;

	foreach my $Amplicon (keys %{$class->amplicons}) {
		$Amplicon = $class->amplicons->{$Amplicon};
		for (my $i = $Amplicon->start + 1; $i <= $Amplicon->end; $i++) {
			my $name = $Amplicon->contig . ":$i";
			$class->ampliconHash->{$name} = [] unless defined($class->ampliconHash->{$name});
			push @{$class->ampliconHash->{$name}}, $Amplicon;
			}
		}
	}













1;
