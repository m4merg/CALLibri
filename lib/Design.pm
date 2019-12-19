package Design;

use strict;
use warnings;
use Dir::Self;

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
	$self->{controls} = [];
	$self->{segments} = [];
	$self->{ampliconHash} = {};
	$self->{seqdic} = {};
	return (bless $self, $class);
	}

sub loadPanel {
	my $class	= shift;
	my $panelFile	= shift;
	$self->{segments}	= define_segments($panelFile);
	$self->{ampliconHash}	= define_amplicons($panelFile);
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

sub define_segments { # simple merge BED file
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
			$current_segment = {contig => $mas[0], start => $mas[1], end => $mas[2], mutations => []};
			next;
			}
		if (($mas[0] eq $current_segment->{contig}) and ($mas[1] < $current_segment->{start})) {
			die "Input BED file is not sorted at position $mas[0]:$mas[1] <-; Use \"sort -k1,1 -k2,2n <input_bed>\" to sort bed file\n";
			}
		if (($mas[0] eq $current_segment->{contig})
				and ($mas[1] > $current_segment->{start})
				and ($mas[1] < $current_segment->{end})) {
			if ($mas[2] > $current_segment->{end}) {
				$current_segment->{end} = $mas[2];
				} else {
				next;
				}
			} else {
			push (@{$segments}, $current_segment);
			$current_segment = {contig => $mas[0], start => $mas[1], end => $mas[2], mutations => []};
			}
		}
	return $segments;
	}

sub define_amplicons {
	my $inputFile   = shift;
	my $header      = $class->seqdic;
	open (my $panel_fh, "<$inputFile") or die "Can not open panel file\n";

	my $amplicon;
	while (<$panel_fh>) {
		my @mas = split/\t/;
		next unless defined $mas[2];
		next unless $mas[1] =~ /^\d+$/;
		next unless $mas[2] =~ /^\d+$/;
		die "Unknown contig name $mas[0] in BED file. Non-concordant with input BAM file\n" unless defined $header->{$mas[0]};
		die "Segment end position ($mas[2]) out of contig length for chromosome $mas[0]\n" if $mas[2] > $header->{$mas[0]};
		for (my $i = $mas[1] + 1; $i <= $mas[2]; $i++) {
			my $name = "$mas[0]:$i";
			$amplicon->{$name} = [] unless defined($amplicon->{$name});
			push @{$amplicon->{$name}}, "$mas[0]:$mas[1]-$mas[2]";
			}
		}

	close $panel_fh;
	return $amplicon;
	}













1;
