package Allele;

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
	my $self = {};
	$self->{reads} = [];
	return (bless $self, $class);
	}

sub add_read {
	my $class	= shift;
	my $read	= shift;
	push (@{$class->{reads}}, $read);
	}

sub readCount {
	my $class	= shift;
	my $info	= shift;
	my $count = 0;
	my $reads = dclone $class->{reads};
	my $opposite_reads = [];
	if (defined($info->{strand})) {
		$reads = [grep {$_->{strand} eq $info->{strand}} @{$reads}];
		}
	if (defined($info->{amplicon})) {
		$reads = [grep {$_->{amplicon} eq $info->{amplicon}} @{$reads}];
		}
	if (defined($info->{BQ})) {
		$reads = [grep {$_->{BQ} <= $info->{BQ}} @{$reads}];
		}
	if (defined($info->{vote})) {
		$opposite_reads = [grep {$_->{vote} ne $info->{vote}} @{$reads}];
		$reads = [grep {$_->{vote} eq $info->{vote}} @{$reads}];
		}
	if (scalar @{$reads} > 0) {
		$count += readCountCalc($reads, 'positive');
		}
	if (scalar(@{$opposite_reads}) > 0) {
		$count += readCountCalc($opposite_reads, 'negative');
		}
        return $count;
	}

sub readCountCalc {
	my $reads = shift;
	my $type  = shift; # Either 'positive' or 'negative'
	$type = 'positive' unless defined $type;
	my $sum;
	foreach my $Read (@{$reads}) {
		if ($type eq 'positive') {
			$sum += 1 - ($Read->{BQ});
			}
		if ($type eq 'negative') {
			$sum += $Read->{BQ}/3;
			}
		}
	return $sum;
	}













1;
