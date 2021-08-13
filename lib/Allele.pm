package Allele;

use strict;
use warnings;
use Dir::Self;

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(max sum);
use Mojo::Base -base;

my $string_edit_calc = __DIR__ . "/lev.py";
our @ISA = qw(Exporter);
our @EXPORT     = qw//;

sub new {
	my $class = shift;
	my $self = {};
	$self->{reads} = [];
	return (bless $self, $class);
	}

sub Sample {
	my $class = shift;
	return $class->{Sample};
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
	#print STDERR "CLASS_HERE1\n";
	#print STDERR Dumper $class;
	#exit;
	#print STDERR "CLASS_HERE2\n";
	#print STDERR Dumper $info;
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
		$count += $class->readCountCalc($reads, 'positive');
		}
	if (scalar(@{$opposite_reads}) > 0) {
		$count += $class->readCountCalc($opposite_reads, 'negative');
		}
        return $count;
	}

sub readCountCalc {
	my $class = shift;
	my $reads = shift;
	my $type  = shift; # Either 'positive' or 'negative'
	$type = 'positive' unless defined $type;
	my $sum;
	foreach my $Read (@{$reads}) {
		#print STDERR Dumper $Read;
		#print STDERR "$type\n";
		if ($type eq 'positive') {
			#$sum += 1 - ($Read->{BQ});
			#print STDERR (1 - $class->error_prob("edit_ops_forw")),"\n";
			#print STDERR (1 - $Read->{BQ}),"\n";
			$sum += 1 - $class->error_prob("edit_ops_forw");
			}
		if ($type eq 'negative') {
			#$sum += $Read->{BQ}/3;
			#print STDERR $class->error_prob("edit_ops_rev"),"\n";
			#print STDERR $Read->{BQ}/3,"\n";
			$sum += $class->error_prob("edit_ops_rev");
			}
		}
	return $sum;
	}

sub error_prob {
	my $class = shift;
	my $type  = shift;
	my $prob  = 1;
	foreach my $OP (@{$class->{$type}}) {
		my $nominator = $class->Sample->{BQMatrix}->{$OP->[0].">".$OP->[1]};
		my $denominator;
		if (($OP->[0] eq '-')or($OP->[1] eq '-')) {
			$denominator = sum(@{$class->Sample->{NCount}})
			} else {
			$denominator = $class->Sample->{NCount}->{$OP->[0]};
			}
		#print STDERR "$nominator\t$denominator\n";
		$prob = $prob * $nominator / $denominator;
		}
	return $prob;
	}
























1;
