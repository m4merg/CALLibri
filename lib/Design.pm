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
	return (bless $self, $class);
	}

sub phred {
	my $class = shift;
	return $class->{phred}
	}

sub probability {
	my $class = shift;
	return $class->{prob};
	}

sub prob {
	my $class = shift;
	return $class->{prob};
	}

sub score {
	my $phred	= shift;
	my $score	= 10 ** (-$phred/10);
	return $score;
	}

sub unscore {
	my $score	= shift;
	return int(-10*log($score)/log(10))
	}















1;
