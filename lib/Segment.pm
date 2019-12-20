package Segment;

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
	my $class	= shift;
	my $self	= shift;
	return (bless $self, $class);
	}

sub start {
	my $class = shift;
	return $class->{start};
	}

sub end {
	my $class = shift;
	return $class->{end};
	}

sub contig {
	my $class = shift;
	return $class->{contig};
	}

sub mutations {
	my $class = shift;
	}












1;
