package Sample;

use strict;
use warnings;
use Dir::Self;
use parent 'Design';

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(max);
use Mojo::Base -base;

our @ISA = qw(Exporter);
our @EXPORT     = qw//;

sub new {
	my $class	= shift;
	my $bamFile	= shift;
	my $self = {};
	$self->{sam} = load_sam($bamFile);
	$self->{header} = load_header($self->{sam});
	return (bless $self, $class);
	}

sub load_header {
	my $sam = shift;
	my $references;
	for (my $i = 0; $i < $sam->header->n_targets; $i++) {
		$references->{$sam->header->target_name->[$i]} = $sam->header->target_len->[$i];
		}
	return $references;
	}

sub sam {
	my $class = shift;
	return $class->{sam};
	}

sub header {
	my $class = shift;
	return $class->{header};
	}

sub load_sam {
	my $inputBam = shift;
	my $sam = Bio::DB::Sam->new(
		-bam   => "$inputBam",
		-expand_flags  => 1,
		-autoindex => 1
		);
	return $sam;
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
