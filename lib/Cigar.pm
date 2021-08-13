package Cigar;

use strict;
use warnings;
use Dir::Self;
use parent 'Bio::Cigar';

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(max);

sub get_shift {
	my $class	= shift;
	my $n		= shift;
	my $cigar = $class->string;
	my $count_cigar	= 0;
	while ($n > 0) {
		if ( $cigar =~ /^(\d+)(S|M|I|D)/ ) {
			my $number = $1 - 1;
			my $letter = $2;
			if ( ( $letter eq "S" ) or ( $letter eq "I" ) ) {
				++$count_cigar;
				++$n;
				}
			if ( $number eq "0" ) {
				$cigar =~ s/^\d+\D//;
				} else {
					$cigar =~ s/^\d+/$number/;
					}
			} elsif ( $cigar =~ /^(\d+)(D)/m ) {
				$cigar =~ s/^\d+D//;
				next;
				}
		--$n;
		}
	return $count_cigar;
	}















1;
