use strict;
use warnings;

my $b = 56/100000000000;
while ($b < 0.000001) {
	$b = $b * 10
	}
print "$b\n";
