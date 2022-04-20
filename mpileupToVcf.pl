use strict;
use warnings;

open (READ, "<$ARGV[0]");

my %changes;
while (<READ>) {
	chomp;
	my @mas = split/\t/;
	my @data = split//, $mas[4];
	my $skip = 0;
	for (my $i = 0; $i < scalar @data; $i++) {
		if ($skip > 0) {
			$skip = $skip - 1;
			next;
			}
		next if $data[$i] eq '.';
		next if $data[$i] eq ',';
		next if $data[$i] eq '>';
		next if $data[$i] eq '<';
		if ($data[$i] eq '^') {
			$skip = 1;
			next;
			}
		if (($data[$i] eq '+')or($data[$i] eq '-')) {
			my $count = '';
			my $alt = '';
			for (my $j = $i + 1; $j < scalar @data; $j++){
				if ($data[$j] =~ /[0-9]/) {
					$count = $count . $data[$j];
					next;
					}
				if ($data[$j] =~ /[ACGTNacgtn]/) {
					$alt = $alt . $data[$j];
					next;
					}
				last;
				}
			$alt = uc($alt);
			if ($data[$i] eq '+') {
				$changes{"$mas[0]:$mas[1]$mas[2]>$mas[2]$alt"} += 1;
				} else {
				$changes{"$mas[0]:$mas[1]$mas[2]$alt>$mas[2]"} += 1;
				}
			$skip = length("$count$alt");
			}
		if ($data[$i] =~ /[ACGTNacgtn]/) {
			$data[$i] = uc($data[$i]);
			my $ref = $mas[2];
			$ref = uc($ref);
			$changes{"$mas[0]:$mas[1]$ref>$data[$i]"} += 1;
			}
		}
	}

close READ;

foreach my $key (sort {$a cmp $b} keys %changes) {
	next if $changes{$key} < 4;
	if ($key =~ /(\S+):(\d+)(\S+)>(\S+)/) {
		print "$1\t$2\t.\t$3\t$4\t$changes{$key}\n";
		}
	}

