use strict;
use warnings;
use Data::Dumper;
use Getopt::Long 'GetOptions';
use Pod::Usage;
use Math::CDF;

sub get_content {
	my $input = shift;

	open (READ, "<$input");
	
	my %changes;
	my %depths;
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
					last if $count eq 0;
					if ($data[$j] =~ /[0-9]/) {
						$count = $count . $data[$j];
						next;
						}
					if ($data[$j] =~ /[ACGTNacgtn]/) {
						$alt = $alt . $data[$j];
						$count = $count - 1;
						next;
						}
					last;
					}
				$alt = uc($alt);
				if ($data[$i] eq '+') {
					$changes{"$mas[0]:$mas[1]$mas[2]>$mas[2]$alt"} += 1;
					$depths{"$mas[0]:$mas[1]$mas[2]>$mas[2]$alt"} = $mas[3];
					} else {
					$changes{"$mas[0]:$mas[1]$mas[2]$alt>$mas[2]"} += 1;
					$depths{"$mas[0]:$mas[1]$mas[2]$alt>$mas[2]"} = $mas[3];
					}
				$skip = length("$count$alt");
				}
			if ($data[$i] =~ /[ACGTNacgtn]/) {
				$data[$i] = uc($data[$i]);
				my $ref = $mas[2];
				$ref = uc($ref);
				$changes{"$mas[0]:$mas[1]$ref>$data[$i]"} += 1;
				$depths{"$mas[0]:$mas[1]$ref>$data[$i]"} = $mas[3];
				}
			}
		}
	
	close READ;
	my $data;
	$data->{changes} = \%changes;
	$data->{depth} = \%depths;

	my %filter;
	foreach my $key (keys %changes) {
		my $score = 1-Math::CDF::ppois($changes{$key}, $depths{$key}*0.0008);
		$filter{$key} = 'FAIL';
		$filter{$key} = 'PASS' if $score < 0.0001;
		$filter{$key} = 'FAIL' if $depths{$key} < 10;
		}
	$data->{filter} = \%filter;
	return $data;
	}

sub print_header {
	my $output = shift;
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
	$year = 1900 + $year;
	$mon  = (length($mon+1) eq 1 ? "0".($mon+1) : "".($mon+1)."");
	$mday = (length($mday) eq 1 ? "0$mday" : "$mday");
	
	print $output "##fileformat=VCFv4.1\n";
	print $output "##fileDate=$year$mon$mday\n";
	print $output "##source=mpileup\n";
	print $output "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Alternative allele count\">\n";
	print $output "##INFO=<ID=DP,Number=.,Type=String,Description=\"Depths\">\n";
	print $output "##FILTER=<ID=PASS,Description=\"PASS\">\n";
	print $output "##FILTER=<ID=FAIL,Description=\"FAIL\">\n";
	print $output "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
	}

sub c_1 {
        my $a = shift;
        my $b = shift;
        $a = ($a =~ /(\S+):(\d+)/ ? $2 : "NA");
        $b = ($b =~ /(\S+):(\d+)/ ? $2 : "NA");
        return ($a <=> $b)
        }

sub c_2 {
        my $a = shift;
        my $b = shift;
        $a = ($a =~ /(\S+):(\d+)/ ? $1 : "NA");
        $b = ($b =~ /(\S+):(\d+)/ ? $1 : "NA");
        return ($a cmp $b)
        }


sub print_content {
	my $data = shift;
	my $output = shift;
	my $limit = shift;
	my %changes = %{$data->{changes}};
	my %depth = %{$data->{depth}};
	my %filter = %{$data->{filter}};
	#foreach my $key (sort {$a cmp $b} keys %changes) {
	foreach my $key (sort {c_2($a, $b) || c_1($a, $b)} keys %changes) {
		next if $changes{$key} < $limit;
		if ($key =~ /(\S+):(\d+)(\S+)>(\S+)/) {
			print $output "$1\t$2\t.\t".uc($3)."\t".uc($4)."\t.\t",$filter{$key},"\tAD=$changes{$key};DP=$depth{$key}\n";
			}
		}
	}

sub run {
	my $options = shift;

	open (my $output, ">".$options->{output});
	my $data = get_content($options->{input});
	print_header($output);
	print_content($data, $output, $options->{limit});
	}

sub option_builder {
	my ($factory) = @_;
	my %opts;
	&GetOptions (
		'h|help'        => \$opts{'h'},
		'input|input=s'   => \$opts{'input'},
		'output|output=s' => \$opts{'output'},
		'limit|limit=s' =>\$opts{'limit'},
		'seed|seed=s'   => \$opts{'seed'}
	);
	pod2usage(0) if($opts{'h'});
	pod2usage(1) if(!$opts{'input'});
	pod2usage(1) if(!$opts{'output'});
	$opts{'limit'} = 4 unless defined $opts{'limit'};
	return \%opts;
	}

{
        my $options = option_builder();
	eval {run($options)};
	if ($@) {
		print STDERR "$@\n";
		} else {
		}
}


__END__

=head1 NAME

mpileupToVcf

=head1 SYNOPSIS

Create VCF file from Samtools mpileup file

Options:

    -input   [REQUIRED] - input mpileup file
    -output  [REQUIRED] - output .vcf file
    -limit   [DEFAULT:4] - minimum count of alternative alleles to output


