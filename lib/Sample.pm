package Sample;

use strict;
use warnings;
use Dir::Self;
use Allele;
use Bio::DB::Sam;

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(min max);
use Mojo::Base -base;
use Try::Tiny;                  #Install
use Bio::Cigar;                 #Install
use Score;
use Cigar;
use List::Util qw/shuffle/;

my $BQ_threshold = 25;
my $string_edit_calc = __DIR__ . "/lev.py";

our @ISA = qw(Exporter);
our @EXPORT     = qw//;
my @knownTags = qw(ShortOverlap);

sub new {
	my $class	= shift;
	my $bamFile	= shift;
	my $self = {};
	$self->{bampath} = $bamFile;
	return (bless $self, $class);
	}

sub init {
	my $class = shift;
	$class->{sam} = load_sam($class->{bampath});
	$class->{header} = load_header($class->{sam});
	$class->{BQRange}->{0} = '0-5';
	$class->{BQRange}->{1} = '0-5';
	$class->{BQRange}->{2} = '0-5';
	$class->{BQRange}->{3} = '0-5';
	$class->{BQRange}->{4} = '0-5';
	$class->{BQRange}->{5} = '0-5';
	$class->{BQRange}->{6} = '6-10';
	$class->{BQRange}->{7} = '6-10';
	$class->{BQRange}->{8} = '6-10';
	$class->{BQRange}->{9} = '6-10';
	$class->{BQRange}->{10} = '6-10';
	$class->{BQRange}->{11} = '11-15';
	$class->{BQRange}->{12} = '11-15';
	$class->{BQRange}->{13} = '11-15';
	$class->{BQRange}->{14} = '11-15';
	$class->{BQRange}->{15} = '11-15';
	$class->{BQRange}->{16} = '16-20';
	$class->{BQRange}->{17} = '16-20';
	$class->{BQRange}->{18} = '16-20';
	$class->{BQRange}->{19} = '16-20';
	$class->{BQRange}->{20} = '16-20';
	$class->{BQRange}->{21} = '21-25';
	$class->{BQRange}->{22} = '21-25';
	$class->{BQRange}->{23} = '21-25';
	$class->{BQRange}->{24} = '21-25';
	$class->{BQRange}->{25} = '21-25';
	$class->{BQRange}->{26} = '26-30';
	$class->{BQRange}->{27} = '26-30';
	$class->{BQRange}->{28} = '26-30';
	$class->{BQRange}->{29} = '26-30';
	$class->{BQRange}->{30} = '26-30';
	$class->{BQRange}->{31} = '31-35';
	$class->{BQRange}->{32} = '31-35';
	$class->{BQRange}->{33} = '31-35';
	$class->{BQRange}->{34} = '31-35';
	$class->{BQRange}->{35} = '31-35';
	foreach my $BQ (36 .. 100) {
		$class->{BQRange}->{$BQ} = '36-100';
		}
	


	}

sub get_value {
   my ( $self, $key ) = @_;
   return $self->{$key} // 0;
}

sub set_value {
   my ( $self, $key, $value ) = @_;
   $self->{$key} = $value;
}

	
sub allele {
	my $class	= shift;
	my $alleleName	= shift;
	if (defined($class->{allele})) {
		if (defined($class->{allele}->{$alleleName})) {
			return $class->{allele}->{$alleleName};
			} else {
			$class->{allele}->{$alleleName} = Allele->new($class);
			}
		} else {
		$class->{allele} = {};
		$class->{allele}->{$alleleName} = Allele->new($class);
		}
	$class->{allele}->{$alleleName}->{Sample} = $class;
	$class->{allele}->{$alleleName}->{alleleName} = $alleleName;
	if ($alleleName =~ /(\S+):(\d+)(\S+)>(\S+)/) {
		$class->{allele}->{$alleleName}->{refA} = $3;
		$class->{allele}->{$alleleName}->{altA} = $4;
		$class->{allele}->{$alleleName}->{edit_ops_forw} = edit_ops($3, $4);
		$class->{allele}->{$alleleName}->{edit_ops_rev} = edit_ops($4, $3);
		}
	return $class->{allele}->{$alleleName};
	}

sub edit_ops {
	my $str1 = shift;
	my $str2 = shift;
	my $edit_ops = `python $string_edit_calc $str1 $str2`;
	chomp $edit_ops;
	$edit_ops = [split/\n/, $edit_ops];
	my $result = [];
	map {if ($_ =~ /(\S+)>(\S+)/) {push @{$result}, [$1, $2]}} @{$edit_ops};
	return $result;
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

sub select_amplicon {
        my $Variation           = shift;
        my $alignment           = shift;
        my $amplicons = $Variation->{amplicons};
	return (undef, undef) unless defined $amplicons;
	my ($ref, $match, $query) = $alignment->padded_alignment;
	
	$match =~ s/^\s+|\s+$//g;
        if (scalar @{$amplicons} eq 0) {
                return undef;
                }
        my $min = 3000000000;
        my $selection;
        foreach my $Amplicon (@{$amplicons}) {
                my $start       = abs($alignment->start - $Amplicon->[0]);
                my $end         = abs($alignment->end - $Amplicon->[1]);
                if ($start + $end < $min) {
                        $min = $start + $end;
                        $selection = $Amplicon;
                        }
                }
	my $overlap = min($alignment->end, $selection->[1]) - max($alignment->start, $selection->[0]);
	my $distance_from_insert = min(
				abs(min($alignment->end, $selection->[1]) - $selection->[1]),
				abs(max($alignment->start, $selection->[0]) - $selection->[0])
				);
	my $seq_length = length($alignment->query->dna);
	#print STDERR "OVERLAP - $overlap\t$seq_length\t$distance_from_insert\n";
	my @tags;
	if ($overlap - $distance_from_insert < 40) {
		push @tags, "ShortOverlap";
		}
	return undef if $overlap < 0.5*length($match);
	#print STDERR "PASS\n";
        my $name = $Variation->{contig} . ":" . $selection->[0] . "-" . $selection->[1];
        return ($name, [@tags]);
        }

sub normalizeBQ {
	my $class	= shift;
	my $sam = $class->sam;
	my $Design = $class->Design;
	my $alignmentCountMax = 500000;
	my $NCount;
	my $BQMatrix;
	my $sumScore = 0;
	BQOUTER: foreach my $segment (shuffle @{$class->Design->segments}) {
		my $sam_segment = $sam->segment($segment->{contig}, $segment->{start}, $segment->{end});
		my @all_alignments = $sam_segment->features;
		foreach my $alignment (@all_alignments) {
			if (defined($alignment->get_tag_values("SUPPLEMENTARY"))) {
				next if $alignment->get_tag_values("SUPPLEMENTARY") eq '1';
				}
			next if $alignment->get_tag_values("UNMAPPED") eq '1';
			next if $alignment->get_tag_values("NOT_PRIMARY") eq '1';
			#print STDERR "",$alignment->qname,"\n";
			#next unless $alignment->qname eq 'KOVMX:09610:14033';
			my ($ref, $match, $query) = $alignment->padded_alignment;
			$ref =~ s/(-+)$/"_" x length($1)/e;
			$ref =~ s/^(-+)/"_" x length($1)/e;
			$query =~ s/(-+)$/"_" x length($1)/e;
			$query =~ s/^(-+)/"_" x length($1)/e;
			$ref = [split//, $ref];
			$match = [split//, $match];
			$query = [split//,$query];
			my @scores = @{$alignment->qscore};
			my $scoresPos = 0;
			for (my $i = 0; $i < scalar (@{$match}); $i++) {
				#print STDERR "",$ref->[$i],"",$match->[$i],"",$query->[$i]," ",$scores[$scoresPos],"\n";
				#next if $ref->[$i] eq '-';
				#next if $query->[$i] eq '-';
				next if $ref->[$i] eq 'N';
				next if $query->[$i] eq 'N';
				next if $ref->[$i] eq '_';
				next if $query->[$i] eq '_';
				$NCount->{$ref->[$i]} += 1;
				next if $ref->[$i] eq $query->[$i];
				$BQMatrix->{"".$ref->[$i].">".$query->[$i].""} += 1;
				$scoresPos += 1 unless $query->[$i] eq '-';
				}
			$alignmentCountMax -= 1;
			#print STDERR "$alignmentCountMax\n";
			last BQOUTER if $alignmentCountMax < 0;
			#if ($alignmentCountMax < 0) {
			#	foreach my $key (keys %{$BQMatrix}) {
			#		$alignmentCountMax += 10000 if $BQMatrix->{$key} < 30;
			#		last if $BQMatrix->{$key} < 30;
			#		}
			#	last BQOUTER if $alignmentCountMax < 0;
			#	}
			}
		}
	#print STDERR Dumper $NCount;
	#print STDERR Dumper $BQMatrix;
	$class->{BQMatrix} = $BQMatrix;
	$class->{NCount} = $NCount;
	}

sub pipeline {
	my $class	= shift;
	my $segment	= shift;
	my $sam = $class->sam;
	
	my $sam_segment = $sam->segment($segment->{contig}, $segment->{start}, $segment->{end});
	return undef unless defined $sam_segment;
	my @all_alignments = $sam_segment->features;
	#print STDERR "STARTED SEGMENT: ", $segment->{contig},"\t", $segment->{start},"\t", $segment->{end},"\n";
	#print STDERR "COUNT: ",scalar(@{$segment->{variations}}),"\n";
	#print STDERR "COUNT: ",scalar(@all_alignments),"\n";
	foreach my $alignment (@all_alignments) {
		foreach my $CandidateVariation (@{$segment->{variations}}) {
			my $index = $CandidateVariation->{index};
			#next if $index ne "chr9:135974142C>CG";
			if (defined($alignment->get_tag_values("SUPPLEMENTARY"))) {
				next if $alignment->get_tag_values("SUPPLEMENTARY") eq '1';
				}
			next if $alignment->get_tag_values("UNMAPPED") eq '1';
			next if $alignment->get_tag_values("NOT_PRIMARY") eq '1';
			
			#next unless $alignment->qname eq '53QOO:02914:04180';
			my $stat = get_stat($CandidateVariation, $alignment);
			#print STDERR Dumper $stat;
			next unless defined($stat);
			my $qscore = $class->get_qscore($alignment, $stat);
			#print STDERR "!",$alignment->qname,"\t",$stat->{oref_add},"\t",$stat->{oalt_add},"\t",$alignment->qual,"\t",Score->new($qscore)->phred,"\n";
			next if (Score->new($qscore)->phred) < $BQ_threshold;
			my $read;my $tags = [];
			#$stat->{name}           = $alignment->qname;
			$read->{BQ}			= $qscore;
			$read->{strand}			= $alignment->strand;
			($read->{amplicon}, $tags)	= select_amplicon($CandidateVariation, $alignment);
			#print "",join("\t", @{$tags}),"\n";
			#print "",($tags ? join("\t", @{$tags}) : "NA"),"\n";
			next unless defined ($read->{amplicon});
			#print STDERR "  -------   ",$CandidateVariation->{ref},"\n";
			#print STDERR "  -------   ",$CandidateVariation->{alt},"\n";
			#!T2QOU:03477:01058
			#print STDERR Dumper $CandidateVariation;
			#print STDERR "",substr($stat->{match}, $stat->{aps}, $stat->{ape}-$stat->{aps}),"\n";
			if ($tags) {
				$read->{tags} = {};
				foreach my $key (@{$tags}) {
					$read->{tags}->{$key} = 1;
					}
				} else{
				$read->{tags} = {};
				}
			foreach my $key (@knownTags) {
				if (not(defined($read->{tags}->{$key}))) {
					$read->{tags}->{$key} = 0;
					}
				}
			if (($stat->{oref} eq ($CandidateVariation->{ref})) and ($stat->{oalt} eq ($CandidateVariation->{alt})) and ($stat->{match_around} =~ /^\|*-\|*$/)) {
					$read->{vote} = 'alt';
					$read->{MQ} = $alignment->qual;
					$class->allele($CandidateVariation->{index})->add_read($read);
					#print STDERR "ALT\n";
					#print STDERR Dumper $stat;
				} elsif ($stat->{match_var} =~ /^\|*$/) {
					$read->{vote} = 'ref';
					$read->{MQ} = $alignment->qual;
					#print STDERR "REF\n";
					#print STDERR Dumper $stat;
					$class->allele($CandidateVariation->{index})->add_read($read);
					} else {
						#print STDERR Dumper $stat;
						#print STDERR "WHAT?\n";
						#print STDERR substr($stat->{ref}, $stat->{aps}-2, $stat->{ape}-$stat->{aps} + 4),"\n";
						#print STDERR substr($stat->{match}, $stat->{aps}-2, $stat->{ape}-$stat->{aps} + 4),"\n";
						#print STDERR substr($stat->{query}, $stat->{aps}-2, $stat->{ape}-$stat->{aps} + 4),"\n";
						#print STDERR $alignment->qname,"\tWHAT?\n";
					}
			}
		}
	}

sub prep {
	my $class	= shift;
	my $segment	= shift;
	my $sam = $class->sam;

	my $sam_segment = $sam->segment($segment->{contig}, $segment->{start}, $segment->{end});
	return undef unless defined $sam_segment;
	my @all_alignments = $sam_segment->features;
	foreach my $alignment (@all_alignments) {
		my ($ref,$matches,$query) = $alignment->padded_alignment;
		parse_matches($alignment->start, $alignment->padded_alignment);
		}
	}

sub parse_matches {
	my $position	= shift;
	my @ref		= split//, shift;
	my @matches	= split//, shift;
	my @alt		= split//, shift;
	my $i = 0;
	#print STDERR "$position\t@ref\t@matches\t@alt\n";exit;
	while ($i < scalar @matches) {
		my $vRef = '';
		my $vAlt = '';
		if ($matches[$i] eq ' ') {
			my $j = 0;
			while ($matches[$i + $j] eq ' ') {
				last unless defined $ref[$i + $j];
				last unless defined $alt[$i + $j];
				$vRef = "$vRef" . $ref[$i + $j];
				$vAlt = "$vAlt" . $alt[$i + $j];
				++$j;
				}
			}
		while (1) {
			++$position if $ref[$i] ne '-';
			++$i;
			last unless defined $matches[$i];
			last if $matches[$i] eq '|';
			}
		}
	exit;
	}

sub get_stat { # see pipeline function
        my $Mutation = shift;
        my $alignment = shift;
        #aps,ape,rps,rpe,qps,qpe - 0-based coordinates
        #aps/ape - start/end position of variant site in alignment
        #rps/rpe - start/end positions of variation site in reference
        my $qps; my $rps; my $aps; my $qpe; my $rpe; my $ape;
        my $opqs; my $opqe;
        $rps = ($Mutation->{position} - $alignment->start);
        $rpe = ($Mutation->{position} - $alignment->start) + length($Mutation->{ref});
        return undef if $rps < 1;
        return undef if ($rps) > ($rpe);

        my ($ref, $match, $query) = $alignment->padded_alignment;
	
	#print STDERR "$ref\n$match\n$query\n";
        my $ts = $ref; $ts =~ s/-//g;
        return undef if ($rpe) > length($ts);

        my $cigar = Cigar->new($alignment->cigar_str);
        try {($qps, $opqs) = $cigar->rpos_to_qpos($rps);};
        return undef unless defined $qps;
        $qpe = $qps + length($Mutation->{alt});
        return undef if $qpe < 0;
        return undef if $qps < 0;
	#print STDERR "",$alignment->cigar_str,"\n$rps\n",$cigar->get_shift($rps),"\n";
        $aps = $rps + $cigar->get_shift($rps);
        $ape = $aps + max(length($Mutation->{alt}),length($Mutation->{ref}));
        my $stat;
        $stat->{oref} = substr($ref, $aps, $ape-$aps);
        $stat->{oalt} = substr($query, $aps, $ape-$aps);
	
	#print STDERR " !--------------------- $aps, $ape\n";
        $stat->{oref_add} = substr($ref, $aps - 2, $ape-$aps + 4);
        $stat->{oalt_add} = substr($query, $aps - 2, $ape-$aps + 4);

        $stat->{oref} =~ s/-//g;
        $stat->{oalt} =~ s/-//g;

        $stat->{match} = $match;
	$stat->{match_var} = substr($stat->{match}, $aps, $ape-$aps);
	$stat->{ref} = $ref;
	$stat->{query} = $query;
        $stat->{aps} = $aps; $stat->{ape} = $ape;
        $stat->{qps} = $qps; $stat->{qpe} = $qpe;
	$stat->{match_around} = ($aps < 1 ? substr($match, 0, $aps) : substr($match, $aps - 1, 1)).'-'.($ape < length($match) ? substr($match, $ape, 1) : "");

        return $stat;
        }

	
sub get_qscore { # see pipeline function
	my $class	= shift;
        my $alignment   = shift;
        my $stat        = shift;
        my @scores = @{$alignment->qscore};
        my $qscore_start = max(0, $stat->{qps} - $class->Design->config->{qscore_averaging_range});
        my $qscore_end   = min($stat->{qpe} - 1 + $class->Design->config->{qscore_averaging_range}, (scalar @scores) - 1);
        my $qscore = 0;
        map {$qscore = $qscore + Score->new($_)->prob} @scores[($qscore_start)..($qscore_end)];
        $qscore = $qscore / ($qscore_end - $qscore_start + 1);
        return $qscore;
        }

sub Design {
	my $class = shift;
	return $class->{Design};
	}












1;
