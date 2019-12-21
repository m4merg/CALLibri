package Sample;

use strict;
use warnings;
use Dir::Self;
use Allele;

use Data::Dumper;
use Storable 'dclone';
use List::Util qw(min max);
use Mojo::Base -base;
use Try::Tiny;                  #Install


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
	return $class->{allele}->{$alleleName};
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
        if (scalar @{$amplicons} eq 0) {
                return undef;
                } elsif ((scalar @{$amplicons}) eq 1) {
                my $name = $Variation->{contig} . ":" . $amplicons->[0]->[0] . "-" . $amplicons->[0]->[1];
                return $name;
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
        my $name = $Variation->{contig} . ":" . $selection->[0] . "-" . $selection->[1];
        return $name;
        }

sub pipeline {
        my $class	= shift;
        my $segment	= shift;
        my $sam = $class->sam;

        my $sam_segment = $sam->segment($segment->{contig}, $segment->{start}, $segment->{end});
        return undef unless defined $sam_segment;
        my @all_alignments = $sam_segment->features;
        foreach my $alignment (@all_alignments) {
                foreach my $CandidateVariation (@{$segment->{variations}}) {
                        if (defined($alignment->get_tag_values("SUPPLEMENTARY"))) {
                                next if $alignment->get_tag_values("SUPPLEMENTARY") eq '1';
                                }
                        next if $alignment->get_tag_values("UNMAPPED") eq '1';
                        next if $alignment->get_tag_values("NOT_PRIMARY") eq '1';

                        my $stat = get_stat($CandidateVariation, $alignment);
                        next unless defined($stat);
                        my $qscore = $class->get_qscore($alignment, $stat);
                        #print "!",$alignment->qname,"\t",$stat->{oref_add},"\t",$stat->{oalt_add},"\t",Score->new($qscore)->phred,"\n";
                        my $read;
                        $read->{name}           = $alignment->qname;
                        $read->{BQ}             = $qscore;
                        $read->{strand}         = $alignment->strand;
                        $read->{amplicon}       = select_amplicon($CandidateVariation, $alignment);
                        if (($stat->{oref} eq ($CandidateVariation->{ref})) and ($stat->{oalt} eq ($CandidateVariation->{alt}))) {
                                        $read->{vote} = 'alt';
                                        $class->allele($CandidateVariation->{index})->add_read($read);
                                } elsif (((length($stat->{oref}) eq length($CandidateVariation->{ref}))and
                                        (length($stat->{oalt}) eq length($CandidateVariation->{ref})))or
                                        (substr($stat->{match}, $stat->{aps}, $stat->{ape}-$stat->{aps}) =~ /^\|*$/)) {
                                        $read->{vote} = 'ref';
                                        $class->allele($CandidateVariation->{index})->add_read($read);
                                        } else {
                                                #print STDERR substr($stat->{match}, $stat->{aps}-2, $stat->{ape}-$stat->{aps} + 4),"\n";
                                                #print STDERR $alignment->qname,"\tWHAT?\n";
                                        }
                        }
                }
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

        my $ts = $ref; $ts =~ s/-//g;
        return undef if ($rpe) > length($ts);

        my $cigar = Cigar->new($alignment->cigar_str);
        try {($qps, $opqs) = $cigar->rpos_to_qpos($rps);};
        return undef unless defined $qps;
        $qpe = $qps + length($Mutation->{alt});
        return undef if $qpe < 0;
        return undef if $qps < 0;
        $aps = $rps + $cigar->get_shift($rps);
        $ape = $aps + max(length($Mutation->{alt}),length($Mutation->{ref}));
        my $stat;
        $stat->{oref} = substr($ref, $aps, $ape-$aps);
        $stat->{oalt} = substr($query, $aps, $ape-$aps);

        $stat->{oref_add} = substr($ref, $aps - 2, $ape-$aps + 4);
        $stat->{oalt_add} = substr($query, $aps - 2, $ape-$aps + 4);

        $stat->{oref} =~ s/-//g;
        $stat->{oalt} =~ s/-//g;

        $stat->{match} = $match;
        $stat->{aps} = $aps; $stat->{ape} = $ape;
        $stat->{qps} = $qps; $stat->{qpe} = $qpe;

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
