package PPH::Profile;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Profile

=head1 DESCRIPTION

Calculate a PSIC matrix from a sequence alignment and also provide
routines to return PSIC scores for a given substitution

(former pph_psic)

=head1 SUBVERSION

 $LastChangedDate: 2012-02-16 22:27:02 -0500 (Thu, 16 Feb 2012) $
 $LastChangedRevision: 391 $
 $LastChangedBy: ivan $


=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
use Carp qw(cluck confess);
use strict;
use warnings;

use Fcntl qw(:flock);
use File::Basename;
use File::Copy;

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION = 0.1;
    @ISA           = qw(Exporter);
    @EXPORT        = qw();
    @EXPORT_OK     = qw();
    %EXPORT_TAGS = ();
    #----------------------------------------
    # DEBUGGING
    #----------------------------------------
    use PPH::Config; # imports: $DEBUG @DEBUGGING_MODE %CONFIG
    if ($DEBUG) {
      # Test if we have Smart::Comments installed
      eval { require Smart::Comments; Smart::Comments->import(-ENV) };
      die "Please install Smart::Comments module or switch off verbose/DEBUG mode\n" if $@;
    }
    #----------------------------------------
}
our @EXPORT_OK;

use PPH::BLAST;
use PPH::Misc;
use PPH::Align;

use constant PSICAA => 'ARNDCQEGHILKMFPSTWYV';
my %AAidx; @AAidx{split //, PSICAA} = 0 .. length(PSICAA) - 1;

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

=head1 ROUTINES

#----------------------------------------------------------------------

=head2 set_PSIC_profile

 Usage   : PPH::Profile->set_PSIC_profile($PROT, $SNP)
 Function: Calculate PSIC scores for the given SNP
           Calculates and loads PSIC matrix if no profile was found
 Args    : Object of type PROT and SNP
 Returns : PSIC profiles are stored as:
           PPH->{Protein}{Pkey}{Profile}{Msa} - protein MSA-based profile
           PPH->{Protein}{Pkey}{Profile}{Mz}  - MultiZ46Way exome-based
                                                profile (optional)

=cut

#----------------------------------------------------------------------
sub set_PSIC_profile {
  my ($PROT, $SNP) = @_;
  #### PPH__Profile__set_PSIC_profile: $PROT, $SNP

  # Attempt to calculate and load protein MSA-based PSIC profile matrix
  # into memory, if it has not been tried already.
  unless (exists $PROT->{Profile}{Msa}) {

    my %PRF;

    # The empty hashref stub indicating profile creation was attempted.
    # Replaced with ref to valid data if processing completes successfully.
    $PROT->{Profile}{Msa} = {};

    my $alnfile = $PROT->{Acc} . '.aln';
    my $prffile = $PROT->{Acc} . '.prf';

    # All calculated files are saved under SCRATCH directory by default. When SAVEKNOWN
    # option is set however, files for known UniProtKB proteins are saved to PRECOMPATH.
    $CONFIG{'SAVEPATH'} = $CONFIG{'SAVEKNOWN'} && $PROT->{UniProt} ? $CONFIG{'PRECOMPATH'} : $CONFIG{'SCRATCH'};
    # Note that we have to lock the whole pipeline BEFORE checking for the files availability
    # to avoid creating a race condition (well, as much as we possibly can).
    my $lockfile = "$CONFIG{SAVEPATH}/lock/$PROT->{Acc}.lock";
    open my $flk, '>', $lockfile or die "ERROR: Can't create lockfile $lockfile: $!\n";
    flock($flk, LOCK_EX) or die "ERROR: Failed to obtain lock for $lockfile: $!\n";
    ####### |     file locked: $lockfile

    my $aln_found = PPH::Misc::locate_file($alnfile, "$CONFIG{PRECOMPATH}/alignments", "$CONFIG{SCRATCH}/alignments");
    my $prf_found = PPH::Misc::locate_file($prffile, "$CONFIG{PRECOMPATH}/profiles", "$CONFIG{SCRATCH}/profiles");
    # Skip MSA pipeline completely if and only if both .aln and .prf files
    # are already present and accounted for
    if ($aln_found && $prf_found) {
      $PRF{Files}{Aln} = $alnfile;
      $PRF{Files}{Prf} = $prffile;
      my @source = qw{ Missing Precomputed Cached };
      push @{ $PRF{$source[$aln_found]} }, 'alignment';
      push @{ $PRF{$source[$prf_found]} }, 'profile';
    } else {
      # Profile calculation is not allowed, give up
      unless ($CONFIG{'PRFPRECOMPONLY'}) {
        # Do the whole shebang: BLAST / MSA / PSIC
        calculate_profile($PROT, \%PRF);
      }
    }
    close($flk);

    # Load complete matrix of profile scores into memory, preferrably from a local
    # file copy, if it is was just computed and saved to the current temp directory.
    $PRF{Scores} = _read_profile($prffile) if
      exists $PRF{Files} && exists $PRF{Files}{Prf} && defined $PRF{Files}{Prf};

    # Processing successfully completed
    $PROT->{Profile}{Msa} = \%PRF if exists $PRF{Scores} && @{ $PRF{Scores} };

  } # end of MSA-based profile creation / loading

  # Now extract the substitution's MSA PSIC scores and load them into SNP object.
  # Note: absence of profile data is not necesseraly a fatal exception
  # (but we should probably provide better diagnostics for such cases).
  my $PRF = $PROT->{Profile}{Msa};
  if (keys %$PRF && @{ $PRF->{Scores} }) {
    $SNP->{Scores}{Msa}{Psic1} = $PRF->{Scores}[$SNP->{Pos}-1][$AAidx{$SNP->{Aa1}}];
    $SNP->{Scores}{Msa}{Psic2} = $PRF->{Scores}[$SNP->{Pos}-1][$AAidx{$SNP->{Aa2}}];
    $SNP->{Scores}{Msa}{Nobs}  = $PRF->{Scores}[$SNP->{Pos}-1][-1];
    #----------------------------------------
    # IAA: Changed PsicD to signed difference of scores instead of abs() difference used before.
    #      Note: PPH::Rules::test() method still uses abs(PsicD) value for compatibility with PPHv1.
    #----------------------------------------
    $SNP->{Scores}{Msa}{PsicD} = $SNP->{Scores}{Msa}{Nobs} ?
        sprintf('%+.3f', $SNP->{Scores}{Msa}{Psic1} - $SNP->{Scores}{Msa}{Psic2}) : undef;
  }

  # Attempt to read precomputed MultiZ46Way exome-based PSIC profile
  # matrix for the SNP into memory, if it has not been tried already.
  # Note: there could be different transcript exons mapped to different
  #       SNP positions, therefore several profiles may be loaded.
  if ($CONFIG{'MULTIZPATH'} && exists $SNP->{Gene} &&
      exists $SNP->{Gene}{TxName} && defined $SNP->{Gene}{TxName} && length $SNP->{Gene}{TxName}) {
    my $txname = $SNP->{Gene}{TxName};
    unless (exists $PROT->{Profile}{Mz}{$txname}) {
      my %PRF;

      # The empty hashref stub indicating profile reading was attempted.
      # Replaced with ref to valid data if processing completes successfully.
      $PROT->{Profile}{Mz}{$txname} = {};

      my $alnfile = $SNP->{Gene}{TxName} . '.aln';
      my $prffile = $SNP->{Gene}{TxName} . '.prf';
      # Both precomputed aln and prf MZ files should always be present
      if (PPH::Misc::locate_file($alnfile, "$CONFIG{MULTIZPATH}/precomputed/aln") &&
          PPH::Misc::locate_file($prffile, "$CONFIG{MULTIZPATH}/precomputed/prf")) {
        $PRF{Files}{Aln} = $alnfile;
        $PRF{Files}{Prf} = $prffile;
        push @{ $PRF{Precomputed} }, 'alignment', 'profile';
        # Load complete profile matrix into memory
        $PRF{Scores} = _read_profile($prffile);
        # Extract CDS protein sequence for the transcript, currently only required for
        # ident_fij_score calculations but we also need it to validate that transcript
        # sequence has the same reference AA at the mapped SNP position as the protein
        # sequence.
        my $txseq = (PPH::Seq::fastacmd($txname, "$CONFIG{GOLDENPATH}/$CONFIG{GENESET}/genes/knownGeneAA"))[4];
        if (defined $txseq && (my $txlen = length $txseq)) {
          $PRF{Tx}{Seq} = $txseq;
          $PRF{Tx}{Len} = $txlen;
          # Complete set of MZ data loaded successfully
          $PROT->{Profile}{Mz}{$txname} = \%PRF;
        } else {
          warn "ERROR: Failed to fetch $txname protein sequence from the available version of $CONFIG{GENESET} database\n";
        }
      } else {
        warn "ERROR: Unable to locate PSIC profile for $txname in $CONFIG{GENESET} MultiZ alignments\n";
      }

    }

    # Now extract the substitution's MZ PSIC scores for the SNP position
    if (exists $SNP->{Gene}{CdnPos} && defined $SNP->{Gene}{CdnPos}) {
      my $PRF = $PROT->{Profile}{Mz}{$txname};
      if (keys %$PRF && @{ $PRF->{Scores} }) {
        # Validate reference AA at the gene SNP position of the transcript sequence
        my $txaa   = substr($PRF->{Tx}{Seq}, $SNP->{Gene}{CdnPos} - 1, 1);
        my $refkey = $CONFIG{'REVERSEDIRECTION'} ? 'Aa2' : 'Aa1';
        my $varkey = $CONFIG{'REVERSEDIRECTION'} ? 'Aa1' : 'Aa2';
        my $valid  = 0;
        if      ($txaa eq $SNP->{$refkey}) {
          $valid = 1;
        } elsif ($txaa eq $SNP->{$varkey}) {
          substr($PRF->{Tx}{Seq}, $SNP->{Gene}{CdnPos} - 1, 1, $SNP->{$refkey});
          warn "WARNING: Replaced reference AA residue ($txaa) with ($SNP->{$refkey}) in $txname protein sequence at position: $SNP->{Gene}{CdnPos}\n";
          $valid = 1;
        }
        if ($valid) {
          $SNP->{Scores}{Mz}{TxName} = $txname;
          $SNP->{Scores}{Mz}{Psic1}  = $PRF->{Scores}[$SNP->{Gene}{CdnPos}-1][$AAidx{$SNP->{Aa1}}];
          $SNP->{Scores}{Mz}{Psic2}  = $PRF->{Scores}[$SNP->{Gene}{CdnPos}-1][$AAidx{$SNP->{Aa2}}];
          $SNP->{Scores}{Mz}{Nobs}   = $PRF->{Scores}[$SNP->{Gene}{CdnPos}-1][-1];
          $SNP->{Scores}{Mz}{PsicD}  = $SNP->{Scores}{Mz}{Nobs} ?
              sprintf('%+.3f', $SNP->{Scores}{Mz}{Psic1} - $SNP->{Scores}{Mz}{Psic2}) : undef;
        } else {
          warn "ERROR: Neither AA1 ($SNP->{Aa1}) nor AA2 ($SNP->{Aa2}) at position ($SNP->{Pos}) in $SNP->{Acc} query protein sequence matches $txname AA residue $txaa at position ($SNP->{Gene}{CdnPos})\n";
        }
      }
    }
  }

  # Choose which set of the PSIC scores we should use
  my ($nobs_mz, $nobs_msa);
  if ($CONFIG{'MULTIZPATH'}) {
    $nobs_msa = $SNP->{Scores}{Msa}{Nobs} if
      keys %{ $PROT->{Profile}{Msa} }     &&
      exists $SNP->{Scores}{Msa};
    $nobs_mz = $SNP->{Scores}{Mz}{Nobs}   if
      keys %{ $PROT->{Profile}{Mz} }      &&
      exists $SNP->{Scores}{Mz};
  } else {
    $nobs_msa = $SNP->{Scores}{Msa}{Nobs} if
      keys %{ $PROT->{Profile}{Msa} }     &&
      exists $SNP->{Scores}{Msa};
  }
  if      (defined $nobs_msa && defined $nobs_mz) {
    $SNP->{Scores}{Selected} = $nobs_mz >= $nobs_msa ? 'Mz' : 'Msa';
  } elsif (defined $nobs_msa) {
    $SNP->{Scores}{Selected} = 'Msa';
  } elsif (defined $nobs_mz) {
    $SNP->{Scores}{Selected} = 'Mz';
  } else {
    $SNP->{Scores}{Selected} = '';
  }

  #-----------------------------------------------------
  # Calculate alignment / tree based substitution scores
  #-----------------------------------------------------
  PPH::Align::ident_fij_score($PROT, $SNP) if $CONFIG{'SUBSTSCORES'} && $SNP->{Scores}{Selected};

  ##### |     set_PSIC_profile returns (SNP): $SNP
  return 1;
}
#----------------------------------------------------------------------


#------------------------------------------------------------------------------

=head2 calculate_profile

 Usage   : $PPH::Profile::calculate_profile($PROT, $PRF)
 Function: calculates the PSIC matrix
 Args    : References to object of type PROT and PRF
 Returns : 'calculated' if successful

=cut

#----------------------------------------
sub calculate_profile {
  my ($PROT, $PRF) = @_;
  #### PPH__Profile__calculate_profile: $PROT, $PRF

  my $blastfile = $PROT->{Acc} . '.blast';
  my $alnfile   = $PROT->{Acc} . '.aln';

  # Do a BLAST search
  # Skip if either one of the BLAST or alignment files is already present
  my @source = qw{ Missing Precomputed Cached };
  my ($aln_found, $blast_found);
  if      ($aln_found = PPH::Misc::locate_file($alnfile, "$CONFIG{PRECOMPATH}/alignments", "$CONFIG{SCRATCH}/alignments")) {
    $PRF->{Files}{Aln} = $alnfile;
    push @{ $PRF->{$source[$aln_found]} }, 'alignment';
  } elsif ($blast_found = PPH::Misc::locate_file($blastfile, "$CONFIG{PRECOMPATH}/blastfiles", "$CONFIG{SCRATCH}/blastfiles")) {
    $PRF->{Files}{Blast} = $blastfile;
    push @{ $PRF->{$source[$blast_found]} }, 'blast';
  } else {
    # Write protein sequence to a file, if it does not exist already
    my $seqfile = $PROT->{Acc} . '.seq';
    unless (-e $seqfile) {
      open my $W, '>', $seqfile or confess "Can't create file: $seqfile";
      print $W ">" . $PROT->{Acc} . "\n" . $PROT->{Seq} . "\n";
      close $W;
    }
    PPH::BLAST::run_BLAST($seqfile, $blastfile, $CONFIG{'BLAST_ALN'}) or return 'BLAST error';
    copy($blastfile, "$CONFIG{SAVEPATH}/blastfiles/")
      or confess "Failed to save file: $blastfile";
    $PRF->{Files}{Blast} = "$CONFIG{SAVEPATH}/blastfiles/$blastfile";
    push @{ $PRF->{Calculated} }, 'blast';
  }

  # MSA pipeline
  # Fall back to BLAST pairwise pseudo-MSA if MSA_LEON_CLUSPACK method fails
  SWITCH:
  {
    # Skip MSA creation if alignment file is already present
    last if $aln_found;

    # Test whether the BLAST search results are too large to handle
    return 'huge BLAST' if -s $blastfile > $CONFIG{'BIGBLAST'};

    my $status = '';

    # Do not even attempt PPH-2 MSA pipleine when USEBLASTMSA config option is set
    unless ($CONFIG{'USEBLASTMSA'}) {
      # IAA: Use different cutoffs for BLAST_MSA_LEON_CLUSPACK PPHv2 pipeline
      #      vs. old PPHv1 BLAST version
      $CONFIG{'MINPSICHITIDE'} = $CONFIG{'MSA_MINPSICHITIDE'};
      $CONFIG{'MAXPSICHITIDE'} = $CONFIG{'MSA_MAXPSICHITIDE'};
      $PRF->{Version} = 2;  # PPHv2 MSA pipeline version is being used

      $status = PPH::Align::Create_aln_BLAST_MSA_LEON_CLUSPACK($PROT, $PRF, $blastfile, $alnfile);

      # Break out here if PPHv2 MSA pipeline was successful
      chomp $status;
      last SWITCH if $status eq 'OK';

      # Fail immediately if BLAST XML parsing reported a fatal error
      if ($status =~ /\bBLAST XML\b/) {
        warn "ERROR: BLAST search failed ($status), giving up\n";
        return $status;
      # Otherwise warn and proceed with PPHv1 BLAST pipeline if PPHv2 failed
      } else {
        chomp $status;
        warn "WARNING: MSA_LEON_CLUSPACK failed ($status), switching to BLAST alignment\n";
      }
    }

    # IAA: Use different cutoffs for PPHv1 BLAST pipeline
    $CONFIG{'MINPSICHITIDE'} = $CONFIG{'BLAST_MINPSICHITIDE'};
    $CONFIG{'MAXPSICHITIDE'} = $CONFIG{'BLAST_MAXPSICHITIDE'};
    $CONFIG{'BLAST_EVALUE'}  = '1e-4';
    $PRF->{Version} = 1;  # PPHv1 MSA pipeline version is being used

    $status = PPH::Align::Create_aln_from_BLAST_m7($PROT, $PRF, $blastfile, $alnfile);

    chomp $status;
    last SWITCH if $status eq 'OK';

    warn "ERROR: BLAST alignment failed ($status), giving up\n";
    return $status;
  }
  #----------------------------------------

  #----------------------------------------
  # PSIC
  #----------------------------------------
  my $prffile = $PROT->{Acc} . '.prf';  # forget about any old .prf file hanging around
  my $exec    = "$CONFIG{PSIC} $CONFIG{PSIC_ARG} $alnfile $CONFIG{MATRIX}/Blosum62.txt $prffile 2>&1";
  ##### |     calculate_profile executes: $exec
  my $ret = `$exec`;
  if (length $ret) {
    chomp $ret;
    my $msg = "PSIC failed: $ret";
    warn "$msg\n";
    return $msg;
  }
  copy($prffile, "$CONFIG{SAVEPATH}/profiles/")
    or confess "ERROR: Failed to save PSIC profile file: $prffile";
  $PRF->{Files}{Prf} = "$CONFIG{SAVEPATH}/profiles/$prffile";
  push @{ $PRF->{Calculated} }, 'profile';
  #----------------------------------------

  return 'calculated';
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 get_aln_scores_for_SNP

 Usage   : $PPH::Profile::get_aln_scores_for_SNP($psicfile, $pos, $var1, $var2)
 Function: Retrieve the PSIC score for a given SNP
 Args    : PSIC filename
           Amino acid position
           Amino acid variation 1 & 2

=cut

#----------------------------------------
sub get_aln_scores_for_SNP {
  my ($psic_file, $pos, $var1, $var2) = @_;

  #### PPH__Profile__get_aln_scores_for_SNP: $psic_file, $pos, $var1, $var2
  open my $F, '<', $psic_file or confess "Can't read file: $psic_file";
  my @lines = <$F>;
  close $F;

  ($#lines + 1 < $pos)
      and warn "Insufficient No. of lines in profile $psic_file for position: $pos"
          and return(undef, undef, 0);

  my @alphabet = (split /\s+/, $lines[1])[1..20];

  ((grep $var1, @alphabet) and (grep $var2, @alphabet))
      or do {
          warn "Unable to locate variant aa1 ($var1) or aa2 ($var2) in profile: $psic_file";
          return(undef, undef, 0);
      };

  my ($num, @scores) = split /\s+/, $lines[$pos + 1];
  my $Nobs           = pop @scores; # Last column in matrix is Nobs

  ($num != $pos)
      and warn "Can't locate corresponding line in profile $psic_file for pos: $pos"
          and return(undef, undef, 0);

  ($#scores != 19)
      and warn "Invalid line in profile $psic_file for pos: $pos"
          and return(undef, undef, 0);

  my %aaidx;
  @aaidx{@alphabet}  = (0..$#alphabet);

  ##### |     get_aln_scores_for_SNP returns: $scores[$aaidx{$var1}], $scores[$aaidx{$var2}], $Nobs
  return($scores[$aaidx{$var1}], $scores[$aaidx{$var2}], $Nobs);
}
#------------------------------------------------------------------------------

# Reads PSIC profile matrix from a prf-formatted file on disk
# Returns array of arrays of PSIC scores
sub _read_profile {
  my $psicfile = shift;
  ##### |     _load_profile: $psicfile
  open my $F, '<', $psicfile or confess "ERROR: Cannot open profile file: $psicfile";
  # Skip header and column names lines
  <$F>; <$F>;
  my @profile;
  my $lineno = 1;
  while (<$F>) {
    my @a = split;
    # First column is Pos
    my $pos = shift @a;
    confess "ERROR: Wrong position number ($pos) at row $lineno of PSIC profile: $psicfile" if $pos != $lineno;
    # Last column is Nobs
    push @profile, \@a;
    $lineno++;
  }
  close $F;
  ##### |     _load_profile returns: @profile
  return @profile ? \@profile : [];
}

#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

1;
