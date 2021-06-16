package PPH::Align;
use strict;
use warnings;
use Carp qw(cluck confess);
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Align

=head1 DESCRIPTION

Module containing the alignment pipeline

=head1 SUBVERSION

 $LastChangedDate: 2012-02-18 13:28:16 -0500 (Sat, 18 Feb 2012) $
 $LastChangedRevision: 393 $
 $LastChangedBy: ivan $



=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

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
    die "Please install Smart::Comments module or switch off verbose mode / DEBUG configuration option\n" if $@;
  }
  #----------------------------------------
  if ($CONFIG{'MULTIZTREE'}) {
    eval { require Bio::TreeIO; };
    die "Please install Bio::TreeIO module or switch off MULTIZTREE configuration option\n" if $@;
  }
}
our @EXPORT_OK;

use PPH::Data;
use PPH::BLAST;

#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII


#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

=head1 ROUTINES

=head2 Create_aln_BLAST_MSA_LEON_CLUSPACK

  Usage   : Create_aln_BLAST_MSA_LEON_CLUSPACK($PROT, $PRF, $blast_file, $aln_file)
  Function: Reads the BLAST output
            Filters HSP according to CONFIG parameters
            Retrieves full-length sequences
            Aligns with MAFFT (uses timeout - returns)
            Improve alignment quality with LEON (uses timeout - ignores step)
            Clusters sequences with CLUSPACK (uses timeout - ignores step)
  Args    : Object of type PROT (Acc, Desc, Seq)
            Object of type PRF
            BLAST output filename (XML)
            Alignment filename
  Returns : 'OK' if successful

=cut

#----------------------------------------
sub Create_aln_BLAST_MSA_LEON_CLUSPACK {
  my ($PROT, $PRF, $blast_file, $aln_file) = @_;
  #### PPH__Align__Create_aln_BLAST_MSA_LEON_CLUSPACK: $PROT, $PRF, $blast_file, $aln_file

  my $method = 'BLAST'; # Reporting in the end which methods were used

  #----------------------------------------
  # Read BLAST XML file
  #----------------------------------------
  my (@HSPs, $blastmsg);
  $blastmsg = PPH::BLAST::parse_XML_BLAST($blast_file, \@HSPs) and return $blastmsg;
  #----------------------------------------

  my $hsps_raw = scalar @HSPs;

  #----------------------------------------
  # Filter alignments
  #----------------------------------------
  if ($PROT->{Len} < $CONFIG{'MINPSICHITLEN'}) {
    ##### |  If query sequence is too short set MINPSICHITLEN to 75% of the query sequence length
    $CONFIG{'MINPSICHITLEN'} = 0.75 * $PROT->{Len};
  }
  (@HSPs) = grep _filter_PSIC_HSP(), @HSPs;
  my $hsps_filtered = scalar @HSPs;
  return "all $hsps_raw sequences filtered out" unless $hsps_filtered;
  #----------------------------------------

  my $facc = $PROT->{Acc}; $facc =~ s/\./_/; # Some programs are sensitive to extra dots in filenames
  #----------------------------------------
  # Temporary files
  #----------------------------------------
  my $msa_in             = $facc . '_msa_in';
  my $msa_out            = $facc . '_msa_ou';
  my $leon_out           = $facc . '_leon___ou';
  my $cluspack_out       = $leon_out . '.clu';
  my $cluspack_out_tfa   = $leon_out . '2.tfa';
  my $cluspack_out_fixed = $leon_out . '_fixed.tfa';
  #----------------------------------------

  #----------------------------------------
  # Retrieve protein sequences
  #----------------------------------------
  open my $FAS, '>', $msa_in
      or confess "PPH::Align::Create_aln_BLAST_MSA_LEON_CLUSPACK: ".
          "Can't write to $msa_in";
  print $FAS ">QUERY\n$PROT->{Seq}\n";  # "QUERY" helps to find the query later

  my $aa_count  = $PROT->{Len};         # total No of amino acids in alignment
  my $nseq      = 0;                    # No of sequences in alignment (QUERY excluded)
  my $minseqlen = $aa_count;            # length of the shortest sequence
  my $maxseqlen = $aa_count;            # length of the longest sequence
  my %hitnums;                          # hash of matched proteins (hit numbers)
  my %species;                          # hash of species names
  # IAA: Sort HSPs by Identity (now calculated by _calc_ident()), descending
  foreach my $hsp (sort { $b->{'Identity'} <=> $a->{'Identity'} } @HSPs) {

    # Parse species name from hit description
    my $tax;
    if (  # UniRef version
          $hsp->{'Hit_def'} =~ /\sTax=(.+?)(\s+\w+?=|$)/ ||
          # Swiss-Prot version
          $hsp->{'Hit_def'} =~ /\sOS=(.+?)(\s+\w+?=|$)/
        ) {
      $tax = $1;
      $tax =~ s/^\s+//;
      $tax =~ s/\s+$//;
    }

    next if ($CONFIG{'NOREFORG'} || $CONFIG{'NOPARALOGS'}) &&
      exists $CONFIG{'REFORGANISM'} && length $CONFIG{'REFORGANISM'} &&
      defined $tax && length $tax && $tax =~ /$CONFIG{REFORGANISM}/o;

    # The filtering below is heavily based on the fact that HSPs are
    # sorted by descending identity level, i.e. the higest scoring HSPs
    # are the first ones encountered.
    if ($CONFIG{'MULTIHSP'} && !$CONFIG{'USEFULLSEQ'}) {
      # The NOPARALOGS filter below is somewhat questionable since it only
      # takes into consideration single highest scoring HSP for each protein.
      # It is possible that in some cases another hit with several lower
      # scoring HSPs might be a better (more representative) selection instead.
      # However, implementing such scoring scheme would be rather hard.
      if (defined $tax && length $tax) {
        next if exists $species{$tax} && $species{$tax} ne $hsp->{'Hit_num'};
        $species{$tax} = $hsp->{'Hit_num'} if $CONFIG{'NOPARALOGS'};
      }
    } else {
      # Skip multiple HSPs leaving only single highest scoring one for each protein hit.
      next if exists $hitnums{$hsp->{'Hit_num'}};
      $hitnums{$hsp->{'Hit_num'}} = 1;
      # In addition, skip HSPs for paralogs based on species name.
      if (defined $tax && length $tax) {
        next if exists $species{$tax};
        $species{$tax} = 1 if $CONFIG{'NOPARALOGS'} && $tax;
      }
    }

    my ($head, $id, $seq);
    # Retrieve full protein sequence
    if ($CONFIG{'USEFULLSEQ'}) {
      ($head, $id, undef, undef, $seq) =
        PPH::Seq::fastacmd($hsp->{'Hit_id'}, $CONFIG{'NRDB_FASTA'});
      unless (defined($head)) {
        warn "Can't retrieve $hsp->{Hit_id} from $CONFIG{NRDB_FASTA}";
        next;
      }
      # IAA: Replace Selenocysteine (U) with Cysteine (C) and Pyrrolysine (O)
      # with Lysine (K) in the sequence since MAFFT as well as many other tools
      # do not recognize these non-standard residue codes.
      $seq =~ tr/UOuo/CKck/;
      # Replace UniRef DB ID with SwissProt (enusres that full ID
      # string is <=30 characters in length, important for Cluspack)
      $head =~ s/^>gnl\|uniref\d+\|/>sp|/;
    # Use BLAST HSP sequence
    } else {
      $id    = $hsp->{'Hit_accession'};
      $head  = $hsp->{'Hit_id'};
      # Remove entry name in case of Swiss-Prot database
      $head =~ s/^(sp\|\S+)\|\S+/$1/;
      # Replace UniRef DB ID with SwissProt (ensures that full ID
      # string is <=30 characters in length, important for Cluspack)
      $head =~ s/^gnl\|uniref\d+\|/sp|/;
      # Append HSP number to protein ID to make it unique (otherwise
      # many parts of the pipeline will have problems)
      $head  = '>' . $head . '#' . $hsp->{'Hsp_num'};
      $head .= ' ' . $hsp->{'Hit_def'} if length $hsp->{'Hit_def'};
      $seq   = $hsp->{'Hsp_hseq'};
    }
    #----------------------------------------

    #----------------------------------------
    # Limit total residue count and total
    # number of sequences in resultant MSA
    #----------------------------------------
    my $seqlen = length($seq);
    # Hard limit max total number of aa residues in all input sequences
    warn("WARNING: BLAST search results for $PROT->{Acc} truncated due to total residue count > $CONFIG{MSA_SEQLEN_LIMIT}\n"),
      last if $aa_count + $seqlen > $CONFIG{'MSA_SEQLEN_LIMIT'};
    $minseqlen = $seqlen if $minseqlen > $seqlen;
    $maxseqlen = $seqlen if $maxseqlen < $seqlen;
    $aa_count += $seqlen;
    # Discard excessive HSP hits
    warn("WARNING: BLAST search results for $PROT->{Acc} truncated due to total number of HSPs > $CONFIG{MSA_NSEQS_LIMIT}\n"),
      last if $nseq >= $CONFIG{'MSA_NSEQS_LIMIT'};
    $nseq++;
    #----------------------------------------

    print $FAS "$head\n$seq\n";
  }
  close $FAS;
  #----------------------------------------

  $method .= "($hsps_raw,$hsps_filtered,$nseq)";

  #----------------------------------------
  # MSA alignment (timeout / limits)
  #----------------------------------------
  my $msa_arg;
  my $msa_strategy;
  # Check for enforced MSA strategy override
  if      (exists $CONFIG{'MSA1_OVERRIDE'} && length $CONFIG{'MSA1_OVERRIDE'}) {
    $msa_arg = $CONFIG{'MSA_' . $CONFIG{'MSA1_OVERRIDE'}};
    $msa_strategy = $CONFIG{'MSA1_OVERRIDE'};
  # Most accurate but slow (L-INS-i)
  } elsif ($maxseqlen <= $CONFIG{'MSA_SLOW_MAXSEQLEN'} && $nseq <= $CONFIG{'MSA_SLOW_MAXNSEQS'}) {
    $msa_arg = $CONFIG{'MSA_SLOW'};
    $msa_strategy = 'SLOW';
  # Fast but still accurate (FFT-NS-i)
  } elsif ($maxseqlen <= $CONFIG{'MSA_NORM_MAXSEQLEN'} && $nseq <= $CONFIG{'MSA_NORM_MAXNSEQS'}) {
    $msa_arg = $CONFIG{'MSA_NORM'};
    $msa_strategy = 'NORM';
  # Very fast but rough alignment (FFT-NS-2)
  } else {
    $msa_arg = $CONFIG{'MSA_FAST'};
    $msa_strategy = 'FAST';
  }

  eval {
    local $SIG{'ALRM'} = sub { die 'TIMEOUT'; } if $CONFIG{'MSA_T_LIMIT'};
    alarm $CONFIG{'MSA_T_LIMIT'} if $CONFIG{'MSA_T_LIMIT'};
    my $cmd = "$CONFIG{MSA} $msa_arg $msa_in 1>$msa_out 2>>$CONFIG{LIMBO}";
    ##### |  Starting MSA: $cmd
    system($cmd) == 0 or die "Exited with status $?\n";
    alarm 0 if $CONFIG{'MSA_T_LIMIT'};
  };

  if ($@ =~ /^TIMEOUT/) {
    return "MSA_$msa_strategy($facc,$nseq,$maxseqlen) TIMEOUT";
  } elsif ($@) {
    return "MSA_$msa_strategy($facc,$nseq,$maxseqlen) ERROR: $@";
  }

  return "MSA_$msa_strategy($facc,$nseq,$maxseqlen) ERROR: empty output file $msa_out"
    if -z $msa_out;

  $method .= "_MSA_$msa_strategy($nseq,$maxseqlen)";

  #----------------------------------------

  #----------------------------------------
  # LEON - alignment improvement
  #----------------------------------------
  my $status = PPH::LEON::run_LEON($msa_out, 'QUERY', $leon_out);
  if ($status eq 'OK') {
    $method .= '_LEON';
  } else {
    chomp $status;
    warn "WARNING: LEON failed ($status), fall back to MSA_$msa_strategy($facc,$nseq,$maxseqlen)\n";
    # Replace original LEON basenames with MSA-based ones
    $leon_out           = $msa_out;
    $cluspack_out       = $leon_out . '.clu';
    $cluspack_out_tfa   = $leon_out . '2.tfa';
    $cluspack_out_fixed = $leon_out . '_fixed.tfa';
  }
  #----------------------------------------

  #----------------------------------------
  # CLUSPACK - sequence clustering
  #----------------------------------------
  ##### CLUSPACK: $leon_out
  # $leon_out: P04798_leon___ou

  eval {
    local $SIG{'ALRM'} = sub { die 'TIMEOUT'; } if $CONFIG{'CLUSPACK_T_LIMIT'};
    alarm $CONFIG{'CLUSPACK_T_LIMIT'} if $CONFIG{'CLUSPACK_T_LIMIT'};
    my $cmd = sprintf("$CONFIG{CLUSPACK} $CONFIG{CLUSPACK_ARG} >>$CONFIG{LIMBO} 2>&1", $leon_out);
    ##### |  Starting cluspack: $cmd
    if (system($cmd)) {
      if      ($? == -1) {
        die "Failed to execute: $!\n";
      } elsif ($? & 127) {
        die sprintf("Died with signal: %d\n", ($? & 127));
      } else {
        die sprintf("Exited with value: %d\n", $? >> 8);
      }
    }
    alarm 0 if $CONFIG{'CLUSPACK_T_LIMIT'};
  };
  #----------------------------------------

  #-----------------------------------------
  # If CLUSPACK fails then return with error
  # unless DIEONCLUSPACKERROR is NOT set and
  # the error type is a nonzero return value
  # from cluspack executable.
  #-----------------------------------------
  if ($@ =~ /^TIMEOUT/) {
    return 'CLUSPACK TIMEOUT';
  } elsif ($@) {
    if ($CONFIG{'DIEONCLUSPACKERROR'}) {
      return "CLUSPACK ERROR: $@";
    } else {
      if ($@ =~ /^Exited with/) {
        warn "CLUSPACK WARNING: $@";
      } else {
        return "CLUSPACK ERROR: $@";
      }
    }
  }
  #----------------------------------------

  # IAA: Recreate properly formatted, identity level-sorted FASTA MAF from the original
  # MSA-produced MAF and the list of protein identifiers in the cluspack .clu file.
  # We cannot use cluspack-produced MAF since cluspack mangles and reshuffles sequences.
  # Note: this file will include QUERY alignment as the top entry but no_of_alignments
  #       returned does not count QUERY.
  my @cluster = _read_cluspack_clu($cluspack_out);
  # Occasionally, clustering of very small MAFs may result in a top cluster with orphaned
  # QUERY as the only cluster member. Fall back to the original LEON MAF in such case.
  # Note: we still need to extract corresponding deflines from the original MSA MAF since
  #       LEON MAF has them stripped.
  my $no_of_alignments;
  if (scalar @cluster) {
    $no_of_alignments = subset_maf($msa_out, \@cluster, $cluspack_out_fixed);
    $method .= "_CLUSPACK";
  } else {
    warn "WARNING: No alignments left after Cluspack, fall back to LEON\n";
    my (undef, $leon_accs) = _read_fasta_maf($leon_out);
    $no_of_alignments = subset_maf($msa_out, $leon_accs, $cluspack_out_fixed);
  }
  my $final_msa = $facc . '.msa.fa';
  unless ($no_of_alignments) {
    warn "No alignments left after LEON, fall back to MSA_$msa_strategy($facc,$nseq,$maxseqlen)\n";
    my (undef, $msa_accs) = _read_fasta_maf($msa_out);
    $no_of_alignments = subset_maf($msa_out, $msa_accs, $final_msa);
  } else {
    #---------------------------------------------------------------------------
    # IAA: Recreate proper multiple alignment from LEON/Cluspack-filtered subset
    #      of the original multiple alignment
    #---------------------------------------------------------------------------
    # Check the MAF size and select strategy accordingly.
    # Probably superfluous since filtered alignments expected
    # to be relatively small, but better be safe than sorry.
    my ($nseq, $aacount, $maxseqlen, $minseqlen) = _countseq_maf($cluspack_out_fixed);
    $nseq--;  # exclude QUERY from sequence count
    my $msa_arg;
    my $msa_strategy;
    # Check for enforced MSA strategy override
    if      (exists $CONFIG{'MSA2_OVERRIDE'} && length $CONFIG{'MSA2_OVERRIDE'}) {
      $msa_arg = $CONFIG{'MSA_' . $CONFIG{'MSA2_OVERRIDE'}};
      $msa_strategy = $CONFIG{'MSA2_OVERRIDE'};
    # Most accurate but slow (L-INS-i)
    } elsif ($maxseqlen <= $CONFIG{'MSA_SLOW_MAXSEQLEN'} && $nseq <= $CONFIG{'MSA_SLOW_MAXNSEQS'}) {
      $msa_arg = $CONFIG{'MSA_SLOW'};
      $msa_strategy = 'SLOW';
    # Fast but still accurate (FFT-NS-i)
    } elsif ($maxseqlen <= $CONFIG{'MSA_NORM_MAXSEQLEN'} && $nseq <= $CONFIG{'MSA_NORM_MAXNSEQS'}) {
      $msa_arg = $CONFIG{'MSA_NORM'};
      $msa_strategy = 'NORM';
    # Very fast but rough alignment (FFT-NS-2)
    } else {
      $msa_arg = $CONFIG{'MSA_FAST'};
      $msa_strategy = 'FAST';
    }
    my $cluspack_msa = $cluspack_out_fixed . '.msa';
    eval {
      local $SIG{'ALRM'} = sub { die 'TIMEOUT'; } if $CONFIG{'MSA_T_LIMIT'};
      alarm $CONFIG{'MSA_T_LIMIT'} if $CONFIG{'MSA_T_LIMIT'};
        my $cmd = "$CONFIG{MSA} $msa_arg $cluspack_out_fixed 1>$cluspack_msa 2>>$CONFIG{LIMBO}";
      ##### |  Starting MSA: $cmd
      system($cmd) == 0 or die "MSA_$msa_strategy($facc,$nseq,$maxseqlen) ERROR status $?";
      alarm 0 if $CONFIG{'MSA_T_LIMIT'};
    };
    if ($@ =~ /^TIMEOUT/) {
      return "MSA_$msa_strategy($facc,$nseq,$maxseqlen) TIMEOUT";
    } elsif ($@) {
      return "MSA_$msa_strategy($facc,$nseq,$maxseqlen) ERROR: $@";
    }
    return "MSA_$msa_strategy($facc,$nseq,$maxseqlen) ERROR: empty output file $cluspack_msa"
      if -z $cluspack_msa;
    my (undef, $final_accs) = _read_fasta_maf($cluspack_msa);
    $no_of_alignments = subset_maf($cluspack_msa, $final_accs, $final_msa);
    $method .= "_MSA_$msa_strategy($no_of_alignments,$maxseqlen)";
  }

  # Load final FASTA MAF into memory
  my ($clu_maf, $clu_acc) = _read_fasta_maf($final_msa);

  # Validate
  my $qseq = $clu_maf->{'QUERY'}{'seq'};
  confess "Empty/missing QUERY sequence in MAF: $final_msa" unless length $qseq;
  confess "QUERY sequence too short in MAF: $final_msa" if length $qseq < $PROT->{Len};

  # Transform FASTA MAF into (pseudo) CLUSTAL format for PSIC calculation
  open my $CLOUT, '>', $aln_file
    or confess "PPH::Align::Create_aln_BLAST_MSA_LEON_CLUSPACK: Can't create file: $aln_file";
  # Print CLUSTAL file header
  print $CLOUT "CLUSTAL $PROT->{Acc} $method\n\n";
  foreach my $accession (@{ $clu_acc }) {
    # QUERY is not included in CLUSTAL-formatted file
    next if $accession eq 'QUERY';
    # Ungap aligned sequence
    my $alignment = $clu_maf->{$accession};
    my (undef, $seq) = _ungap_HSP_alignment($PROT->{'Seq'}, $qseq, $alignment->{'seq'});
    # Output CLUSTAL-formatted description and sequence
    print $CLOUT _formatDesc4PSIC_aln($alignment->{'tag'}, $alignment->{'desc'}), $seq, "\n";
  }
  close($CLOUT);

  copy($final_msa, "$CONFIG{SAVEPATH}/alignments/$PROT->{Acc}.msa.fa")
    or cluck "Can't save FASTA MSA file: $PROT->{Acc}.msa.fa";

  copy($aln_file, "$CONFIG{SAVEPATH}/alignments/")
    or cluck "Can't save ALN MSA file: $aln_file";

  $PRF->{Files}{Aln} = "$CONFIG{SAVEPATH}/alignments/$aln_file";
  push @{ $PRF->{Calculated} }, 'alignment';

  return 'OK';
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 Create_aln_from_BLAST_m7

 Usage   : PPH::Align->Create_aln_from_BLAST_m7($PROT, $PRF, $blast_file, $aln_file)
 Function: Reads BLAST output and creates a alignment file which is readable
           for PSIC
           It filters for too short / too similar / too weak homologs
           NOTE This implementation allows multiple HSPs per Hit (sequence)
 Args    : Object of type SNP
           BLAST filename
           Alignment filename
 Returns : "OK" if successful

=cut

#----------------------------------------
sub Create_aln_from_BLAST_m7 {
  my ($PROT, $PRF, $blast_file, $aln_file) = @_;
  #### PPH Align Create_aln_from_BLAST_m7: $PROT, $PRF, $blast_file, $aln_file

  #----------------------------------------
  # Read the XML file
  #----------------------------------------
  my (@HSPs, $blastmsg);
  $blastmsg = PPH::BLAST::parse_XML_BLAST($blast_file, \@HSPs) and return $blastmsg;
  #----------------------------------------

  my $hsps_raw = scalar @HSPs;

  #----------------------------------------
  # Filter alignments
  #----------------------------------------
  if ($PROT->{Len} < $CONFIG{'MINPSICHITLEN'}) {
    ##### |   Setting MinPsicLen since sequence too short to 75% of sequence length
    $CONFIG{MINPSICHITLEN} = 0.75 * $PROT->{Len};
  }
  (@HSPs) = grep _filter_PSIC_HSP(), sort { $b->{'Identity'} <=> $a->{'Identity'} } @HSPs;

  my $no_of_alignments = scalar @HSPs;
  return "all $hsps_raw sequences filtered out" unless $no_of_alignments;
  #----------------------------------------

  #----------------------------------------
  # Prepare alignment sequences
  # Print to file
  #----------------------------------------

  my $qseq = $PROT->{Seq};
  my $qlen = $PROT->{Len};
  ##### |      Open: $aln_file
  open my $CLOUT, '>', $aln_file
    or confess "PPH::Align::Create_aln_from_BLAST_m7: Can't create file: $aln_file";
  # IAA: Note, no_of_alignments is putative, see validation test below
  print $CLOUT "CLUSTAL $PROT->{Acc} BLAST($hsps_raw,$no_of_alignments)\n\n";
  my $fa_qseq = '-' x $qlen;
  my @fa_seqs;
  # IAA: Sort HSPs by Identity (now calculated by _calc_ident()), descending
  foreach my $hsp (sort { $b->{'Identity'} <=> $a->{'Identity'} } @HSPs) {

    #----------------------------------------
    # Ungap the query / truncate hit
    #----------------------------------------
    my ($q_hsp, $h_hsp) = _ungap_HSP_alignment(
      substr(
        $qseq,
        $hsp->{'Hsp_query-from'} - 1,
        $hsp->{'Hsp_query-to'} - $hsp->{'Hsp_query-from'} + 1
      ),
      $hsp->{'Hsp_qseq'},
      $hsp->{'Hsp_hseq'});
    #----------------------------------------

    #----------------------------------------
    # Format the sequence for PSIC alignment
    #----------------------------------------
    my $Desc = _formatDesc4PSIC_aln($hsp->{Hit_id}, $hsp->{Hit_def});
    my $Aseq = _position_seq_in_aln($qlen,
                                    $hsp->{'Hsp_query-from'},
                                    $h_hsp);
    #----------------------------------------

    #----------------------------------------
    # Final test - late for debuggin reason
    #----------------------------------------
    unless (length($Aseq) == $qlen) {
      warn join "\n",
          "$hsp->{Hit_id}: length of alignment ",
              length($Aseq)." != $qlen (query_length)",
                  ($hsp->{'Hsp_query-from'} - 1)."\t".
                      ($hsp->{'Hsp_query-to'} -
                            $hsp->{'Hsp_query-from'} + 1),
                      $qseq, $q_hsp, $h_hsp, $Aseq;
      next;
    }
    #----------------------------------------

    #----------------------------------------
    # Format sequences for MSA alignment
    #----------------------------------------
    my $Mdef = '>' . $hsp->{Hit_id} . '#' . $hsp->{Hsp_num};
    $Mdef .= ' ' . $hsp->{Hit_def} if length $hsp->{Hit_def};
    push @fa_seqs, $Mdef . "\n" . $Aseq . "\n";

    ##### |      $Desc $Aseq
    print $CLOUT $Desc, $Aseq, "\n";
  }
  # foreach @HSPs
  #------------------------------
  close $CLOUT;

  # Note: this creates MSA with ungapped query sequence, essentially simply
  # a FASTA-formatted version of the .aln alignment since assembling a true
  # gapped query from BLAST HSPs is a daunting task not worth the time spent.
  my $fa_file = $PROT->{Acc} . '.msa.fa';
  open my $FAOUT, '>', $fa_file
    or confess "PPH::Align::Create_aln_from_BLAST_m7: Can't create file: $fa_file";
  print $FAOUT ">QUERY\n$qseq\n";
  print $FAOUT join('', @fa_seqs);
  close $FAOUT;

  copy($fa_file, "$CONFIG{SAVEPATH}/alignments/$PROT->{Acc}.msa.fa")
    or cluck "Can't save FASTA MSA file: $PROT->{Acc}.msa.fa";
  copy($aln_file, "$CONFIG{SAVEPATH}/alignments/")
    or cluck "Can't save ALN MSA file: $aln_file";
  $PRF->{Files}{Aln} = "$CONFIG{SAVEPATH}/alignments/$aln_file";
  push @{ $PRF->{Calculated} }, 'alignment';

  return 'OK';
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head3 _filter_PSIC_HSP

  Usage   : (@HSP) = grep _filter_PSIC_HSP($pos), @HSP
  Function: Filters a HSP using definitions at CONFIG:
            MINPSICHITLEN, MINPSICHITIDE, MAXPSICHITIDE
  Args    : pos - OPTIONAL - HSP must contain position in query
            HSP (BLAST)
  Returns : TRUE if OK

=cut

#----------------------------------------
sub _filter_PSIC_HSP {
  my $pos = shift;

  # IAA: Fixed incorrect length calculation (was missing + 1)
  my $length = $_->{'Hsp_query-to'} - $_->{'Hsp_query-from'} + 1;

  ($length          < $CONFIG{'MINPSICHITLEN'}) and return 0;
  ($_->{'Identity'} > $CONFIG{'MAXPSICHITIDE'}) and return 0;
  ($_->{'Identity'} < $CONFIG{'MINPSICHITIDE'}) and return 0;

  if (defined($pos)) {
    ###### test: (pos2hit($_, $pos))
    pos2hit($_, $pos) or return 0;
  }

  #### |     _filter_PSIC_HSP returns: "$_->{Hit_accession}#$_->{Hsp_num} Length=$length Identity=$_->{Identity}"
  return 1;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head3 _ungap_HSP_alignment

 Usage   : _ungap_HSP_alignmnet($qseq, $q_hsp, $h_hsp)
 Function: Removes gaps from query sequence in BLAST HSP alignment
 Args    : Unmasked Query sequence (starting at HSP query start)
           Query sequence (masked)
           Sbjct sequence
 Returns : Ungapped Query sequence
           Trimmed  Sbjct sequence

=cut

#----------------------------------------
sub _ungap_HSP_alignment {
    my ($q_orig, $q_seq, $h_seq) = @_;

    chomp($q_orig); chomp($q_seq); chomp($h_seq);
    ##### PPH Align _ungap_HSP_alignment:
    #####  input: $q_orig, $q_seq, $h_seq
    #----------------------------------------
    # The easy part: simply remove gaps
    # from the query sequence and remove
    # corresponding pos in subject sequence
    #----------------------------------------
    my $old_pos = 0;
    while ((my $pos = index($q_seq, '-', $old_pos)) > -1) {
      $old_pos = $pos;
      #----------------------------------------
      # Upps, what happens if there is a X in
      # the query sequence... => just get to
      # the next position
      # IAA: this makes zero sense to me and also this code
      #      does not seem to be doing what was intended anyway,
      #      so let's disable it for the moment
      #----------------------------------------
#       if (substr($q_orig, $pos, 1) eq 'X') {
#         $old_pos++;
#       }
      substr($q_seq, $pos, 1, ''); substr($h_seq, $pos, 1, '');
    }
    #----------------------------------------

    #----------------------------------------
    # IAA: the "tricky" code below is:
    #      a) redundant since initial BLAST is run
    #         with -F F (no filtering) option
    #      b) not clear why we would want masked (X)
    #         residues in the query to be unmasked again
    #      c) even if we do want them unmasked why not to
    #         take the whole original unaligned query
    #         sequence instead?
    # In light of the questions above I am disabling the
    # "unmask query" code below completely.
    #
    # More tricky: remove excessive 'X's by
    # comparing original query seq with given
    # query seq in hsp
    # (BLAST sometimes adds too many 'X' to
    #  the query sequence)
    #                       chunk_orig_start
    #                             |
    #  q_orig VQKLHDFLAHSSEESEETSSPPRLAMNSK
    #  q_seq  VQKLHDFLAHXXXXXXXXXPPRLAMXXX
    #  h_seq  VQKLHDFLAHSEESEETCSSPRLVMXXX
    #                   |        |
    #                 X_start chunk_start
    #
    # - make chunks of query (q_chunks)
    # - locate query chunk in original seq
    # - locate X_start, bit_start

    #----------------------------------------

#     my @q_chunks = split /X+/, $q_seq; # make chunks
#     my $index_start = (@q_chunks) ? length($q_chunks[0]) : 0;
#
#     shift @q_chunks;		# ignore first chunk
#         foreach my $q_chunk (@q_chunks) {
# 	    my $X_start     = index($q_seq, "X",          $index_start);
#             my $chunk_start = index($q_seq, "X".$q_chunk, $index_start) + 1;
#
#             #----------------------------------------
#             # OK, assume that no additional X were
#             # added
#             # if this fails search for the seq chunk
#             # from the pos of the chunk before
#             # Might be a bit problematic with
#             # repetitive sequences as only the first
#             # sequence after the last sequence chunk
#             # will be found
#             #----------------------------------------
#             my $chunk_orig_start = index($q_orig, $q_chunk, $chunk_start);
#             unless ($chunk_orig_start == 0) {
#                 $chunk_orig_start  = index($q_orig, $q_chunk, $index_start);
#             }
#             #----------------------------------------
#
#             my $X_diff           = $chunk_start - $chunk_orig_start;
#
#             ###### |      $q_orig
#             ###### |      $q_seq
#             ###### |      $h_seq
#             ###### |      $q_chunk
#             ###### |      $index_start $X_start $chunk_start $chunk_orig_start $X_diff
#
#             #----------------------------------------
#             # Too many X in q_hsp - shorten them
#             #----------------------------------------
#             if ($X_diff > 0) {
#                 substr($q_seq, $X_start, $X_diff, "");
#                 substr($h_seq, $X_start, $X_diff, "");
#             }
#             #----------------------------------------
#             # Too few so expand!
#             #----------------------------------------
#             elsif ($X_diff < 0) {
#                 substr($q_seq, $X_start, 1, "X" x (1+abs($X_diff)));
#                 substr($h_seq, $X_start, 1, "X" x (1+abs($X_diff)));
#             }
#             #----------------------------------------
#             $index_start = $chunk_orig_start + length($q_chunk);
#         }
    #----------------------------------------

    ###### return: $q_seq, $h_seq
    die 'PPH Align _ungap_HSP_alignment: lengths? ' .
      length($q_orig) . ' ' . length($q_seq) . ' ' . length($h_seq)."\n" .
	    "'$q_orig'\n'$q_seq'\n'$h_seq'"
	      if (length($q_orig) != length($q_seq)) or
           (length($q_orig) != length($h_seq));

    return ($q_seq, $h_seq);
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head3 _position_seq_in_aln

  Usage   : ($aseq) = _position_seq_in_aln($alength, $start, $aseq)
  Function: Simply positions a sequence in an alignment filling missing
            regions with gaps '-'
  Args    : Alignment length
            Start pos of aseq in alignment ($hsp->{'Hsp_query-from'})
            Alignment sequence             ($hsp->{'Hsp_hseq'})
  Returns : Positioned sequence

=cut

#----------------------------------------
sub _position_seq_in_aln {
    my ($alength, $start, $seq) = @_;
    my $aseq = '-' x $alength;
    substr($aseq, $start - 1, length($seq), $seq);
    return $aseq;
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head3 _formatDesc4PSIC_aln

 Usage   : ($desc) = format4PSIC_aln($acc, $def)
 Function: nice formatting of description line 4 PSIC
           - if possible extract taxonomic info and format it
           (This format is readable for PSIC)
 Args    : acc    : accession  ($hsp->{Hit_id})
           def    : definition ($hsp->{Hit_def})
 Returns : Description (correct length)

=cut

#----------------------------------------
sub _formatDesc4PSIC_aln {
  my ($acc, $def) = @_;

  # Try to extract tax info
  my $tax;
  if (defined $def && length $def) {
    chomp($def);
    # NCBI BLAST format
    if      (($tax) = $def =~ / \[(\w+ \w+)\]/) {
      $def =~ s/ \[\w+ \w+\]//;
    # UniRef format
    } elsif (($tax) = $def =~ /\sTax=(.+?)(\s+\w+?=|$)/) {
      $def =~ s/(\sTax=.+?)(\s+\w+?=|$)/$2/;
    # Swiss-Prot format
    } elsif (($tax) = $def =~ /\sOS=(.+?)(\s+\w+?=|$)/) {
      $def =~ s/(\sOS=.+?)(\s+\w+?=|$)/$2/;
    }
    if (defined $tax && length $tax) {
      # Abbreviate tax if too long
      if (length($tax) > 16) {
        my ($g, $s) = split /\s+/, $tax, 2;
        $tax = substr($g, 0, 1) . '.';
        $tax .= " $s" if defined $s && length $s;
      }
      $def = "$tax $def";
    }
  } else {
    $def = '';
  }

  # Remove HSP numbers from tags (if added in multi-HSP mode)
  $acc =~ s/^(\S+)#\d+/$1/;

  # Total description width = 70 columns (67 chars + 3 trailing spaces)
  #return (sprintf '%-21.21s %-45.45s   ', $acc, $def);
  return (sprintf '%-67.67s   ', "$acc $def");
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

=head2 pos2hit

 Usage   : PPH::Align::pos2hit($hsp, $pos)
 Function: Finds in the HSP hit the corresponding AA at query position
 Args    : Hashref to XML-HSP
           Positon of AA in the complete sequence
 Returns : Position in the hit
           Amino acid of query and hit

=cut

#----------------------------------------
sub pos2hit {
    my ($hsp, $pos) = @_;

    #### PPH__Align__pos2hit...

    ####### |       pos, Hsp-query-from, Hsp-query-to: $pos, $hsp->{'Hsp_query-from'}, $hsp->{'Hsp_query-to'}
    return () unless $pos >= $hsp->{'Hsp_query-from'} &&
                     $pos <= $hsp->{'Hsp_query-to'};

    my $Qlen = length $hsp->{Hsp_qseq};
    my $Hlen = length $hsp->{Hsp_hseq};
    confess "ERROR: pos2hit: Query ($Qlen) and hit ($Hlen) sequence lengths do not match"
      unless $Qlen == $Hlen;

    my $gapcha = $CONFIG{BLASTGAP};
    # Find AA in HSP, complicated because of gaps
    my $Qpos = $hsp->{'Hsp_query-from'}; # pos in query
    my $Hpos = $hsp->{'Hsp_hit-from'};   # pos in hit
    my $Apos;                            # pos in alignment (zero-based scale)
    for ($Apos=0; $Apos<$Qlen; $Apos++) {
      substr($hsp->{Hsp_qseq}, $Apos, 1) eq $gapcha or $Qpos++;
      substr($hsp->{Hsp_hseq}, $Apos, 1) eq $gapcha or $Hpos++;
      last if $Qpos == $pos;
    }
    return () if $Apos == $Qlen; # pos not withing the alignment range

    $Apos++;  # adjust to match Qpos / Hpos

    my $Qaa = substr($hsp->{Hsp_qseq}, $Apos, 1);
    my $Haa = substr($hsp->{Hsp_hseq}, $Apos, 1);
    return () if $Qaa eq $gapcha || $Haa eq $gapcha;

    ##### PPH__Align__pos2hit returns: $Hpos, $Qaa, $Haa
    return $Hpos, $Qaa, $Haa;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 filter_3D_HSPs

 Usage   : @HSPs = grep PPH::Align::filter_HSPs($pos), @HSPs
 Function: Filters given HSPs (Hashes with XML keys) according to rule in
           CONFIG
           - This procedure concentrates on the given position
             so it might be better to use another filtering
             procedure in order to ensure recycable BLASTs!
 Args    : Position of aa in Query
           Vias grep: array of HSPs (Hashes with XML keys)
 Returns : Adds following keys to HSP and returns filtered HSPs
           Hit_pos  - Corresponding position of query aa
           Query_aa - The query aa at pos
           Hit_aa   - The hit   aa at pos

=cut

#----------------------------------------
sub filter_3D_HSP {
  my $pos = shift;
  ###### filter_3D_HSP: $pos

  #----------------------------------------
  # General filter - close similarity
  #----------------------------------------
  ####### |       hsp: $_->{'Hit_id'}, $_->{'Hsp_align-len'}, $_->{'Hsp_gaps'}, $_->{'Identity'}
  ($_->{'Hsp_align-len'} < $CONFIG{MIN3DHITLEN})  and return 0;
  ($_->{'Hsp_gaps'}      > $CONFIG{MAX3DHITGAPS}) and return 0;
  ($_->{'Identity'}      < $CONFIG{MIN3DHITIDE})  and return 0;

  # Skip positional mapping / filtering if pos argument is missing
  return 1 unless defined $pos;

  #----------------------------------------
  # Map the position to the hit sequence
  #----------------------------------------
  ####### |       pos: $pos
  my ($Hpos_of_aa, $Qaa, $Haa) = pos2hit($_, $pos);
  defined($Haa) or return 0;

  #----------------------------------------
  # Can AA in hit differ from AA in query?
  #----------------------------------------
  !$CONFIG{MAP3D2MISMATCH}
    and $Qaa ne $Haa
    and $Qaa ne $CONFIG{UNKNOWNRESIDUE}
    and $Haa ne $CONFIG{UNKNOWNRESIDUE}
    and return 0;

  $_->{Hit_pos}  = $Hpos_of_aa;
  $_->{Query_aa} = $Qaa;
  $_->{Hit_aa}   = $Haa;
  return 1;
  ###### |      filter_3D_HSP returns: "OK"
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 ident_fij_score

 Usage   : my ($status) =
              PPH::Align->ident_fij_score($SNP, $PSIC)
 Function: This routines finds all occurance of amino acids at the SNP position
           It multiplies
           a) the sequence identity of the closest sequence containing
              amino acid f_aa (aa2 / aa1, the other is called q_aa) with
           b) the frequency of the q_aa -> f_aa substitution derived
              from the BLOSUM62 frequency matrix
 Args    : hash ref containing SNP info (that includes PSIC alignment info)
 Returns : $SNP->{PSIC}{Aln} =
           {
             Ident -> {f_aa} = ident                  # max identity of f_aa
             IdentQmin       = ident                  # min identity of q_aa
             IdentPmax       = max(fqj * identaaj)    # max P found for any aaj
             IdentPSNP       = max(fq_aa_f_aa * ident_f_aa)
                                                      # max freq * ident of f_aa
           }

=cut

#------------------------------------------------------------------------------

# IAA: Changed to reflect the fact that there is no query in MAF now.
#      Note that this function depends on the MAF being sorted by
#      descending identity level vs. query sequence.
#      Also note, that BLOSUM62N now holds normalized frequencies.
#      01/29/2009: Heavily modified to correct numerous errors in calculations
#                  in the original version.
{
my @ALNKEY;
my @ALNSEQ;
my @AIDENT;
my @ATAX;
my $TREE;
sub ident_fij_score {
  my ($PROT, $SNP) = @_;
  ### PPH__Align__ident_fij_score: $PROT, $SNP

  my $selected = $SNP->{Scores}{Selected};
  return 0 unless $selected;

  my $refaa = $SNP->{Aa1};
  my $mutaa = $SNP->{Aa2};
  my ($query, $qlen, $pos, $alnfile, $alnkey, $alnidx);
  if      ($selected eq 'Msa') {
    $query   = $PROT->{Seq};
    $qlen    = $PROT->{Len};
    $pos     = $SNP->{Pos};
    $alnfile = $PROT->{Profile}{Msa}{Files}{Aln};
    $alnkey  = $SNP->{Acc};
    $alnidx  = 0;
  } elsif ($selected eq 'Mz') {
    my $tx   = $SNP->{Gene}{TxName};
    $query   = $PROT->{Profile}{Mz}{$tx}{Tx}{Seq};
    $qlen    = $PROT->{Profile}{Mz}{$tx}{Tx}{Len};
    $pos     = $SNP->{Gene}{CdnPos};
    $alnfile = $PROT->{Profile}{Mz}{$tx}{Files}{Aln};
    $alnkey  = $tx;
    $alnidx  = 1;
  } else {
   warn "WARNING: Unsupported alignment version selected: $selected\n";
   return 0;
  }

  warn("WARNING: ident_fij_score: Substitution position ($pos) outside query range: [1..$qlen]\n"),
    return 0 if $pos < 1 || $pos > $qlen;

  # Query sequence is only used for calculating ident levels
  substr($query, $pos - 1, 1, $mutaa) if $CONFIG{'REVERSEDIRECTION'};

  my %ALN;

  # Do we want to calculate tree-based scores?
  my $treeflag = $CONFIG{'MULTIZTREE'} && $selected eq 'Mz' ? 1 : 0;

  # Calculate identity scores and cache them together with the current alignment
  unless (defined $ALNKEY[$alnidx] && $ALNKEY[$alnidx] eq $alnkey && $ALNSEQ[$alnidx]) {
    # Clear buffer from prevous alignment data
    delete $ALNSEQ[$alnidx]; delete $AIDENT[$alnidx];
    @ATAX = () if $treeflag;
    open my $F, '<', $alnfile
        or do { warn "WARNING: ident_fij_score: Can't open alignment file: $alnfile\n"; return 0 };
    # Skip first two CLUSTAL header lines but parse first line to determine the
    # version of MSA pipeline used to create the file (unless was already set).
    $_ = <$F>;
    if ($selected eq 'Msa' && !exists $PROT->{Profile}{Msa}{Version}) {
      if (/ BLAST\(\d+,\s*\d+\)\s*$/) {
        $PROT->{Profile}{Msa}{Version} = 1;
      } else {
        $PROT->{Profile}{Msa}{Version} = 2;
      }
    }
    $_ = <$F>;
    while (<$F>) {
      chomp;
      # columns 1-70 are description
      push @{ $ALNSEQ[$alnidx] }, substr($_, 70);
      push @{ $AIDENT[$alnidx] }, _calc_ident($query, $ALNSEQ[$alnidx][-1]);
      if ($treeflag) {
        my $desc = (split(' ', $_, 2))[0];
        push @ATAX, (split(/_/, $desc, 2))[1];
      }
    }
    close($F);
    $ALNKEY[$alnidx] = $alnkey;
  }

  # Load tree used to calculate tree-based scores (currently only defined for multiz alignments)
  if ($treeflag && !defined $TREE) {
    my $treeio = new Bio::TreeIO(-file=>$CONFIG{'MULTIZTREE'}, -format=>'newick');
    $TREE = $treeio->next_tree;	# original tree instance
  }
  my ($tree, $top_node);
  if ($treeflag) {
    $tree = $TREE->clone();								# working tree copy
    $top_node = $TREE->find_node('hg19');	# human leaf
  }

  my @f;
  my ($nseqs, $nobs, $nvars, $nsubs) = (0, 0, 0, 0);
  my %nsites;
  my @tdist;
  my $ident;
  my %ataxa;
  for (my $i=0; $i<@{$ALNSEQ[$alnidx]}; $i++) {
    my $alnlen = length $ALNSEQ[$alnidx][$i];
    # AA at the SNP position in the alignment
    warn("WARNING: ident_fij_score: Substitution position ($pos) outside alignment sequence ($i) range: [1..$alnlen]\n"), next
      if $pos > $alnlen || $pos < 1;
    my $s_aa  = substr($ALNSEQ[$alnidx][$i], $pos - 1, 1);
    warn("WARNING: ident_fij_score: Position ($pos) undefined for alignment sequence: $i\n"), next
      unless $s_aa;
    $nseqs++;
    # Note that gaps in alignment are completely igonred!
    next if $s_aa eq $CONFIG{CLUSTALGAP};
    next unless exists $PPH::Data::BLOSUM62N{$s_aa};
    my $tax;
    if ($treeflag) {
      $tax = $ATAX[$i];
      $ataxa{$tax} = 1;
    }
    $nobs++;
    # Upon the loop exit, $ident will hold identity level for the last (i.e., most divereged) sequence in alignment
    $ident = $AIDENT[$alnidx][$i];
    # Skip all sequences in alignment that are conserved at the SNP pos (vs. query)
    next if $s_aa eq $refaa;
    # Past this point we have a substitution in the alignment (vs. query)
    $ALN{Sites} .= $s_aa;
    if ($treeflag) {
      my $current_node = $TREE->find_node($tax);
      push @tdist, sprintf("%.3f", $TREE->distance(-nodes => [$top_node, $current_node]));
    }
    $nsites{$s_aa} = 1;
    $nsubs++;
    # Minimum and maximum sequence identity level for each substitution occurring at the SNP pos (not used)
    $ALN{IdentMin}{$s_aa} = $ident if !exists($ALN{IdentMin}{$s_aa}) || $ALN{IdentMin}{$s_aa} > $ident;
    $ALN{IdentMax}{$s_aa} = $ident if !exists($ALN{IdentMax}{$s_aa}) || $ALN{IdentMax}{$s_aa} < $ident;
    #
    # This should really be called IdentQmax since this is a maximum identity level
    # (or closest divergence distance) at which first substitution (vs. query) is
    # observed in the alignment at the SNP position
    $ALN{IdentQmin} = $ident unless exists($ALN{IdentQmin});
    my $fij   = $PPH::Data::BLOSUM62N{$s_aa}{$mutaa};
    my $idfij = $ident * $fij;
    ####### |     ident_fij_score analyzing: $i, $refaa, $s_aa, $ident, $fij
    # Maximum of { P(a1,a2) * sequence identity } for the variant aa (mutaa, aa2) if it occures at the SNP pos
    if ($s_aa eq $mutaa) {
      $ALN{IdentPSNP} = $idfij if !exists($ALN{IdentPSNP}) || $ALN{IdentPSNP} < $idfij;
      $nvars++;
    }
    # Maximum of { P(a1,a2) * sequence identity } for all aa pairs at the SNP pos among all sequences in alignment
    $ALN{IdentPmax} = $idfij if !exists $ALN{IdentPmax} || $ALN{IdentPmax} < $idfij;
  }
  $ALN{Nseqs}  = $nseqs;
  $ALN{Nsubs}  = $nsubs;
  $ALN{Nvars}  = $nvars;
  $ALN{Nsites} = scalar keys %nsites;
  if ($treeflag) {
    foreach my $node ($tree->get_leaf_nodes) {
      my $id = $node->id;
      next if $id eq 'hg19';	# keep human branch that is never present in .aln file
      unless ($ataxa{$id}) {
        $tree->remove_Node($id);
        $tree->contract_linear_paths(1);
      }
    }
    $ALN{Ntree} = $tree->total_branch_length;
    $ALN{Tdist} = join(',', @tdist);
  } else {
    $ALN{Ntree} = '';
    $ALN{Tdist} = '';
  }
  $SNP->{Scores}{Aln} = \%ALN;
  ###### |     ident_fij_score returns: %ALN
  return 1;
}
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Create a subset MAF containing only sequences matching a list of protein
# accessions; sort sequences in the output MAF by their identity level
# vs. QUERY (descending)
#
# Parameter(s) passed:
#   input MAF filename, ref to a list of protein accessions, output MAF filename
#   optional strip_description flag, description part of defline is stripped if set
#   optional skip_query flag, QUERY entry is not included in the output MAF if set
# Value(s) returned:
#   number of entries included in the output MAF, not counting
#   QUERY (regardless of skip_query flag value)
#
sub subset_maf {
  my ($input_ma_file, $subset, $output_ma_file, $strip_description, $skip_query) = @_;

  # Read original raw MAF
  my ($ali, $keys) = _read_fasta_maf($input_ma_file);

  my @acclist = @$keys;
  confess "No QUERY entry in MAF: $input_ma_file" unless exists $ali->{'QUERY'};
  my $query = $ali->{'QUERY'};
  confess "QUERY sequence missing from MAF: $input_ma_file" unless length $query->{seq};
  for (my $i=0; $i<@acclist; $i++) {
    my $acc = $acclist[$i];
    if ($acc eq 'QUERY') {
      $ali->{$acc}{ident} = 1.0;
    } else {
      $ali->{$acc}{ident} = _calc_ident($query->{seq}, $ali->{$acc}{seq});
      # Remove any previous identity level annotation(s) if present in input MAF
      $ali->{$acc}{desc} =~ s/ Identity\s*=\s*[01]?\.(\d+)//g;
    }
  }

  my @cluster = @$ali{ @$subset };

  open(FOUT, ">$output_ma_file") or confess "Can't create output file: $output_ma_file";
  print(FOUT ">QUERY\n", $query->{seq}, "\n") unless $skip_query;
  my $no_of_alignments = 0;
  foreach my $alignment (sort {$b->{ident} <=> $a->{ident}} @cluster) {
    # Skip QUERY if it was included in the subset list of accessions
    next if $alignment->{tag} eq 'QUERY';
    $no_of_alignments++;
    print(FOUT '>', $alignment->{tag});
    unless ($strip_description) {
      print(FOUT ' ', $alignment->{desc}) if length $alignment->{desc};
      print(FOUT ' ', sprintf('Identity=%.3f', $alignment->{ident}))
        if exists $alignment->{ident} && defined $alignment->{ident};
    }
    print(FOUT "\n", $alignment->{seq}, "\n");
  }
  close(FOUT);

  return $no_of_alignments;
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Read and parse a Cluspack .clu file with the list of sequence clusters
#
# Parameter(s) passed:
#   filename
# Value(s) returned:
#   list of accession numbers of the proteins which belong to the cluster
#   containing QUERY sequence (excluding QUERY identifier itself)
#
sub _read_cluspack_clu {
  my $fname = shift;

  my @cluster;
  my $queryflag = 0;
  open(FIN, $fname) or confess "Can't open Cluspack clusters file: $fname";
  <FIN>;  # Skip header line (e.g. "Number of clusters : 5")
  while (<FIN>) {
    chomp;
    next if /^\s*$/;
    # New cluster encountered (or a block of unclustered sequences)
    if (/^Cluster\s+\d+\s*;\s+size=\d+\s*$/ ||
        /^unclustered\s+points\s+\d+\s*;\s*size=\d+\s*$/) {
      # Return a previous cluster containing query if already encountered
      # Note: @cluster list can be empty, check on caller's side!
      return @cluster if $queryflag;
      @cluster = ();
      next;
    }
    # Parse protein FASTA tag
    my $acc = ( _parse_defline($_) )[2];
    # Query sequence is in this cluster
    if ($acc eq 'QUERY') {
      $queryflag = 1;
    } else {
      push @cluster, $acc;
    }
  }
  close(FIN);

  # Only a single cluster was present in the .clu file or
  # QUERY-containing cluster was the last one in the file
  if ($queryflag) {
    return @cluster;  # @cluster list can be empty, check on caller's side!
  } else {
    confess "No query sequence found among Cluspack clusters: $fname";
  }
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Read and parse a multiple alignment file in FASTA format
#
# Parameter(s) passed:
#   filename
# Value(s) returned:
#   ref to hash of alignment data
#   ref to list of protein accession numbers in MAF order
#
sub _read_fasta_maf {
  my $fname = shift;

  #   Amino acid codes accepted in alignments:
  #
  #   A  alanine               O  pyrrolysine (not an IUPAC/IUBMB standrad yet)
  #   B  aspartate/asparagine  P  proline
  #   C  cystine               Q  glutamine
  #   D  aspartate             R  arginine
  #   E  glutamate             S  serine
  #   F  phenylalanine         T  threonine
  #   G  glycine               U  selenocysteine
  #   H  histidine             V  valine
  #   I  isoleucine            W  tryptophan
  #   J  leucine/isoleucine    Y  tyrosine
  #   K  lysine                Z  glutamate/glutamine
  #   L  leucine               X  any
  #   M  methionine            *  translation stop
  #   N  asparagine            -  gap of indeterminate length
  #
  my $AACODE  = 'ABCDEFGHIKLMNOPQRSTUVWYZX*-';

  open(FIN, $fname) or confess "Can't open MAF: $fname";
  my $maf;
  { local $/ = undef;
    $maf = <FIN>;
  }
  close(FIN);
  my @alignments = split /^>/m, $maf;
  shift @alignments;  # remove leading empty field
  confess "Empty MAF or incorrect format: $fname" unless @alignments;
  my $defline;
  my @acclist;
  my %alihash;
  my $alicount = 0;
  foreach my $alignment (@alignments) {
    $alicount++;
    $defline = substr($alignment, 0, index($alignment, "\n") + 1, '');
    my ($tag, $db, $acc, $name, $desc) = _parse_defline($defline);
    $alignment =~ s/\s+//sg;  # sequence
    confess "Non-unique protein accession ($acc) in alignment ($alicount): $fname"
      if exists $alihash{$acc};
    # Validate sequence
    confess "Illegal sequence character in alignment ($alicount): $fname"
      if $alignment =~ /[^$AACODE]/o;
    $alihash{$acc}{tag}  = $tag;
    $alihash{$acc}{db}   = $db;
    $alihash{$acc}{name} = $name;
    $alihash{$acc}{desc} = $desc;
    $alihash{$acc}{seq}  = $alignment;
    push @acclist, $acc;
  }

  return \%alihash, \@acclist;
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Parse FASTA definition line in NCBI BLAST format
# (very simplistic version)
#
# Parameter(s) passed:
#   FASTA definition line (as string)
# Value(s) returned:
#   list of 5 strings:
#     tag  - complete unparsed FASTA tag, aka protein ID, e.g. 'sp|Q3LFU0|CP1A1_BALAC'
#     db   - primary database code, e.g. 'sp'
#     acc  - primary protein accession number, e.g. 'Q3LFU0'
#     name - primary protein db entry name, e.g. 'CP1A1_BALAC'
#     desc - FASTA description
# Note: all return values are optional except acc and set to empty
#       strings if missing from defline
#
sub _parse_defline {
  my $defline = shift;
  confess "Missing FASTA definition line" unless defined $defline && length $defline;
  chomp($defline);
  $defline =~ s/^>//;
  my ($tag, $desc) = split /\s+/, $defline, 2;
  $tag =~ s/^\s+//;
  confess "Malformed or missing identifier in FASTA defline: $defline"
    unless defined $tag && length $tag;
  # Get primary db name, protein accession number and name, e.g.:
  #
  #  Database Name                         Identifier Syntax
  #
  #  GenBank                               gb|accession|locus
  #  EMBL Data Library                     emb|accession|locus
  #  DDBJ, DNA Database of Japan           dbj|accession|locus
  #  NBRF PIR                              pir||entry
  #  Protein Research Foundation           prf||name
  #  SWISS-PROT                            sp|accession|entry name
  #  Brookhaven Protein Data Bank          pdb|entry|chain
  #  Patents                               pat|country|number
  #  GenInfo Backbone Id                   bbs|number
  #  General database identifier           gnl|database|identifier
  #  NCBI Reference Sequence               ref|accession|locus
  #  Local Sequence identifier             lcl|identifier
  my ($db, $acc, $name) = ('', '', '');
  my @f = split /\|/, $tag;
  # Simple unformatted tag
  if      (@f == 1) {
    $acc  = $f[0];
  # NCBI GenBank format, e.g., gi|129295|sp|P01013|OVAX_CHICK (db entry name is optional!)
  } elsif (@f > 3 && $f[0] eq 'gi') {
    $db   = $f[2];
    $acc  = $f[3];
    $name = $f[4] if @f > 4;
  # Some databases do not use accession numbers, e.g.: pir||entry
  } elsif (@f == 3 && $f[1] eq '' && length $f[2]) {
    $db   = $f[0];
    $acc  = $f[2];
  # Generic database, e.g.: gnl|uniref100|B4DTD3
  } elsif ($f[0] eq 'gnl') {
    $db   = $f[1];
    $acc  = $f[2];
  # Local database, e.g.: lcl|B4DTD3
  # or a non-standard variant used by PPH internally: sp|Q9R207
  } elsif (@f == 2 && ($f[0] eq 'lcl' || $f[0] eq 'sp')) {
    $db   = $f[0];
    $acc  = $f[1];
  # Most common case, e.g.: sp|accession|entry_name
  } else {
    $db   = $f[0];
    $acc  = $f[1];
    $name = $f[2];
  }
  confess "Empty accession number in FASTA defline: $defline"
    unless length $acc;
  # Description is optional
  if (defined $desc && length $desc) {
    $desc =~ s/\s+$//;
  } else {
    $desc = '';
  }
  return $tag, $db, $acc, $name, $desc;
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Calculate identity level for BLAST query vs, BLAST match protein sequences
#
# Parameter(s) passed:
#   seq1 - _aligned_ query sequence
#   seq2 - BLAST hit sequence
# Value(s) returned:
#   fraction of identical residues
#
sub _calc_ident {
  my ($seq1, $seq2) = @_;
  my $len1 = length $seq1;
  my $len2 = length $seq2;
  # IAA 11/06/2010: Bail out with a fatal error if lengths of two
  #                 sequences do not match.
  confess "ERROR: _calc_ident() Lengths of sequences do not match (len1=$len1, len2=$len2)\n" if $len1 != $len2;
  # Translate all unknown and ambiguous aa codes into gaps.
  # Actually, this works by translating everything that is
  # _not_ a valid aa code (or is not already a gap).
  $seq1 =~ tr/ACDEFGHIKLMNOPQRSTVWY-/-/c;
  $seq2 =~ tr/ACDEFGHIKLMNOPQRSTVWY-/-/c;
  my ($valid, $ident, $gaps) = (0, 0, 0);
  for my $i (0 .. $len1-1) {
    my $a1 = substr($seq1, $i, 1);
    my $a2 = substr($seq2, $i, 1);
    # IAA 03/25/09: Completely ignore (skip) all gaps in query sequence (seq1).
    #               This has also an effect of skipping all masked and/or ambiguous
    #               residues in query (see tr/// lines above).
    next if $a1 eq '-';
    if ($a2 ne '-') {
      $valid++;
      $ident++ if $a1 eq $a2;
    } else {
      $gaps++;
    }
  }
  # IAA 03/25/09: Weight identity rate by gaps (effectively counting gaps as mismatches).
  $valid += $gaps;
  if ($valid) { return $ident / $valid; }
  else        { return 0; }
}

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

# Count number and length of sequences in MAF
# Parameter(s) passed:
#   filename (MAF in FASTA format)
# Value(s) returned:
#   nseq      - number of sequences in alignment
#   aacount   - total residues in all sequences
#   maxseqlen - maximum sequence length (sans gaps)
#   minseqlen - minimum sequence length (sans gaps)
#
sub _countseq_maf {
  my $mafile = shift;

  my ($nseq, $aacount, $maxseqlen, $minseqlen) = (0, 0, 0, 0);
  my ($alignments, undef) = _read_fasta_maf($mafile);
  foreach my $pacc (keys %{ $alignments }) {
    my $alignment = $alignments->{$pacc};
    # Count everything except spaces and gaps
    my $seqlen = $alignment->{seq} =~ tr/ -//c;
    $minseqlen = $seqlen if $minseqlen > $seqlen;
    $maxseqlen = $seqlen if $maxseqlen < $seqlen;
    $aacount  += $seqlen;
    $nseq++;
  }

  return $nseq, $aacount, $maxseqlen, $minseqlen;
}

#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
1;
