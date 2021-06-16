package PPH;
use Carp qw(cluck confess);
use warnings;
use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH - main PolyPhen-2 module

=head1 SYNOPSIS


 my $PPH = PPH->new('my_sequences.fa');  # Filename argument is optional
 $pph->analyze($SNP);
 $pph->interpret_results($SNP);
 print $pph->return_header(), "\n";
 print $pph->print_result(), "\n";


=head1 DESCRIPTION

Access to PolyPhen routines

=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2009-11-28 00:15:09 -0500 (Sat, 28 Nov 2009) $
 $LastChangedRevision: 254 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION = 0.1;
    @ISA           = qw(Exporter);
    @EXPORT        = qw(&new);
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

use File::Copy;
use XML::Simple;

use PPH::Seq;
use PPH::Pfam;
use PPH::Gene;
use PPH::Profile;
use PPH::3D;
use PPH::Rules;
# use PPH::BLAT;

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

# Header for tabulated output: 50 columns
our @tabheader = qw(
  o_acc
  o_pos
  o_aa1
  o_aa2
  rsid
  acc
  pos
  aa1
  aa2
  nt1
  nt2
  prediction
  based_on
  effect
  site
  region
  PHAT
  dScore
  Score1
  Score2
  MSAv
  Nobs
  Nstruct
  Nfilt
  PDB_id
  PDB_pos
  PDB_ch
  ident
  length
  NormASA
  SecStr
  MapReg
  dVol
  dProp
  B-fact
  H-bonds
  AveNHet
  MinDHet
  AveNInt
  MinDInt
  AveNSit
  MinDSit
  Transv
  CodPos
  CpG
  MinDJxn
  PfamHit
  IdPmax
  IdPSNP
  IdQmin
);
if ($CONFIG{'ALNSCORES'}) { # optional alignment-based scores: 8 columns
  push @tabheader,
		'Len',
    'Nseqs',
    'Nsites',
    'Nsubs',
    'Nvars',
    'Ntree',
    'Sites',
    'Tdist'
  ;
}

# Tabulated output format
our @tabformat = (
  # Original query values
  '%-20s', #  1    o_acc
  '%6s',   #  2    o_pos
  '%5s',   #  3    o_aa1
  '%5s',   #  4    o_aa2
  # Mapped query
  '%-10s', #  5    rsid
  '%-10s', #  6    acc
  '%6s',   #  7    pos
  '%3s',   #  8    aa1
  '%3s',   #  9    aa2
  '%3s',   # 10    nt1
  '%3s',   # 11    nt2
  # PolyPhen-1 annotations
  '%18s',  # 12    prediction
  '%20s',  # 13    based_on
  '%10s',  # 14    effect
  # Features
  '%8s',   # 15    site
  '%8s',   # 16    region
  '%8s',   # 17    PHAT
  '%6s',   # 18    dScore
  '%6s',   # 19    Score1
  '%6s',   # 20    Score2
  '%4s',   # 21    MSAv
  '%6s',   # 22    Nobs
  '%8s',   # 23    Nstruct
  '%6s',   # 24    Nfilt
  '%6s',   # 25    PDB_id
  '%7s',   # 26    PDB_pos
  '%6s',   # 27    PDB_ch
  '%6s',   # 28    ident
  '%6s',   # 29    length
  '%7s',   # 30    NormASA
  '%6s',   # 31    SecStr
  '%6s',   # 32    MapReg
  '%6s',   # 33    dVol
  '%6s',   # 34    dProp
  '%6s',   # 35    B-fact
  '%8s',   # 36    H-bonds
  '%8s',   # 37    AveNHet
  '%8s',   # 38    MinDHet
  '%8s',   # 39    AveNInt
  '%8s',   # 40    MinDInt
  '%8s',   # 41    AveNSit
  '%8s',   # 42    MinDSit
  '%6s',   # 43    Transv
  '%6s',   # 44    CodPos
  '%3s',   # 45    CpG
  '%8s',   # 46    MinDJxn
  '%12s',  # 47    PfamHit
  '%8s',   # 48    IdPmax
  '%8s',   # 49    IdPSNP
  '%8s'    # 50    IdQmin
);
if ($CONFIG{'ALNSCORES'}) { # optional alignment-based scores
  push @tabformat,
    '%6s', # 51    Len
    '%6s', # 52    Nseqs
    '%6s', # 53    Nsites
    '%6s', # 54    Nsubs
    '%6s', # 55    Nvars
    '%6s', # 56    Ntree
    '%s',  # 57    Sites
    '%s',  # 58    Tdist
  ;
}

#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#----------------------------------------------------------------------

=head1 ROUTINES

=head2 new

 Usage:    $pph->new();
           $pph->new('sequences.fa')
 Function: Initialize Object PPH
 Args    : FASTA file containing the protein sequences for the SNPs
           (optional)
 Returns : PPH object

=cut

#----------------------------------------
sub new {
  my ($class, $FASTA_file) =  @_;
  ### [<now>] Starting new PPH instance...
  my $self = {};
  bless($self, $class);
  defined $FASTA_file and -e $FASTA_file and do {
    #----------------------------------------
    # Do a simple test if formatting of
    # uniprot.seq file is like:
    # >id
    # seq
    # >id
    # ...
    #----------------------------------------
    open my $F, '<', "$CONFIG{UNIPROT}.seq"
        or confess "Can't open file $CONFIG{UNIPROT}.seq";
    <$F>; <$F>; $_ = <$F>;
    /^>/ or confess "Wrong format of $CONFIG{UNIPROT}.seq\n".
        "Sequence should be in one line only\n$_\n";
    close $F;
    #----------------------------------------
    $self->{'FASTA_file'} = $FASTA_file;
  };
  return $self;
}
#----------------------------------------------------------------------

#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

#======================================================================

=head2 prepareInput

 Usage   : PPH->prepareInput($SNP)
 Function: This routine checks if protein of SNP is contained in a local
           UniProtKB copy, otherwise retrieves sequence from a user-supplied
           sequence file
 Args    : Object of type SNP
 Returns : Empty string if successful; SNP object will contain all the info

=cut

#========================================
# IAA: Changes to allow sequence extraction from user-supplied sequence database
#      to fail with a warning (instead of fatal error); in such case an attempt is
#      made to extract sequence with the matching protein ID from UniProt/Swiss-Prot
sub prepareInput {
  my ($self, $SNP) = @_;

  # We ONLY have a dbSNP rsID submitted
  my $dbsnp = $SNP->{Query}{rsID};

  if ($dbsnp) {
    ($SNP->{Acc}, $SNP->{Pos}, $SNP->{Aa1}, $SNP->{Aa2}) = PPH::Seq::snp2acc($SNP->{Query}{rsID});
    return "Reference dbSNP $SNP->{Query}{rsID} was not found in the available version of UniProtKB"
      unless length $SNP->{Acc};
  }

  my $qacc = $SNP->{Acc}; # Original query protein accession (before mapping, except for dbSNP rsID)
  $SNP->{Pkey} = $qacc;   # This is a hash key for the protein associated with this SNP

  # Only attempt to map protein accession and extract its sequence if it has not been tried already.
  # If this attempt fails then we will store an empty hash as a stub indicating that further attempts
  # on the same protein will be futile. Use the original query protein accession for Protein hash
  # key since this is the only one we know at this point.
  unless (exists $self->{Protein}{$qacc}) {

    my %PROT;

    # The empty hashref stub indicates that mapping has been attempted. It will be replaced by
    # hashref to valid data only if processing completes successfully.
    $self->{Protein}{$qacc} = {};

    $PROT{UniProt}   = 0; # Assume initially that the protein is not in UniProtKB
    $PROT{Canonical} = 0; # Assume initially that the protein is non-canonical (e.g., an isoform)

    # Attempt to extract sequence from the user-supplied sequence database if
    # specified with -s option (skip if this is a dbSNP SNP).
    if (exists $self->{FASTA_file} && !$dbsnp) {
      ($PROT{Acc}, $PROT{Name}, $PROT{Desc}, $PROT{Seq}) = PPH::Seq::fetch_FASTA($qacc, $self->{FASTA_file});
    }

    # Sequence successfully extracted, now attempt to match it to one of UniProtKB sequences
    if ($CONFIG{'MAP2UNIPROT'} && defined $PROT{Seq} && length $PROT{Seq}) {

      my ($acc, $name, $desc, $seq) = PPH::Seq::grep_in_FASTA($PROT{Seq}, "$CONFIG{UNIPROT}.seq");
      # Exact (sub)sequence match found in UniProtKB
      if (defined $acc) {
        # Offset from the user-submitted sequence into matched UniProtKB sequence.
        $PROT{Offset}    = index($seq, $PROT{Seq});
        # Replace user-submitted sequence with the matched UnoProtKB one
        # ToDo: Should we also save the user-submitted sequence?
        $PROT{Seq}       = $seq;
        $PROT{Name}      = $name;
        $PROT{Acc}       = $acc;
        $PROT{Desc}      = $desc;
        $PROT{UniProt}   = 1;
        $PROT{Canonical} = 1 unless $acc =~ /-\d+$/;
      }

    # Either user sequence file was not specified, or sequence failed to be extracted from it,
    # or this is a dbSNP refSNP. Attempt to match and extract sequence by a protein accession.
    } else {

      # Check if the protein ID submitted is from one of the supported alien databases
      # (UniProtKB, RefSeq, GenBank, Ensembl) and map it to a primary UniProtKB accession.
      # For a dbSNP SNP, copy SNP accession to protein since it's already been mapped.
      if ($dbsnp) {
        $PROT{Acc} = $SNP->{Acc};
      } else {
        $PROT{Acc} = PPH::Seq::id2acc($qacc);
        $PROT{Acc} = PPH::Seq::xdb2acc($qacc) unless length $PROT{Acc};
      }

      return "Unable to locate protein entry $qacc in the available version of UniProtKB"
        unless length $PROT{Acc};

      # Fetch sequence from the built-in UniProtKB database by its accession
      (undef, undef, undef, $PROT{Desc}, $PROT{Seq}) =
        PPH::Seq::fastacmd($PROT{Acc}, $CONFIG{'UNIPROT'});

      return "Failed to extract $qacc query sequence from the available version of UniProtKB"
        unless exists $PROT{Seq} && length $PROT{Seq};

      # At this point the submitted variant has been successfully mapped to a UniProtKB entry.
      $PROT{UniProt}      = 1;
      $PROT{Canonical}    = 1 unless $PROT{Acc} =~ /-\d+$/;
    }

    # IAA: Replace Selenocysteine (U) with Cysteine (C) and Pyrrolysine (O)
    # with Lysine (K) in the sequence since MAFFT as well as many other tools
    # do not recognize these non-standard residue codes.
    $PROT{Seq} =~ tr/UOuo/CKck/;
    $PROT{Len} = length $PROT{Seq};

    $self->{Protein}{$qacc} = \%PROT; # processing successfully completed

  } # end of protein accession mapping / sequence extraction

  my $PROT = $self->{Protein}{$qacc};

  # Do not attempt further processing if the protein sequence has not been validated
  return "Processing of query entry $qacc aborted due to previous error(s)" unless keys %$PROT;
  # Update SNP protein accession to match mapped protein accession
  $SNP->{Acc} = $PROT->{Acc};
  # Adjust the SNP position if necessary
  $SNP->{Pos} += $PROT->{Offset} if $PROT->{Offset};
  # Validate AA substitution submitted against query protein sequence
  return "Invalid variation position ($SNP->{Pos}) specified for $SNP->{Acc} query sequence, valid positions are [1..$PROT->{Len}]"
    if $SNP->{Pos} < 1 || $SNP->{Pos} > $PROT->{Len};
  # Query sequence residue at the position of substitution
  my $qres = substr($PROT->{Seq}, $SNP->{Pos} - 1, 1);
  return "Neither AA1 ($SNP->{Aa1}) nor AA2 ($SNP->{Aa2}) in input matches $SNP->{Acc} query sequence residue ($qres) at position ($SNP->{Pos})"
    if $SNP->{Aa1} ne $qres && $SNP->{Aa2} ne $qres;
  # Default substitution direction is:
  #        (aa1==qres)->aa2
  # but if REVERSEDIRECTION is set then we will look instead at:
  #   aa1->(aa2==qres)
  # This is sometimes useful for analysing mutations in reference
  # organism (e.g., human) lineage
  my $refkey  = $CONFIG{'REVERSEDIRECTION'} ? 'Aa2' : 'Aa1';
  if ($SNP->{$refkey} ne $qres) {
    if ($CONFIG{'SWAPRESIDUES'}) {
      warn "WARNING: Swapped input residues AA1 ($SNP->{Aa1}) and AA2 ($SNP->{Aa2}) for $SNP->{Acc} query sequence at position ($SNP->{Pos})\n";
      # Swap AA residues
      ($SNP->{Aa1}, $SNP->{Aa2}) = ($SNP->{Aa2}, $SNP->{Aa1});
      # Also swap nucleotides when present
      if (exists $SNP->{Nt1} && exists $SNP->{Nt2}) {
        ($SNP->{Nt1}, $SNP->{Nt2}) = ($SNP->{Nt2}, $SNP->{Nt1});
        # Recalculate CpG if it was present in input
        if (exists $SNP->{CpG}) {
          if (exists $SNP->{Flanks}) {
            $SNP->{CpG} = PPH::Gene::CpG(substr($SNP->{Flanks}, 0, 1), $SNP->{Nt1}, $SNP->{Nt2}, substr($SNP->{Flanks}, 1, 1));
          } else {
            delete $SNP->{CpG};
          }
        }
      }
    } else {
      return "Wrong amino acid residue ($qres) found instead of ($SNP->{$refkey}) in $SNP->{Acc} query sequence at position ($SNP->{Pos})";
    }
  }

  # If this is a canonical Swiss-Prot protein and no dbSNP rsID was supplied in input
  # then try to locate dbSNP SNP matching this substitution in Swiss-Prot annotations.
  if ($CONFIG{'MAP2DBSNPS'} && !$dbsnp && $PROT->{Canonical}) {
    my $aapair = join('', sort($SNP->{Aa1}, $SNP->{Aa2}));
    my $rsid = PPH::Seq::acc2snp($SNP->{Acc}, $SNP->{Pos}, $aapair);
    $SNP->{rsID} = $rsid if defined $rsid && length $rsid;
  }

  return '';  # Success
}
#======================================================================


#======================================================================

=head2 analyze

 Usage   : PPH->analyze($SNP)
 Function: This is the main routine of PolyPhen:

             o It validates the SNP object by:
               a) Checking if SNP is not a synonymous substitution
               b) The protein sequence is given either by UniProt
                  (preferrably) or by the sequence file (PolyPhen-2
                  checks if sequence is nevertheless in UniProt)
             o It finds UniProt annotations for the SNP (e.g. transmembrane
               region) or which are in proximity to the SNP residue
             o It calculates the PSIC Profile scores for the
               substitution based on a sequence alignment
             o It checks structural parameters and checks for possible
               interactions with other UniProt-annotated sites

  Args    : Object of type SNP
  Returns : undef if successful

=cut

#========================================
sub analyze {
    my ($self, $SNP) = @_;
    ### PPH__analyze: $SNP

    my $PROT = $self->{Protein}{$SNP->{Pkey}};

    #----------------------------------------
    # Check for UniProt annotations
    ### check: $CONFIG{CALCUNIPROT} and $PROT->{Canonical}
    if ($CONFIG{'CALCUNIPROT'} and $PROT->{'Canonical'}) {
      ### <now> CALCUNIPROT started...
      PPH::Seq::get_uniprot_annotation($PROT, $SNP);
      ### <now> CALCUNIPROT finished...
    }

    #----------------------------------------
    # Check for Pfam annotations
    ### check: $CONFIG{FINDPFAM} and $PROT->{Canonical}
    if ($CONFIG{'FINDPFAM'} and $PROT->{'Canonical'}) {
      ### <now> FINDPFAM started...
      PPH::Pfam::find_pfam($PROT, $SNP);
      ### <now> FINDPFAM finished...
    }

    #----------------------------------------
    # Find gene / transcript encoding the query protein
    ### check: $CONFIG{MAPGENE} && ! exists $SNP->{Gene}
    if ($CONFIG{'MAPGENE'} && ! exists $SNP->{Gene}) {
      ### <now> MAPGENE started...
      PPH::Gene::find_gene($PROT, $SNP);
      ### <now> MAPGENE finished...
    }

    #----------------------------------------
    # Calculate PSIC scores
    ### check: $CONFIG{CALCPRF}
    if ($CONFIG{'CALCPRF'}) {
      ### <now> CALCPRF started...
      PPH::Profile::set_PSIC_profile($PROT, $SNP);
      ### <now> CALCPRF finished...
    }

    #----------------------------------------
    # Calculate 3D
    ### check: $CONFIG{CALC3D}
    if ($CONFIG{'CALC3D'}) {
      ### <now> CALC3D started...
      PPH::3D::set_3D_parameters($PROT, $SNP);
      ### <now> CALC3D finished...
    }

    #### PPH__analyze returns: $SNP
    return 1;
}
#======================================================================


#======================================================================

=head2 interpret_results

 Usage   : PPH->interpret_results($SNP)
 Function: This routine evaluates all parameters given for the SNP
           object and stores them in $SNP->{Rules}
 Args    : Object of type SNP
 Returns : TRUE if successful - SNP will contain the prediction

=cut

#========================================

sub interpret_results {
    my ($self, $SNP) = @_;
    ### PPH__interpret_results: $SNP
    $SNP->{Rules} = PPH::Rules::test($self->{Protein}{$SNP->{Pkey}}, $SNP);
    #### Upon return: $SNP->{Rules}
    return 1;
}
#======================================================================

#======================================================================

=head2 return_header

 Usage   : PPH->return_header([-format=>1])
 Args    : Optional -format flag to return formatted column headers
 Returns : Array containing all column headers as either raw or
           formatted strings depending on the format flag

=cut

#========================================
sub return_header {
  my ($self, %p);
  if (ref($_[0])) { ($self, %p) = @_; } # object method call
  else            { %p = @_;          } # package subroutine call
  if ($p{'-formatted'} || $p{'formatted'}) {
    my @formatted = ();
    for (my $i=0; $i<@tabheader; $i++) {
      push @formatted, sprintf($tabformat[$i], $tabheader[$i]);
    }
    return @formatted;
  } else {
    return @tabheader;
  }
}
#======================================================================


#======================================================================

=head2 format_result

 Usage   : PPH->format_result($SNP)
 Function: Prints values of SNP object which were filled by PolyPhen
 Args    : Object of type SNP
 Returns : Array of values computed from PolyPhen

=cut

#========================================
sub format_result {
    my ($self, $SNP) = @_;
    ##### PPH__format_result

    my @out;
    #----------------------------------------
    # Original query input values.
    push @out, (
      _value($SNP, 'Query', 'rsID') || _value($SNP, 'Query', 'Acc'),
      _value($SNP, 'Query', 'Pos'),
      _value($SNP, 'Query', 'Aa1'),
      _value($SNP, 'Query', 'Aa2')
    );
    #----------------------------------------

    #----------------------------------------
    # Mapped input values plus nucleotides
    push @out, (
      _value($SNP, 'rsID'),
      _value($SNP, 'Acc'),
      _value($SNP, 'Pos'),
      _value($SNP, 'Aa1'),
      _value($SNP, 'Aa2'),
      _value($SNP, 'Gene', 'Nt1'),
      _value($SNP, 'Gene', 'Nt2')
    );
    #----------------------------------------

    #----------------------------------------
    # Predictions
    push @out, (
      _value($SNP, 'Rules', 'Prediction'),
      _value($SNP, 'Rules', 'Basis'),
      _value($SNP, 'Rules', 'Effect')
    );
    #----------------------------------------

    #----------------------------------------
    # UniProt annotations if available
    #
    # Ternary logic for nominal attributes:
    #   annotations seeked but not found    - NO
    #   annotations not checked or missing  - <empty>
    #   annotation found                    - <value>
    #
    if ($CONFIG{'CALCUNIPROT'} and $self->{Protein}{$SNP->{Pkey}}{'Canonical'}) {
      my $features;
      if (keys %{ $SNP->{Features} }) {
        # One (first) feature per each type
        $features->{Site}   = $SNP->{Features}{Sites}[0]   if exists $SNP->{Features}{Sites};
        $features->{Region} = $SNP->{Features}{Regions}[0] if exists $SNP->{Features}{Regions};
        $features->{PHAT}   = $SNP->{Features}{PHAT}       if exists $SNP->{Features}{PHAT};
      }
      my $site   = _value($features, 'Site');
      my $region = _value($features, 'Region');
      push @out, (
        length $site   ? $site   : 'NO',
        length $region ? $region : 'NO',
        _value($features, 'PHAT')
      );
    } else {
      push @out, ('') x 3;
    }
    #----------------------------------------

    #----------------------------------------
    # General alignment info
    my $selected = $SNP->{Scores}{Selected};
    my ($scores, $MSAv);
    if      ($selected eq 'Msa') {
      $scores = $SNP->{Scores}{Msa};
      $MSAv   = $self->{Protein}{$SNP->{Pkey}}{Profile}{Msa}{Version};
    } elsif ($selected eq 'Mz') {
      $scores = $SNP->{Scores}{Mz};
      $MSAv   = 3;
    } else {
      $MSAv   = '';
    }
    push @out, (
      _value($scores, 'PsicD'),
      _value($scores, 'Psic1'),
      _value($scores, 'Psic2'),
      $MSAv,
      _value($scores, 'Nobs'),
      _value($self->{Protein}{$SNP->{Pkey}}, 'Structure', 'HitsAll'),
      _value($SNP, 'Structure', 'HitsMapped') ? $SNP->{Structure}{HitsMapped} : ''
    );
    #----------------------------------------

    #----------------------------------------
    # Structure part
    # - get info for first hit only
    #----------------------------------------
    my ($map, $cont);
    $map = $SNP->{Structure}{Maps}[0] if keys %{ $SNP->{Structure} };
    unless (defined $map) {
      push @out, ('') x 18;
    } else {
      my $nhsp = $CONFIG{CONTALLHITS} ? scalar @{ $SNP->{Structure}{Maps} } : 1;
      $cont = $SNP->{Structure}{BestContact} if exists $SNP->{Structure}{BestContact};
      push @out, (
        _value($map, 'Hsp', 'PDB_id'),
        _value($map, 'Hsp', 'PDB_pos'),
        _value($map, 'Hsp', 'PDB_chain'),
        sprintf('%.2f', _value($map, 'Hsp', 'Identity')),
        _value($map, 'Hsp', 'Hsp_align-len'),
        _value($map, 'Params', 'NormASA') ?
           sprintf('%.3f', $map->{'Params'}{'NormASA'}) : '',
        _value($map, 'Params', 'SecStr'),
        _value($map, 'Params', 'MapReg'),
        _value($map, 'Params', 'dVol'),
        _value($map, 'Params', 'dProp'),
        _value($map, 'Params', 'NormB'),
        _value($map, 'Params', 'Hbonds') ?
           scalar @{ $map->{Params}{Hbonds} } : '',
        _value($cont, 'Number', 'Heteroatoms') ?
           sprintf('%.2g', $cont->{Number}{Heteroatoms} / $nhsp) : '',
        _value($cont, 'Heteroatoms', 'Dist'),
        _value($cont, 'Number', 'InterChain') ?
           sprintf('%.2g', $cont->{Number}{InterChain} / $nhsp) : '',
        _value($cont, 'InterChain', 'Dist'),
        _value($cont, 'Number', 'CritSites') ?
           sprintf('%.2g', $cont->{Number}{CritSites} / $nhsp) : '',
        _value($cont, 'CritSites', 'Dist')
      );
    }
    #----------------------------------------

    #----------------------------------------
    # Adding nucleotide sequence-based info:
    #   Transv
    #   CodPos (frame)
    #   CpG
    #   MinDJxn
    #----------------------------------------
    push @out, _value($SNP, 'Gene', 'Transv');
    push @out, _value($SNP, 'Gene', 'Frame');
    push @out, _value($SNP, 'Gene', 'CpG');
    my $dx1 = _value($SNP, 'Gene', 'JXdon');
    my $dx2 = _value($SNP, 'Gene', 'JXacc');
    $dx1    = 10000 if $dx1 eq '';
    $dx2    = 10000 if $dx2 eq '';
    my $dx  = abs($dx1) < abs($dx2) ? $dx1 : $dx2;
    if ($dx == 10000) {
      push @out, '';
    } else {
      push @out, sprintf('%+d', $dx);
    }

    #----------------------------------------
    # Adding Pfam (first domain found only)
    #----------------------------------------
    if ($CONFIG{'FINDPFAM'} and $self->{Protein}{$SNP->{Pkey}}{'Canonical'} and $self->{Protein}{$SNP->{Pkey}}{Pfam}) {
      push @out, exists $SNP->{Features}{Domains} ? $SNP->{Features}{Domains}[0] : 'NO';
    } else {
      push @out, '';
    }

    #----------------------------------------
    # Adding alignment-based scores
    # IdPmax
    # IdPSNP
    # IdPQmin
    #----------------------------------------
    if (_value($SNP, 'Scores', 'Aln')) {
      my $val = _value($SNP, 'Scores', 'Aln', 'IdentPmax');
      push @out, ($val eq '') ? '' : sprintf('%.3f', $val * 100);
      $val = _value($SNP, 'Scores', 'Aln', 'IdentPSNP');
      push @out, ($val eq '') ? '' : sprintf('%.3f', $val * 100);
      $val = _value($SNP, 'Scores', 'Aln', 'IdentQmin');
      push @out, ($val eq '') ? '' : sprintf('%.2f', $val * 100);
    } else {
      push @out, ('') x 3;
    }

    #----------------------------------------
    # Optional alignment-based scores
    # Len
    # Nseqs
    # Nsites
    # Nsubs
    # Nvars
    # Ntree
    # Sites
    # Tdist
    #----------------------------------------
    if ($CONFIG{'ALNSCORES'}) {
			my $val;
			$val = _value($SNP, 'Len');
			push @out, ($val eq '') ? '' : sprintf('%d', $val);
      if (_value($SNP, 'PSIC', 'Aln')) {
          $val = _value($SNP, 'PSIC', 'Aln', 'Nseqs');
          push @out, ($val eq '') ? '' : sprintf('%d', $val);
          $val = _value($SNP, 'PSIC', 'Aln', 'Nsites');
          push @out, ($val eq '') ? '' : sprintf('%d', $val);
          $val = _value($SNP, 'PSIC', 'Aln', 'Nsubs');
          push @out, ($val eq '') ? '' : sprintf('%d', $val);
          $val = _value($SNP, 'PSIC', 'Aln', 'Nvars');
          push @out, ($val eq '') ? '' : sprintf('%d', $val);
          $val = _value($SNP, 'PSIC', 'Aln', 'Ntree');
          push @out, ($val eq '') ? '' : sprintf('%.3f', $val);
          $val = _value($SNP, 'PSIC', 'Aln', 'Sites');
          push @out, ($val eq '') ? '' : $val;
          $val = _value($SNP, 'PSIC', 'Aln', 'Tdist');
          push @out, ($val eq '') ? '' : $val;
      } else {
					push @out, ('') x 7;
      }
    }

    # Let's ignore the above formatting and start from scratch ;)
    map { if (defined) { s/^\s+//; s/\s+$//; } else { $_ = ''; } } @out;

    my @formatted = ();
    for (my $i=0; $i<@out; $i++) {
      push @formatted, sprintf($tabformat[$i], $out[$i]);
    }

    # Add comments if present
    push @formatted, '# ' . $SNP->{Query}{Comments} if exists $SNP->{Query}{Comments};

    return @formatted;
}

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP

#----------------------------------------

=head3 _value

 Usage   : _value($hashref, @keys)
 Function: This small routine returns the value of a multidim hash
           - it carefully checks if each of the keys really exists
             and otherwise returns undef
 Args    : Hashref
           List of keys
 Result  : Value of innermost key

=cut

#----------------------------------------
sub _value {
    my ($hr, @keys) = @_;
    scalar(@keys) or return defined $hr ? $hr : '';
    my $k = shift @keys;
    exists $hr->{$k} ? _value($hr->{$k}, @keys) : '';
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head2 pph2xml

 Usage   : pph2xml($$);
           $PPH->snp2xml(\*STDOUT);
           PPH::pph2xml($PPH, \*F);
 Function: Returns XML output of the data structure of PPH object
 Args    : PPH object
           Ref to open file handle or glob
 Returns : TRUE if successful;
 IAA     : Changed to use XMLout()

=cut

#----------------------------------------
sub pph2xml {
  my ($self, $FH) = @_;
  #### pph2xml: $self, $FH
  confess 'Object is not of type PPH' unless ref $self eq 'PPH';
  # File::Temp filehandle is incompatible with XMLout
  delete $self->{FASTA_file} if exists $self->{FASTA_file};
  # We do not need full complexity of PPH internals in the XML version
  foreach my $pkey (keys %{ $self->{Protein} }) {
    my $PROT = $self->{Protein}{$pkey};
    if (exists $PROT->{Profile}) {
      delete $PROT->{Profile}{Msa}{Scores} if exists $PROT->{Profile}{Msa} && exists $PROT->{Profile}{Msa}{Scores};
      foreach my $tx (keys %{ $PROT->{Profile}{Mz} }) {
        delete $PROT->{Profile}{Mz}{$tx}{Scores} if exists $PROT->{Profile}{Mz}{$tx}  && exists $PROT->{Profile}{Mz}{$tx}{Scores};
      }
    }
    delete $PROT->{Structure}{Hsps} if exists $PROT->{Structure} && exists $PROT->{Structure}{Hsps};
    delete $PROT->{Structure} if exists $PROT->{Structure} && ! keys %{ $PROT->{Structure} };
    delete $PROT->{Features}{CritSites} if exists $PROT->{Features} && ! keys %{ $PROT->{Features}{CritSites} };
    delete $PROT->{Features} if exists $PROT->{Features} && ! keys %{ $PROT->{Features} };
    if (exists $PROT->{Gene}) {
      foreach my $txname (keys %{ $PROT->{Gene} }) {
        my $TX = $PROT->{Gene}{$txname};
        next unless keys %$TX;
        delete $TX->{id} if exists $TX->{id};
        delete $TX->{TxSeq} if exists $TX->{TxSeq};
        delete $TX->{Cds2Tx} if exists $TX->{Cds2Tx};
        delete $TX->{chrExStarts} if exists $TX->{chrExStarts};
        delete $TX->{chrExEnds} if exists $TX->{chrExEnds};
      }
      delete $PROT->{Gene} if exists $PROT->{Gene} && ! keys %{ $PROT->{Gene} };
    }
  }
  foreach my $SNP (@{ $self->{Variant} }) {
    delete $SNP->{Features} if exists $SNP->{Features} && ! keys %{ $SNP->{Features} };
    if ($SNP->{Structure} && exists $SNP->{Structure}{Maps}) {
      foreach my $map (@{ $SNP->{Structure}{Maps} }) {
        delete $map->{Hsp}{Hsp_qseq} if exists $map->{Hsp} && $map->{Hsp}{Hsp_qseq};
        delete $map->{Hsp}{Hsp_hseq} if exists $map->{Hsp} && $map->{Hsp}{Hsp_hseq};
        delete $map->{Hsp}{Hsp_midline} if exists $map->{Hsp} && $map->{Hsp}{Hsp_midline};
      }
    }
  }
  XMLout($self, OutputFile=>$FH, RootName=>'PPH', XMLDecl=>1, AttrIndent=>1, KeyAttr=>[], SuppressEmpty=>1);
  return 1;
}
#----------------------------------------------------------------------

sub xml2pph {
  my ($self, $FH) = @_;
  #### xml2pph: $self, $FH
  $self = XMLin($FH, KeyAttr=>[], ForceArray=>[ 'Variant', 'Prediction', 'Avail', 'Maps', 'Domains' ]);
  bless($self);
  return $self;
}

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
1;
