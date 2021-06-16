package PPH::3D;
use strict;
use warnings;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::3D - routines for analyzing structures - former pph_3d.pl

=head1 DESCRIPTION

Depends on the PDB database...

=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2011-09-30 16:17:16 -0400 (Fri, 30 Sep 2011) $
 $LastChangedRevision: 370 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use Carp qw(cluck confess);
use Fcntl qw(:flock);
use File::Copy;
use Storable;
use DBI;
use PPH::Data;
use PPH::Misc;
use File::Basename;

BEGIN {
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

#GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

=head1 ROUTINES

=head2 set_3D_parameters

 Usage   : PPH::3D::set_3D_parameters($SNP)
 Function: Find homologous protein structures for the given SNP and retrieve
           structural information
 Args    : Object of type SNP
 Returns : TRUE if successful
           $SNP->{Struct} will contain all information

=cut

#----------------------------------------
sub set_3D_parameters {
  my ($PROT, $SNP) = @_;
  #### PPH__3D__set_3D_parameters

  # Attempt a BLAST search against the database of PDB protein sequences and
  # load resultant BLAST HSPs into memory, if it has not been tried already.
  unless (exists $PROT->{Structure}) {

    my %STR;
    # The empty hashref stub indicating profile creation was attempted.
    # Replaced with ref to valid data if processing completes successfully.
    $PROT->{Structure} = {};

    my $blastfile = $PROT->{Acc} . '.pdb.blast';

    # Locate BLAST file or run BLAST if the file cannot be found.
    # Note that we have to lock the whole pipeline BEFORE checking
    # for the file availability to avoid race conditions.
    $CONFIG{'SAVEPATH'} = $CONFIG{'SAVEKNOWN'} && $PROT->{UniProt} ? $CONFIG{'PRECOMPATH'} : $CONFIG{'SCRATCH'};
    my $lockfile = "$CONFIG{SAVEPATH}/lock/$PROT->{Acc}.lock";
    open my $flk, '>', $lockfile or die "ERROR: Can't create lockfile $lockfile: $!\n";
    flock($flk, LOCK_EX) or die "ERROR: Failed to obtain lock for $lockfile: $!\n";
    ####### |     file locked: $lockfile
    unless (PPH::Misc::locate_file($blastfile,
      "$CONFIG{PRECOMPATH}/structures", "$CONFIG{SCRATCH}/structures")) {
      # Calculation of 3D not allowed
      close($flk), return(0) if $CONFIG{'PRECOMPONLY3D'};
      # Write protein sequence to a file, if it does not exist already
      my $seqfile = "$PROT->{Acc}.seq";
      unless (-e $seqfile) {
        open my $W, '>', $seqfile or confess "Can't create file: $seqfile";
        print $W ">" . $PROT->{Acc} . "\n" . $PROT->{Seq} . "\n";
        close $W;
      }
      PPH::BLAST::run_BLAST($seqfile, $blastfile, $CONFIG{'BLAST_3D'}) or do {
        cluck 'WARNING: BLAST_3D: Multiple accession numbers';
        close($flk);
        return 0;
      };
      copy($blastfile, "$CONFIG{SAVEPATH}/structures")
        or warn "Failed to save file: $blastfile\n";
    }
    close($flk);

    # Parse BLAST XML file and return list of HSPs.
    # The keys correspond to BLAST XML plus Hit-id (Id) and Hit-desc (Desc).
    # Returns non-empty string in case of XML parsing error or no hits found
    my @HSPs;
    PPH::BLAST::parse_XML_BLAST($blastfile, \@HSPs) and return 0;
    $STR{HitsAll} = scalar @HSPs;

    # Filter HSPs according to CONFIG settings
    @HSPs = grep PPH::Align::filter_3D_HSP(), @HSPs;

    # Abort if no valid hits found
    if (@HSPs) {
      $STR{HitsFiltered} = scalar @HSPs;
    } else {
      return 0;
    }

    # Sort HSPs by identity if requested
    @HSPs = sort {$b->{Identity} <=> $a->{Identity}} @HSPs if $CONFIG{'SORTBYIDE'};
    # Load filtered HSPs into memory
    $STR{Hsps} = \@HSPs;

    $PROT->{Structure} = \%STR; # processing successfully completed

  }

  # Do not attempt further processing if the protein structure was not available
  return 0 unless keys %{ $PROT->{Structure} };

  # Map filtered HSPs to the SNP position in query protein sequence and extract
  # various positional 3D parameters from structural databases and files.
  my $hsp_no    = 1;
  my $hsp_count = 0;
  my $dbname = $CONFIG{'PDB_ID_SUFFIX'};
  my %best_contact; # best contact in all HSPs
  foreach my $hsp (@{ $PROT->{Structure}{Hsps} }) {

    # Map SNP position to HSP (FASTA PDB) sequence
    my ($Hpos, $Qaa, $Haa) = PPH::Align::pos2hit($hsp, $SNP->{Pos});
    defined($Haa) or next;

    # Can AA in the HSP differ from AA in the query sequence?
    !$CONFIG{MAP3D2MISMATCH}
      and $Qaa ne $Haa
      and $Qaa ne $CONFIG{UNKNOWNRESIDUE}
      and $Haa ne $CONFIG{UNKNOWNRESIDUE}
      and next;

    $hsp_count++; # total count of HSPs successfully mapped to SNP position

    $hsp->{Hit_pos}  = $Hpos;
    $hsp->{Query_aa} = $Qaa;
    $hsp->{Hit_aa}   = $Haa;

    # Remove gi identifier prepended by some versions of BLAST
    $hsp->{Hit_id} =~ s/^gi\|[^|]+\|//;

    my %Map;
    $Map{Hsp} = $hsp;

    # Check if Id format is pdb|id|chain and modify Id and Desc if necessary
    unless ($hsp->{Hit_id} =~ /^$dbname\|\S+\|\S+/) {
      if ($hsp->{Hit_def} =~ /^$dbname\|\S+\|\S+\s+/) {
        ($hsp->{Hit_id}, $hsp->{Hit_def}) = split /\s+/, $hsp->{Hit_def}, 2;
      } else {
        warn "Problems identifying database $CONFIG{PDB_ID_SUFFIX} in '$hsp->{Hit_id}' or '$hsp->{Hit_def}'";
        next;
      }
    }

    # Set PDB_id and PDB_Chain
    (undef, $hsp->{PDB_id}, $hsp->{PDB_chain}) = split /\|/, $hsp->{Hit_id};

    # Find PDB atom position corresponding to the SNP position in the HSP
    $hsp->{PDB_pos} = pos2pdb($hsp->{PDB_id},
                              $hsp->{PDB_chain},
                              $hsp->{Hit_pos});
    ##### |    PDB_pos: $hsp->{PDB_pos}

    if ($hsp_no == 1 || $CONFIG{STRUCTALLHITS}) {
      #----------------------------------------
      # DSSP info
      #----------------------------------------
      my $DSSP = get_DSSP($PROT, $SNP->{Aa1}, $SNP->{Aa2}, $hsp);
      ##### |    set_3D_parameters DSSP: $DSSP
      defined $DSSP and keys %$DSSP and $Map{Params} = $DSSP;
      #----------------------------------------

      #----------------------------------------
      # H-bonds
      #----------------------------------------
      ##### |     Hit_pos or Hit_aa?
      my @HBonds = $CONFIG{CALCHBONDS} ? get_HBonds($PROT, $hsp) : ();
      scalar @HBonds and $Map{Params}{Hbonds} = \@HBonds;
      #----------------------------------------

      #----------------------------------------
      # B-Factor
      #----------------------------------------
      $Map{Params}{NormB} = get_BFactor($hsp);
      #----------------------------------------
    }
    # end of StructAllHits

    # Check for possible SNP contacts to residues in UniProt-annotated sites in the protein
    if ($hsp_no == 1 || $CONFIG{CONTALLHITS}) {

      ##### |    This needs testing
      my %ValidContacts;
      # Map annotations of critical sites to structure
      foreach my $type (keys %{ $PROT->{Features}{CritSites} }) {
        foreach my $site (@{ $PROT->{Features}{CritSites}{$type} }) {
          # 1st map site pos in query (Swiss-Prot) sequence to hit pos (FASTA PDB)
          my ($pos, undef, undef) = PPH::Align::pos2hit($hsp, $site);
          defined($pos) or next;
          # 2nd map hit (FASTA PDB) pos to pos in PDB (PDB ATOM residue)
          $pos = pos2pdb($hsp->{PDB_id},
                         $hsp->{PDB_chain},
                         $pos);
          defined($pos) or next;
          # concatenate pos & chain
          push @{ $Map{AsMapped} }, ($pos . $hsp->{PDB_chain});
        }
      }

      #----------------------------------------
      # Get contacts of SNP in structure
      # (if there are sites also look for intrachain contacts - last argument)
      #----------------------------------------
      my $FoundContacts = get_contacts($hsp, exists $Map{AsMapped});
      #----------------------------------------

      #----------------------------------------
      # Check for critical Sites
      #----------------------------------------
      if (keys %$FoundContacts and exists $Map{AsMapped}) {
        ###### Map{AsMapped}: $Map{AsMapped}
        foreach my $c (@{ $FoundContacts->{INT} }) {
          grep { $c->{Pos}.$c->{Chain} eq $_ } @{ $Map{AsMapped} } or next;
          push  @{ $ValidContacts{CritSites} }, $c;
          (!defined($best_contact{CritSites}{Dist})
            or $c->{Dist} < $best_contact{CritSites}{Dist})
            and $best_contact{CritSites} = Storable::dclone($c);
          $best_contact{Number}{CritSites}++;
        }
      }
      #----------------------------------------

      #----------------------------------------
      # HeteroAtoms
      #----------------------------------------
      foreach my $c (@{ $FoundContacts->{HET} }) {

        #----------------------------------------
        # Ignore contacts to non biol ligands
        #----------------------------------------
        next if grep { $_ eq $c->{Name} } @PPH::Data::Non_biol_ligands;
        #----------------------------------------

        push @{ $ValidContacts{Heteroatoms} }, $c;
        $best_contact{Number}{Heteroatoms}++;

        #----------------------------------------
        # Need to clone it here to get a clean
        # data structure without refs
        #----------------------------------------
        (!exists($best_contact{Heteroatoms})
          or $c->{Dist} < $best_contact{Heteroatoms}{Dist})
          and $best_contact{Heteroatoms} = Storable::dclone($c);
      }
      #----------------------------------------

      #----------------------------------------
      # InterChain
      #----------------------------------------
      if (exists $FoundContacts->{EXT} and @{ $FoundContacts->{EXT} }) {

        #----------------------------------------
        # Grep the contacts
        #----------------------------------------
        @{ $ValidContacts{InterChain} } = @{ $FoundContacts->{EXT} };
        #----------------------------------------

        #----------------------------------------
        # Check for the closest contact
        #----------------------------------------
        foreach my $c (@{$ValidContacts{InterChain}}) {
          (!exists($best_contact{InterChain})
            or $c->{Dist} < $best_contact{InterChain}{Dist})
            and $best_contact{InterChain} = Storable::dclone($c);
          $best_contact{Number}{InterChain}++;
        }
        #----------------------------------------
      }
      # end interchain
      #----------------------------------------

      $Map{Contacts} = \%ValidContacts if keys %ValidContacts;
    }
    # end of ContAllHits

    # Note: all successfully filtered 3D HSP maps are saved if either one
    # or both of STRUCTALLHITS and/or CONTALLHITS options are set; otherwise
    # only map for a single top hit is saved (this is the default behavior).
    if ($hsp_no == 1 || $CONFIG{STRUCTALLHITS} || $CONFIG{CONTALLHITS}) {
      push @{ $SNP->{Structure}{Maps} }, \%Map;
    }

    $hsp_no++;
  } # end of foreach @HSPs

  $SNP->{Structure}{HitsMapped}  = $hsp_count;
  $SNP->{Structure}{BestContact} = \%best_contact if keys %best_contact;

  ###### |     set_3D_parameters returns: $PROT->{Structure}, $SNP->{Structure}
  return 1;
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head2 get_DSSP

  Usage   : PPH::3D->get_DSSP($aa1, $aa2, $hsp)
  Function: Retrieve secondary structure information for a given substitution
  Args    : Amino acid variation 1 & 2
            Object of type HSP (from BLAST)
  Returns : Hashref containing DSSP information

=cut

#----------------------------------------
sub get_DSSP {
    my ($PROT, $aa1, $aa2, $hsp) = @_;
    #### PPH__3D__get_DSSP: $aa1, $aa2, $hsp

    my ($pdbid, $pdbpos, $chain) =
        ($hsp->{PDB_id}, $hsp->{PDB_pos}, $hsp->{PDB_chain});

    my $pdbid_short = lc $pdbid;
    my $dssp_file   = $pdbid_short . '.dssp';
    my $dssp_filez  = $dssp_file . '.Z';

    # If uncompressed file is found in the local DSSP database then file
    # locking is not required since local database mirror is read-only.
    unless (PPH::Misc::locate_file($dssp_file, $CONFIG{'DSSP'})) {

      $CONFIG{'SAVEPATH'} = $CONFIG{'SAVEKNOWN'} && $PROT->{UniProt} ? $CONFIG{'PRECOMPATH'} : $CONFIG{'SCRATCH'};
      my $lockfile = "$CONFIG{SAVEPATH}/structures/$dssp_file";
      open my $flk, '>>', $lockfile or die "ERROR: Can't create lockfile $lockfile: $!\n";
      flock($flk, LOCK_EX) or die "ERROR: Failed to obtain lock for $lockfile: $!\n";
      ####### |     file locked: $lockfile

      #----------------------------------------
      # Already computed (uncompressed)
      #----------------------------------------
      if    (PPH::Misc::locate_file($dssp_file,
             "$CONFIG{PRECOMPATH}/structures", "$CONFIG{SCRATCH}/structures")) {
        # nothing to do, $dssp_file is set to go
      }
      #----------------------------------------
      # Already computed (compessed)
      #----------------------------------------
      elsif (PPH::Misc::locate_file($dssp_filez,
             $CONFIG{'DSSP'}, "$CONFIG{PRECOMPATH}/structures", "$CONFIG{SCRATCH}/structures")) {
        my $rc = system("gunzip -c $dssp_filez >$dssp_file 2>/dev/null");
        warn("WARNING: get_DSSP: gunzip failed for: $dssp_filez\n"), close($flk), return if $rc;
      }
      #----------------------------------------
      # Needs to be computed
      #----------------------------------------
      else {
        my $pdbfile = _get_pdb_file($pdbid);
        warn("WARNING: get_DSSP: File can not be found nor calculated for: $pdbid\n"), close($flk), return
          unless defined $pdbfile;

        warn("WARNING: get_DSSP: $CONFIG{DSSPCMBI} failed for: $pdbfile, $dssp_file\n"), close($flk), return
          if system "$CONFIG{DSSPCMBI} $pdbfile $dssp_file >/dev/null 2>&1";

        copy($dssp_file, "$CONFIG{SAVEPATH}/structures/")
          or cluck "ERROR: get_DSSP: Can't save file: $dssp_file";
      }

      close($flk);
    }

    warn("WARNING: get_DSSP: File empty: $dssp_file\n"), return unless -s $dssp_file;

    ##### |    open: $dssp_file
    open my $F, '<', $dssp_file
        or do { warn "WARNING: get_DSSP: Can't open file: $dssp_file\n"; return };

    #----------------------------------------
    # Ignore Header
    #----------------------------------------
    while (<$F>) { /RESIDUE\s+AA\s+STRUCTURE\s+BP1\s+BP2\s+ACC/ and last }

    #----------------------------------------
    # Now parse DSSP - use substr bc. not
    # all lines contain all elements!
    # For details of the DSSP v3 format, see:
    #   http://swift.cmbi.ru.nl/gv/dssp/DSSP_3.html
    #----------------------------------------
    #01234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789
    #     |     | |  |                  |                                                                   |     |
    #     RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA
    #11721   51 Z D  S    S+     0   0  119     -2,-0.5     2,-0.3     1,-0.2    -1,-0.1   0.723  73.9  10.2 -83.1 -22.3   27.4   -4.2    4.8
    #11722   52 Z S        -     0   0    6    -15,-0.2     2,-0.3    13,-0.0   -15,-0.2  -0.986  56.2-152.0-158.0 158.4   25.5   -6.6    2.7
    #11723   53 Z Y  E     -P 11706   0l  19    -17,-2.5   -17,-2.7    -2,-0.3     2,-0.3  -0.982  14.4-164.1-138.5 148.2   24.1  -10.1    2.5
    #11724   54 Z V  E >   -P 11705   0l  14     -2,-0.3     3,-1.0   -19,-0.2     4,-0.2  -0.988  48.8 -17.4-139.8 144.9   23.3  -12.2   -0.7
    #1tzn.dssp:
    #10137 1038 o S              0   0  163      0, 0.0     2,-0.3     0, 0.0   180,-0.0   0.000 360.0 360.0 360.0-122.2  -24.3   39.3   81.3
    #10138 1039 o n        -     0   0   50      2,-0.0   178,-0.3   179,-0.0     2,-0.3  -0.708 360.0-170.3-106.6 159.0  -23.6   43.1   81.3
    #10139 1040 o R  B     -t 10316   0S  89    176,-2.4   178,-1.0    -2,-0.3     3,-0.1  -0.844  21.8-162.3-138.8 174.2  -20.5   45.0   80.3
    #10140 1041 o R  S    S+     0   0  140     -2,-0.3     2,-0.8   176,-0.2    36,-0.2   0.581  71.5  76.7-130.1 -32.3  -19.0   48.5   80.4
    #10141 1042 o A        +     0   0   46     36,-0.1     2,-0.4    34,-0.1    98,-0.3  -0.774  57.3 173.1 -92.8 109.2  -16.1   48.5   78.0
    #10142 1043 o F        -     0   0    2     -2,-0.8    36,-0.6    34,-0.2     2,-0.6  -0.893  21.4-153.0-116.2 145.4  -17.2   48.7   74.4
    #10143 1044 o D  E     -uv1017810240T   4     96,-1.5    98,-2.8    -2,-0.4     2,-0.5  -0.948  13.5-158.3-119.3 109.0  -15.0   49.1   71.3
    #10144 1045 o L  E     -uv1017910241T   0     34,-1.3    36,-1.6    -2,-0.6     2,-0.5  -0.772   1.5-160.2 -93.8 127.0  -16.8   50.7   68.3
    #10145 1046 o Y  E     -uv1018010242T   0     96,-1.9    98,-3.6    -2,-0.5     2,-0.8  -0.912   2.7-158.0-106.7 124.9  -15.4   50.2   64.8
    #10146 1047 o F  E     -uv1018110243T   0     34,-3.4    36,-3.2    -2,-0.5     2,-0.7  -0.707  10.8-177.9-107.3  83.1  -16.5   52.7   62.2
    #10147 1048 o V  E     -uv1018210244T   0     96,-1.2    98,-3.7    -2,-0.8     2,-0.5  -0.698  14.9-164.1 -79.4 116.0  -16.1   51.1   58.8
    #10148 1049 o L  E     -uv1018310245T   0     34,-3.8    36,-2.9    -2,-0.7     2,-1.0  -0.911  19.4-136.8-113.8 129.1  -17.0   53.8   56.2
    #10149 1050 o D  E     +u 10184   0T   0     96,-2.7    36,-0.2    -2,-0.5    12,-0.1  -0.678  24.9 174.1 -79.8 103.4  -17.8   53.3   52.6
    #10150 1051 o K        +     0   0   42     34,-3.1    64,-2.3    -2,-1.0    36,-0.1  -0.153  31.4 135.4-103.2  39.6  -16.0   56.1   50.8
    #10151 1052 o S  S >  S-     0   0    0     62,-0.3     3,-1.5    34,-0.1    65,-0.2  -0.468  70.2-109.7 -85.2 158.4  -16.9   54.8   47.3
    #----------------------------------------
    while (<$F>) {
        #----------------------------------------
        my $ntmp   = substr($_, 5,6); $ntmp =~ s/ //g; # PDB Residue name
        my $ctmp   = substr($_,11,1);                   # PDB Chain ID
        $ctmp eq ' '     and $ctmp =  $CONFIG{ONECHAIN};
        $ntmp eq $pdbpos and $ctmp eq $chain or next;
        #----------------------------------------

        #----------------------------------------
        my $aa     = substr($_,13,1);  $aa  =~tr/a-z/C/; # Cystins to C (?StS?)
        my $secstr = substr($_,16,1);                     # Secondary Structure
        ($secstr eq ' ') and $secstr = $CONFIG{NOSECSTR};
        #----------------------------------------

        # BP1+BP2 fields can overflow in case of extra long chains: check if
        # fixed line length exceeds standard 136 columns + newline and offset
        # all fields beyond BP2 by the amount of extra columns found.
        my $len    = length $_;
        my $offset = $len > 137 ? $len - 137 : 0;

        #----------------------------------------
        my $acc    = substr($_, 35+$offset,3); $acc =~ s/ //g; # Water-accessible surface area
        my $phi    = substr($_,103+$offset,6); $phi =~ s/ //g; # PHI
        my $psi    = substr($_,109+$offset,6); $psi =~ s/ //g; # PSI
        #----------------------------------------

        #----------------------------------------
        my $noracc = ($acc / $PPH::Data::MaxAccessibility{$aa1});
        my $propi  = _Propensity_interval($noracc);
        my $dprop  = abs($PPH::Data::Propensities{$aa1}->[$propi] -
                             $PPH::Data::Propensities{$aa2}->[$propi]);
        my $dvol   = int($PPH::Data::AA_Volume{$aa2} -
                             $PPH::Data::AA_Volume{$aa1});
        #----------------------------------------

        close($F);
        my $mapreg = _map_PhiPsi($phi, $psi);
        ##### |    DSSP returns: $aa, $secstr, $acc, $noracc, $dprop, $dvol, $phi, $psi, $mapreg
        return {
          DSSPaa  => $aa,
          SecStr  => $secstr,
          ASA     => $acc,
          NormASA => $noracc,
          dProp   => $dprop,
          dVol    => $dvol,
          Phi     => $phi,
          Psi     => $psi,
          MapReg  => $mapreg
        };
    }
    close($F);
    ##### |    DSSP returns undef for: $pdbid, $pdbpos, $chain
    return;
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head3 _map_PhiPsi

 Usage   : _map_PhiPsi($phi, $psi)
 Function: Map the given phi psi angles to areas
 Args    : Phi & Psi angles
 Returns : Mapped region from PPH::Data::PhiPsiMap

=cut

#----------------------------------------
sub _map_PhiPsi {
    my ($phi, $psi) = @_;

  (($phi == 360) or ($psi == 360)) and return '?';

  SWITCH: {
        ($phi >=  179.9) and do { $phi =  179; last SWITCH};
        ($psi >=  179.9) and do { $psi =  179; last SWITCH};
        ($phi <= -179.9) and do { $phi = -179; last SWITCH};
        ($psi <= -179.9) and do { $psi = -179; last SWITCH};
    }
    my $CELL = 10;		# Grads in one Map Cell ?StS?
    my $i    = int((180 + $phi)/$CELL);
    my $j    = int((180 - $psi)/$CELL);

    return substr($PPH::Data::PhiPsiMap[$j], $i, 1);


}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head3 _Propensity_interval

 Usage   : _Propensity_interval
 Function: Return propensisty interval
 Args    : Accessibility
 Returns : Interval 0-6

=cut

#----------------------------------------
sub _Propensity_interval {
    my $Accessibility = @_;

    ##### PPH__3D__Propensity_interval...

    ($Accessibility == 0 ) and return 0;
    ($Accessibility  < 5 ) and return 1;
    ($Accessibility  < 15) and return 2;
    ($Accessibility  < 30) and return 3;
    ($Accessibility  < 50) and return 4;
    ($Accessibility  < 75) and return 5;
    return 6;
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head2 get_HBonds

 Usage   : PPH::3D->get_HBonds($hsp)
 Function: Check if SNP is involved in an hbond
 Args    : Object of type HSP (XML BLAST)
 Returns : Array containing bond information

=cut

#----------------------------------------
sub get_HBonds {
    my ($PROT, $hsp) = @_;
    my ($PDB_id, $Qpos, $Qchain) =
        ($hsp->{PDB_id}, $hsp->{PDB_pos}, $hsp->{PDB_chain});
    #### PPH__3D__get_HBonds: $PDB_id, $Qpos, $Qchain

    my $hb2file = lc($PDB_id) . '.hb2';

    $CONFIG{'SAVEPATH'} = $CONFIG{'SAVEKNOWN'} && $PROT->{UniProt} ? $CONFIG{'PRECOMPATH'} : $CONFIG{'SCRATCH'};
    my $lockfile = "$CONFIG{SAVEPATH}/structures/$hb2file";
    open my $flk, '>>', $lockfile or die "ERROR: Can't create lockfile $lockfile: $!\n";
    flock($flk, LOCK_EX) or die "ERROR: Failed to obtain lock for $lockfile: $!\n";
    ####### |     file locked: $lockfile

    #----------------------------------------
    # Creating unless the hb2 file exists
    #----------------------------------------
    unless (PPH::Misc::locate_file($hb2file,
      "$CONFIG{PRECOMPATH}/structures", "$CONFIG{SCRATCH}/structures")) {
        #----------------------------------------
        # Find pdb file
        #----------------------------------------
        my $PDBid_short = $PDB_id;
        my $file        = _get_pdb_file($PDB_id);
        defined $file and -e $file or do {
          warn "WARNING: get_HBonds: Can't find PDB file for: $PDB_id\n";
          close($flk);
          return;
        };
        #----------------------------------------
        #----------------------------------------
        # Execute hbplus
        #----------------------------------------
        ##### hbplus: $file
        open my $EX, '-|', "$CONFIG{HBPLUS} $file 2>&1" or do {
          unlink $hb2file;
          confess "ERROR: get_HBonds: $CONFIG{HBPLUS} failed for: $file";
        };
        unless (grep /hydrogen bonds found/, do { local $/; <$EX> }) {
          warn "WARNING: get_HBonds: $CONFIG{HBPLUS} failed to locate any hydrogen bonds for: $PDB_id\n";
          close($flk);
          return;
        }
        close $EX;
        # Now save created hb2 file
        copy($hb2file, "$CONFIG{SAVEPATH}/structures/")
          or cluck "WARNING: get_HBonds: Can't save file: $hb2file";
    }
    #----------------------------------------
    close($flk);

    #----------------------------------------
    # Now parse
    # Example:
    # <---DONOR---> <-ACCEPTOR-->    atom                        ^
    # c    i                          cat <-CA-CA->   ^        H-A-AA   ^      H-
    # h    n   atom  resd res      DA  || num        DHA   H-A  angle D-A-AA Bond
    # n    s   type  num  typ     dist DA aas  dist angle  dist       angle   num
    #
    # D0001-MET N    -0035-HOH O   2.62 MH  -2 -1.00 135.3  1.80  -1.0  -1.0     1
    #----------------------------------------

    warn("WARNING: get_HBonds: File empty: $hb2file\n"), return unless -s $hb2file;

    ##### |    open: $hb2file
    open my $HB, '<', $hb2file or do {
      warn "WARNING: get_HBonds: Can't open file: $hb2file\n";
      return;
    };

    # Ignore header
    while (<$HB>) { /angle\s+num\s*$/ and last; }

    ##### |    Bonds within AA allowed? Then no switch!

    #----------------------------------------
    # If the Qchain is only one chain DSSP
    # contains as sign '_' which has to be
    # converted to '-'
    #----------------------------------------
    $Qchain eq $CONFIG{'ONECHAIN'} and $Qchain = '-';
    #----------------------------------------

    #----------------------------------------
    # Now extract the relevant hbonds
    #----------------------------------------
    my @Bonds;
    while (my $l = <$HB>) {

      #----------------------------------------
      # Check if hbplus is ok
      #----------------------------------------
      length($l) != 76 and do {
        warn "WARNING: get_HBonds: File $hb2file seems to have incorrect format:\n$l\n";
        return;
      };
      #----------------------------------------

      #----------------------------------------
      # Extract the relevant info - split won't
      #                             work here
      #----------------------------------------
      my ($chain1,$pos1,$aa1,$atom1) = (substr($l, 0,1), substr($l, 1,4),
                                        substr($l, 6,3), substr($l,10,3));
      my ($chain2,$pos2,$aa2,$atom2) = (substr($l,14,1), substr($l,15,4),
                                        substr($l,20,3), substr($l,24,3));
      my ($Distance, $Type)          = (substr($l,28,4), substr($l,33,2));
      #----------------------------------------

      SWITCH: {
            ($Qchain eq $chain1 and $Qpos eq $pos1)
                and ($Type eq 'SM' or $Type eq 'SS')
                    and do {
                        push @Bonds, "$atom1 -> $chain2 $pos2 $aa2 $atom2 ".
                            "$Distance $Type";
                        last SWITCH;
                    };
            ($Qchain eq $chain2 and $Qpos eq $pos2)
                and ($Type eq 'MS' or $Type eq 'SS')
                    and do {
                        push @Bonds, "$atom2 -> $chain1 $pos1 $aa1 $atom1 ".
                            "$Distance $Type";
                        last SWITCH;
                    };
      }
    }
    close $HB;
    ##### |    get_HBonds returns (n hbonds): scalar(@Bonds)

    return @Bonds;
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head2 get_BFactor

 Usage   : PPH::3D->get_BFactor($hsp)
 Function: Calculate the B-factor for a certain position
 Args    : Object of type HSP (XML BLAST)
 Returns : Returns the output of the C-Program bfact

=cut

#----------------------------------------
sub get_BFactor {
    my ($hsp) = @_;
    my ($PDB_id, $AApos, $chain) =
        ($hsp->{PDB_id}, $hsp->{PDB_pos}, $hsp->{PDB_chain});
    #### PPH__3D__get_BFactor: $PDB_id, $AApos, $chain

    my $PDBid_short = $hsp->{PDB_id};
    my $file        = _get_pdb_file($hsp->{PDB_id});
    defined $file or do { warn "WARNING: get_BFactor: Can't find file for PDB entry: $PDB_id\n"; return ''; };

    my $ret = `$CONFIG{BFACT} $file $AApos $chain 2>&1`;

    ($ret =~ /NormedB:\s+(\S+)\b/ && $ret !~ /Error/) and return $1;

    chomp $ret; $ret =~ s/\n/ /; $ret =~ s/\s+/ /;
    warn "WARNING: bfact($file, $AApos, $chain) failed: $ret\n";

    return '';
}
#----------------------------------------------------------------------

#------------------------------------------------------------------------------------------

=head2 pos2pdb

 Usage   : PPH::3D->pos2pdb($id, $chain, $pos)
 Function: Map a sequence position to PDB structure;
           this is necessary since the structure
           often lacks N-terminus and other parts of
           the sequence
 Args    : PDB structure identifier
           PDB chain ID
           FASTA sequence position
 Returns : PDB postion (residue identifier) if successful

=cut

#------------------------------------------------------------------------------------------

#   Table: fragm
#   fstart / fend correspond to sequence ranges or single positions in FASTA sequence
#   astart / aend correspond to ranges or single residue codes in matching PDB ATOM records
#   scale is 1-based for ranges / positions
#
#  pdbid  chain   fstart    fend   astart    aend
#   102m      A        1     154        0     153
#   103l      A        1      34        1      34
#   103l      A       35      35      40A     40A
#   103l      A       36      36      40B     40B
#   103l      A       37      37      40C     40C
#   103l      A       38     159       41     162
#
{
my $dbfile = "$CONFIG{PDB2FASTA}/$CONFIG{PDB_ID_SUFFIX}.sqlite";
my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
my $select = q{
SELECT CASE
  WHEN fstart < fend AND ?3 BETWEEN fstart AND fend
  THEN ?3 - fstart + astart
  WHEN fstart = fend AND fstart = ?3
  THEN astart
END AS pdbpos
FROM fragm
WHERE pdbid = ?1 AND chain = ?2 AND pdbpos NOT NULL
};
my ($dbh, $sth);
sub pos2pdb {
  my ($id, $chain, $pos) = @_;
  #### PPH__3D__pos2pdb: $id, $chain, $pos
  # Connect to database and prepare select statement if this has not been done yet
  unless (defined $dbh) {
    $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs);
    ####### |     SQL statement: $select
    $sth = $dbh->prepare($select);
  }
  # Search table
  ###### |     pdb2pos executes: $id, $chain, $pos
  $sth->execute($id, $chain, $pos);
  # Fetch results (should be just a single row)
  my $pdbpos = ($sth->fetchrow_array)[0];
  $sth->finish;
  confess "ERROR: pos2pdb: Failed to locate sequence position ($pos) in PDB entry: $id$chain"
    unless defined $pdbpos && length $pdbpos;
  ###### |     pdb2pos returns: $pdbpos
  return $pdbpos;
}
}

#------------------------------------------------------------------------------------------

#----------------------------------------------------------------------

=head2 get_contacts

 Usage   : PPH::3D->get_contacts($hsp, $intra_chain_flag)
 Function: Find Contacts in structure for given position
 Args    : Object of type HSP (XML BLAST)
           Argument for C-Program dist wether to check for intra chain contacts
 Returns : Hashref with different contact types
           { INT | EXT => [ { Pos, Name, Dist} ] }

=cut

#----------------------------------------
sub get_contacts {
    my ($hsp, $intra_chain_flag) = @_;
    my ($PDB_id, $pos, $chain) =
        ($hsp->{PDB_id}, $hsp->{PDB_pos}, $hsp->{PDB_chain});
    #### PPH__3D__get_contatcs (PDB_id, pos, chain, intra_chain_flag): $PDB_id, $pos, $chain, $intra_chain_flag

    my $file = _get_pdb_file($PDB_id);
    defined($file) or do {
      warn "WARNING: PDB file '$PDB_id' can't be found"; return {}
    };

    my $PDBid_short = $PDB_id;

    #----------------------------------------

    #----------------------------------------
    # Most time chain is in uppercase
    # - there are structure such as 1d3i which
    #   contain chains with lower case:
    # - so first try chain and then lc of it
    #----------------------------------------

    my $flag = $intra_chain_flag ? ' -i' : '';

    my $exec = $CONFIG{DIST}."$flag -r $pos ".uc($chain).' '.$CONFIG{CONTTHRESH}." $file 2>&1";
    ##### |    exec: $exec
    open my $DIST, '-|', $exec or do {
        ### Error executing dist: $exec
        return {}
    };

    my $ret = do { local $/; <$DIST> };
    if (grep /Error|Usage/, $ret) {
        close $DIST;
        $exec = $CONFIG{DIST}."$flag -r $pos ".lc($chain).' '.$CONFIG{CONTTHRESH}." $file 2>&1";
        ##### |    retry exec: $exec
        open $DIST, '-|', $exec or do {
            ### Error executing dist: $exec
            return {}
        };
        $ret = do { local $/; <$DIST> };
        if (grep /Error|Usage/, $ret) {
            close $DIST;
            ### Error running dist: $exec
            return {}
        };
    }

    close $DIST;

    # 1thf   INT     1D MET  O  4     m >...< m     3D ALA  N  17       3.538
    my %ret;
    foreach (split /^/, $ret) {
      my ($type, $atom1, $pos2, $aa2, $atom2, $dist) =
          (split /\s+/, $_)[1,4,9,10,11,13];
      my $chain2;
      ($pos2, $chain2) = $pos2 =~ /(^-?\d+)(\w)/;
      push @{ $ret{$type} }, {
        Pos   => $pos2,
        Chain => $chain2,
        Atom1 => $atom1,
        Atom2 => $atom2,
        Name  => $aa2,
        Dist  => $dist
      };
    }
    ##### |    contacts returns: %ret
    return \%ret;
}
#----------------------------------------------------------------------

#----------------------------------------------------------------------

=head3 _get_pdb_file

 Usage   : _get_pdb_file($pdbid)
 Function: returns PDB filename
           - it unzips PDB if necessary and marks this file as GARBAGE
 Args    : PDB identifier
 Returns : filename

=cut

#----------------------------------------
sub _get_pdb_file {
  # pdbid used in PDB filenames is always lowercase
  my $pdbid = lc shift;
  ##### PPH__3D__get_pdb_file: $pdbid

  my $pdbfile   = "$pdbid.pdb";
  my $wwpdbfile = "$CONFIG{PDB}/$CONFIG{PDB_PREFIX}$pdbid.$CONFIG{PDB_SUFFIX}";
  ####### |     lookup wwpdb file: $wwpdbfile

  # Return if we have uncompressed file in the current working directory
  return $pdbfile if -e $pdbfile;

  # Return if we have an uncompressed file available in the local PDB mirror
  return $wwpdbfile if -e $wwpdbfile;

  # Look up compressed file in the local PDB mirror
  my ($zipped, $zfound);
  foreach my $ext (qw(.gz .Z)) { $zipped = $wwpdbfile.$ext; $zfound = 1, last if -e $zipped; }
  warn("WARNING: _get_pdb_file: File was not found for PDB entry: $pdbid\n"), return unless $zfound;

  # Create uncompressed copy of the PDB file in the current directory
  my $status = system "gunzip -c $zipped >$pdbfile 2>/dev/null";
  warn("WARNING: _get_pdb_file: Failed to gunzip file: $zipped\n"), return unless $status == 0;

  return $pdbfile;
}

#----------------------------------------------------------------------

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
1;
