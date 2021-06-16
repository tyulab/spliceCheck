package PPH::Rules;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Rules

=head1 DESCRIPTION

Return a prediction for the effect of a non-synonymous substitution -
rule based

=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2011-09-20 21:30:07 -0400 (Tue, 20 Sep 2011) $
 $LastChangedRevision: 368 $
 $LastChangedBy: ivan $



=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


use Carp qw(cluck);
use Data::Dumper;
use strict;
use warnings;

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

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR


=head1 ROUTINES

=head2 test

 Usage   : $PPH::Rules->test(PROT, SNP)
 Function: Rule based prediction of the effect of a non-synonymous
           substitution
 Args    : Objects of type PROT and SNP
 Returns : Hashref of prediction
           {
             Prediction
             Avail
             Basis
             Effect
             Data
           }

=cut

#----------------------------------------------------------------------
sub test {
  my ($PROT, $SNP) = @_;
  ### PPH__Rules__test...

  my $features;
  if (keys %{ $SNP->{Features} }) {
    $features->{Site}   = $SNP->{Features}{Sites}[0]   if exists $SNP->{Features}{Sites}   && @{ $SNP->{Features}{Sites} };
    $features->{Region} = $SNP->{Features}{Regions}[0] if exists $SNP->{Features}{Regions} && @{ $SNP->{Features}{Regions} };
    $features->{PHAT}   = $SNP->{Features}{PHAT}       if exists $SNP->{Features}{PHAT};
  }
  my $selected  = $SNP->{Scores}{Selected};
  my $alignment = 'alignment' if length $selected;
  $alignment   .= '_mz' if defined $alignment && $selected eq 'Mz';
  my $delta     = length $selected && defined $SNP->{Scores}{$selected}{PsicD} ?
                  abs($SNP->{Scores}{$selected}{PsicD}) : undef;
  my $structure = exists $SNP->{Structure} && keys %{ $SNP->{Structure} } && $SNP->{Structure}{HitsMapped} ?
                  $SNP->{Structure}{HitsMapped} : undef;
  my $struct_params  = defined $structure ? $SNP->{Structure}{Maps}[0]{Params} : undef;
  my $struct_contact = defined $structure && exists $SNP->{Structure}{BestContact} ?
                       $SNP->{Structure}{BestContact} : undef;

  my @avail;
  defined $features  and push @avail, 'FT';
  defined $delta     and push @avail, 'alignment';
  defined $structure and push @avail, 'structure';

  #--------------------------------------------------------------------
  # Worst case : no prediction is possible ('unknown')
  #--------------------------------------------------------------------
  scalar @avail or return
    {
      Prediction => $PPH::Data::Predictions[0],
      Avail      => undef,
      Basis      => undef,
      Effect     => undef,
      Data       => undef
    };
  #---------------------------------------------------------------------

  #---------------------------------------------------------------------
  # We have a UniProt site features
  #---------------------------------------------------------------------
  if (defined $features) {
    #### features defined...
    #----------------------------------------
    # Site features
    #----------------------------------------
    exists $features->{Site}
      and ($features->{Site} ne 'CARBOHYD')
        and return
          {
            Prediction => $PPH::Data::Predictions[3],
            Avail      => [ @avail ],
            Basis      => 'sequence annotation',
            Effect     => grep( { /$features->{Site}/ } keys %PPH::Data::BONDS ) ? '1.2' : '2.2',
            Data       => 'Site type: ' . $features->{Site}
          };
    #----------------------------------------

    #----------------------------------------
    # Region features
    #----------------------------------------
    exists $features->{Region}
      and $features->{Region} eq 'TRANSMEM'
        and return
          {
            Prediction => $features->{PHAT} < 0 ?
                          $PPH::Data::Predictions[2] : $PPH::Data::Predictions[1],
            Avail      => [ @avail ],
            Basis      => exists $features->{RegionStatus} &&
                          defined $features->{RegionStatus} &&
                          $features->{RegionStatus} eq 'predicted' ?
                            'sequence prediction' :
                            'sequence annotation',
            Effect     => $features->{PHAT} < 0 ? '2.2.2' : undef,
            Data       => "PHAT matrix element difference: " . $features->{PHAT}
          };
    #----------------------------------------
  }
  #---------------------------------------------------------------------

  #---------------------------------------------------------------------
  # if nothing else...
  defined $delta
    or return
      {
        Prediction => $PPH::Data::Predictions[0],
        Avail      => [ @avail ],
        Basis      => undef,
        Effect     => undef,
        Data       => undef
      };
  #---------------------------------------------------------------------

  #---------------------------------------------------------------------
  # Structure - contacts to important residues or substrates / cofactors
  #---------------------------------------------------------------------
  if ($delta > 1 and defined $structure and defined $struct_contact) {
    ###### Delta lt 1 (delta, structure, struct_contact): $delta, $structure, $struct_contact
    #----------------------------------------
    # A) Heteroatom (e.g.substrate / cofactor)
    #----------------------------------------
    exists $struct_contact->{Heteroatoms}        and
     exists $struct_contact->{Heteroatoms}{Dist} and
     $struct_contact->{Heteroatoms}{Dist} <= 3   and
      return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '2.2.3',
          Data       => "PSIC score difference: $delta; ligand name: $struct_contact->{Heteroatoms}{Name}, distance: $struct_contact->{Heteroatoms}{Dist} A"
        };
    #----------------------------------------

    #----------------------------------------
    # B) Critical site (e.g. catalytical site)
    #----------------------------------------
    exists $struct_contact->{CritSites}       and
    exists $struct_contact->{CritSites}{Dist} and
    $struct_contact->{CritSites}{Dist} <= 3   and
      return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '2.1',
          Data       => "PSIC score difference: $delta; functional site: $struct_contact->{CritSites}{Name} ($struct_contact->{Number}{CritSites}), distance: $struct_contact->{CritSites}{Dist} A"
        };
    #----------------------------------------
  }
  #---------------------------------------------------------------------

  #---------------------------------------------------------------------
  # PSIC - sequence conservation plus structure parameters if available
  #---------------------------------------------------------------------
  $delta <= 0.5 and return
    {
      Prediction => $PPH::Data::Predictions[1],
      Avail      => [ @avail ],
      Basis      => $alignment,
      Effect     => undef,
      Data       => "PSIC score difference: $delta"
    };
  #----------------------------------------

  #----------------------------------------
  $delta > 0.5 and $delta <= 1.5 and do {
    ###### Delta gt 0.5 and le 1.5 (delta, NormASA, dProp, dVol): $delta, $struct_params->{NormASA}, $struct_params->{dProp}, $struct_params->{dVol}
    #----------------------------------------
    # Is SNP on surface?
    #----------------------------------------
    if (
    defined $structure               and
    exists $struct_params->{NormASA} and
    $struct_params->{NormASA} <= 0.05) {
      #----------------------------------------
      # Propensity changes of aa variants
      #----------------------------------------
      $struct_params->{dProp} >= 1 and return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '1.1.1',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; hydrophobicity change: $struct_params->{dProp}"
        };
      #----------------------------------------

      #----------------------------------------
      # Vol changes of aa variants
      #----------------------------------------
      $struct_params->{dVol} >= 80 and return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '1.1.2',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; volume change: $struct_params->{dVol}"
        };

      $struct_params->{dVol} <= -80 and return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [@avail],
          Basis      => 'structure',
          Effect     => '1.1.3',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; volume change: $struct_params->{dVol}"
        };
      #----------------------------------------
    }
    #----------------------------------------

    #----------------------------------------
    # Lower threshold for SNP accesibility
    #----------------------------------------
    if (
    defined $structure                and
    exists $struct_params->{NormASA}  and
    $struct_params->{NormASA} >  0.05 and
    $struct_params->{NormASA} <= 0.15) {
      #----------------------------------------
      # Propensity changes
      #----------------------------------------
      $struct_params->{dProp} >= 0.75 and return
        {
          Prediction => $PPH::Data::Predictions[2],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '1.1.1',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; hydrophobicity change: $struct_params->{dProp}"
        };
      #----------------------------------------

      #----------------------------------------
      # Vol changes
      #----------------------------------------
      $struct_params->{dVol} >= 60 and return
        {
          Prediction => $PPH::Data::Predictions[2],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '1.1.2',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; volume change: $struct_params->{dVol}"
        };

      $struct_params->{dVol} <= -60 and return
        {
          Prediction => $PPH::Data::Predictions[2],
          Avail      => [@avail],
          Basis      => 'structure',
          Effect     => '1.1.3',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; volume change: $struct_params->{dVol}"
        };
    }
    #----------------------------------------

    return {
      Prediction => $PPH::Data::Predictions[1],
      Avail      => [ @avail ],
      Basis      => $alignment,
      Effect     => undef,
      Data       => "PSIC score difference: $delta"
    };

  };  # end ($delta > 0.5 and $delta <= 1.5)
  #----------------------------------------

  #----------------------------------------
  $delta > 1.5 and $delta <= 2.0 and do {
    ###### Delta gt 1.5 and le 2.0 (delta, NormASA, dProp, dVol): $delta, $struct_params->{NormASA}, $struct_params->{dProp}, $struct_params->{dVol}

    #----------------------------------------
    # Accessbility
    #----------------------------------------
    if (
    defined $structure                and
    exists $struct_params->{NormASA}  and
    $struct_params->{NormASA} <= 0.05) {
      #----------------------------------------
      # Propensity
      #----------------------------------------
      $struct_params->{dProp} >= 1 and return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '1.1.1',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; hydrophobicity change: $struct_params->{dProp}"
        };
      #----------------------------------------
      # Vol change
      #----------------------------------------
      $struct_params->{dVol} >= 80 and return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [ @avail ],
          Basis      => 'structure',
          Effect     => '1.1.2',
          Data       => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; volume change: $struct_params->{dVol}"
        };

      $struct_params->{dVol} <= -80 and return
        {
          Prediction => $PPH::Data::Predictions[3],
          Avail      => [@avail],
          Basis      => 'structure',
          Effect     => '1.1.3',
          Data	     => "PSIC score difference: $delta; normed accessibility: $struct_params->{NormASA}; volume change: $struct_params->{dVol}"
        };
    }
    #----------------------------------------

    return
        {
          Prediction => $PPH::Data::Predictions[2],
          Avail      => [@avail],
          Basis      => $alignment,
          Effect     => undef,
          Data       => "PSIC score difference: $delta"
        };
  };  # end ($delta > 1.5 and $delta <= 2.0)
  #------------------------------

  #------------------------------
  $delta > 2.0 and return
    {
      Prediction  => $PPH::Data::Predictions[3],
      Avail       => [ @avail ],
      Basis       => $alignment,
      Effect      => undef,
      Data        => "PSIC score difference: $delta"
    };
  #------------------------------

  #------------------------------
  #
  warn "ERROR: PPH::Rules::test: Should not have happened -- all prediction rules failed!\n";
  return
    {
      Prediction => $PPH::Data::Predictions[0],
      Avail      => [ @avail ],
      Basis      => undef,
      Effect     => undef,
      Data       => undef
    };
  #------------------------------
}
#----------------------------------------------------------------------

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
1;
