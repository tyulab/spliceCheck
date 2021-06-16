#!/usr/bin/env perl

use Getopt::Std;
use File::Basename;

# PPHv2 features (50 total) vs. WEKA features
# (Note: PPHv2 input Pos is zero-based, WEKA output Pos are 1-based)
#
#       PPHv2                   WEKA
#  --------------- -----------------------------------
#  Pos  Name        Name              Pos33 Pos11 Pos8
#  --------------- -----------------------------------
#    0  o_acc
#    1  o_pos
#    2  o_aa1
#    3  o_aa2
#    4  rsid
#    5  acc
#    6  pos
#    7  aa1
#    8  aa2
#    9  nt1         nt1                   1
#   10  nt2         nt2                   2
#   11  prediction
#   12  based_on
#   13  effect
#   14  site        site                  3
#   15  region      region                4
#   16  PHAT        phat                  5
#   17  dScore      score_delta           6     1    1
#   18  Score1      score1                7     2    2
#   19  Score2      score2                8
#   20  MSAv        msa_ver
#   21  Nobs        num_obs               9     3    3
#   22  Nstruct
#   23  Nfilt
#   24  PDB_id
#   25  PDB_pos
#   26  PDB_ch
#   27  ident       ali_ide              21
#   28  length      ali_len              22
#   29  NormASA     acc_normed           23     9
#   30  SecStr      sec_str              24
#   31  MapReg      map_region           25
#   32  dVol        delta_volume         10     4
#   33  dProp       delta_prop           26
#   34  B-fact      b_fact               27    10
#   35  H-bonds
#   36  AveNHet     het_cont_ave_num     28
#   37  MinDHet     het_cont_min_dist    29
#   38  AveNInt     inter_cont_ave_num   30
#   39  MinDInt     inter_cont_min_dist  31
#   40  AveNSit
#   41  MinDSit
#   42  Transv      transversion         11
#   43  CodPos
#   44  CpG         cpg                  12
#   45  MinDJnc
#   46  PfamHit     pfam_hit             13     5    4
#   47  IdPmax      id_p_max             14     6    5
#   48  IdPSNP      id_p_snp             15
#   49  IdQmin      id_q_min             16     7    6
#
# The following 6 extra features have no direct equivalents in PPH2 output:
#
#   50  cpg_var1_var2                    17
#   51  cpg_transition                   18     8    7
#   52  charge_change                    19
#   53  hydroph_change                   20
#   54  delta_volume_new                 32
#   55  delta_prop_new                   33    11    8

# Same mappings as above encoded as array of arrays
my @PPH2WEKA = (
  [ 'o_acc',      '',                     '',           0,  0,  0 ],
  [ 'o_pos',      '',                     '',           0,  0,  0 ],
  [ 'o_aa1',      '',                     '',           0,  0,  0 ],
  [ 'o_aa2',      '',                     '',           0,  0,  0 ],
  [ 'snp_id',     '',                     '',           0,  0,  0 ],
  [ 'acc',        '',                     '',           0,  0,  0 ],
  [ 'pos',        '',                     '',           0,  0,  0 ],
  [ 'aa1',        '',                     '',           0,  0,  0 ],
  [ 'aa2',        '',                     '',           0,  0,  0 ],
  [ 'nt1',        'nt1',                  '',           1,  0,  0 ],
  [ 'nt2',        'nt2',                  '',           2,  0,  0 ],
  [ 'prediction', '',                     '',           0,  0,  0 ],
  [ 'based_on',   '',                     '',           0,  0,  0 ],
  [ 'effect',     '',                     '',           0,  0,  0 ],
  [ 'site',       'site',                 'anyyes',     3,  0,  0 ],
  [ 'region',     'region',               '',           4,  0,  0 ],
  [ 'PHAT',       'phat',                 'numeric',    5,  0,  0 ],
  [ 'dScore',     'score_delta',          'numeric',    6,  1,  1 ],
  [ 'Score1',     'score1',               'numeric',    7,  2,  2 ],
  [ 'Score2',     'score2',               'numeric',    8,  0,  0 ],
  [ 'MSAv',       '',                     '',           0,  0,  0 ],
  [ 'Nobs',       'num_obs',              'numeric',    9,  3,  3 ],
  [ 'Nstruct',    '',                     '',           0,  0,  0 ],
  [ 'Nfilt',      '',                     '',           0,  0,  0 ],
  [ 'PDB_id',     '',                     '',           0,  0,  0 ],
  [ 'PDB_pos',    '',                     '',           0,  0,  0 ],
  [ 'PDB_ch',     '',                     '',           0,  0,  0 ],
  [ 'ident',      'ali_ide',              'numeric',   21,  0,  0 ],
  [ 'length',     'ali_len',              'numeric',   22,  0,  0 ],
  [ 'NormASA',    'acc_normed',           'numeric',   23,  9,  0 ],
  [ 'SecStr',     'sec_str',              '',          24,  0,  0 ],
  [ 'MapReg',     'map_region',           '',          25,  0,  0 ],
  [ 'dVol',       'delta_volume',         'numeric',   10,  4,  0 ],
  [ 'dProp',      'delta_prop',           'numeric',   26,  0,  0 ],
  [ 'B-fact',     'b_fact',               'numeric',   27, 10,  0 ],
  [ 'H-bonds',    'h_bonds',              'numeric',    0,  0,  0 ],
  [ 'AveNHet',    'het_cont_ave_num',     'numeric',   28,  0,  0 ],
  [ 'MinDHet',    'het_cont_min_dist',    'numeric',   29,  0,  0 ],
  [ 'AveNInt',    'inter_cont_ave_num',   'numeric',   30,  0,  0 ],
  [ 'MinDInt',    'inter_cont_min_dist',  'numeric',   31,  0,  0 ],
  [ 'AveNSit',    '',                     '',           0,  0,  0 ],
  [ 'MinDSit',    '',                     '',           0,  0,  0 ],
  [ 'Transv',     'transversion',         'yesno',     11,  0,  0 ],
  [ 'CodPos',     '',                     '',           0,  0,  0 ],
  [ 'CpG',        'cpg',                  '',          12,  0,  0 ],
  [ 'MinDJnc',    '',                     '',           0,  0,  0 ],
  [ 'PfamHit',    'pfam_hit',             'anyyes',    13,  5,  4 ],
  [ 'IdPmax',     'id_p_max',             'misszero',  14,  6,  5 ],
  [ 'IdPSNP',     'id_p_snp',             'misszero',  15,  0,  0 ],
  [ 'IdQmin',     'id_q_min',             'misszero',  16,  7,  6 ],
  # Mappings for generated features (skip if PPH2 name is empty)
  # cpg_var1_var2() takes care of both cpg_var1_var2 and cpg_transition
  [ \&cpg_var1_var2,    'cpg_var1_var2',    '',        17,  0,  0 ],
  [ '',                 'cpg_transition',   '',        18,  8,  7 ],
  [ \&charge_change,    'charge_change',    'numeric', 19,  0,  0 ],
  [ \&hydroph_change,   'hydroph_change',   'numeric', 20,  0,  0 ],
  # delta_volume_new takes care of both delta_volume_new and delta_prop_new
  [ \&delta_volume_new, 'delta_volume_new', 'numeric', 32,  0,  0 ],
  [ '',                 'delta_prop_new',   'numeric', 33, 11,  8 ]
);

# Lists of values for all nominal attributes
my %ATTRIBUTE = (
  'nt1'            => [ qw{A C G T} ],
  'nt2'            => [ qw{A C G T} ],
  'site'           => [ qw{NO YES} ],
  'region'         => [ qw{NO TRANSMEM INTRAMEM COMPBIAS REPEAT COILED SIGNAL PROPEP} ],
  'sec_str'        => [ qw{HELIX SHEET OTHER} ],
  'map_region'     => [ qw{ALPHA BETA OTHER} ],
  'transversion'   => [ qw{NO YES} ],
  'cpg'            => [ qw{NO YES} ],
  'pfam_hit'       => [ qw{NO YES} ],
  'cpg_var1_var2'  => [ qw{NO YES_S_F YES_S_Y YES_S_C YES_F_L YES_C_W YES_S_L YES_L_F YES_S_W YES_W_C YES_L_I YES_L_V YES_P_S YES_P_L YES_P_H YES_P_R YES_P_T YES_P_A YES_H_Y YES_H_N YES_H_D YES_R_C YES_R_L YES_R_P YES_R_H YES_R_S YES_R_G YES_H_Q YES_P_Q YES_Q_K YES_Q_E YES_R_Q YES_L_M YES_Q_H YES_R_W YES_T_I YES_T_N YES_T_S YES_I_M YES_N_K YES_S_R YES_T_K YES_T_R YES_M_I YES_T_M YES_K_N YES_V_F YES_V_L YES_V_I YES_A_S YES_A_P YES_A_T YES_A_V YES_A_D YES_A_G YES_D_Y YES_D_H YES_D_N YES_G_C YES_G_R YES_G_S YES_D_E YES_A_E YES_E_Q YES_E_K YES_V_M YES_E_D YES_G_W} ],
  'cpg_transition' => [ qw{NO YES_TRANSITION YES_TRANSVERSION} ],
  'class'          => [ qw{NEUTRAL DELETERIOUS} ]
);

my %VOLUME = (
  'A' => 88, 'C' => 108, 'D' => 111, 'E' => 138, 'F' => 190,
  'G' => 60, 'H' => 153, 'I' => 167, 'K' => 168, 'L' => 167,
  'M' => 163,'N' => 114, 'P' => 112, 'Q' => 144, 'R' => 173,
  'S' => 89, 'T' => 116, 'V' => 140, 'W' => 227, 'Y' => 193
);

my %PROPENSITY = (
  'A' => [  0.59,   0.06,  -0.03,  -0.20,  -0.23,  -0.18,   0.18 ],
  'R' => [ -2.99,  -1.14,   0.02,   0.33,   0.47,   0.25,  -0.31 ],
  'N' => [ -1.48,  -0.55,  -0.27,   0.06,   0.26,   0.45,   0.39 ],
  'D' => [ -1.35,  -0.69,  -0.33,   0.10,   0.28,   0.42,   0.50 ],
  'C' => [  0.61,   0.52,   0.51,  -0.02,  -0.75,  -1.80,  -2.00 ],
  'Q' => [ -1.55,  -0.58,  -0.25,   0.06,   0.40,   0.44,  -0.17 ],
  'E' => [ -2.13,  -0.98,  -0.35,  -0.01,   0.40,   0.58,   0.28 ],
  'G' => [  0.31,  -0.23,  -0.07,  -0.21,  -0.12,   0.05,   0.68 ],
  'H' => [ -1.09,  -0.11,   0.25,   0.38,   0.10,  -0.18,  -0.54 ],
  'I' => [  0.87,   0.54,   0.09,  -0.04,  -0.79,  -1.29,  -2.00 ],
  'L' => [  0.59,   0.61,   0.28,  -0.11,  -0.76,  -1.21,  -1.42 ],
  'K' => [ -4.06,  -2.28,  -1.22,  -0.08,   0.60,   0.75,   0.39 ],
  'M' => [  0.53,   0.61,  -0.01,  -0.11,  -0.40,  -1.17,  -0.62 ],
  'F' => [  0.47,   0.64,   0.36,  -0.01,  -0.86,  -1.31,  -1.75 ],
  'P' => [ -0.98,  -0.46,  -0.28,  -0.03,   0.11,   0.41,   0.74 ],
  'S' => [ -0.29,  -0.21,   0.04,  -0.11,  -0.01,   0.21,   0.49 ],
  'T' => [ -0.49,  -0.16,  -0.11,   0.10,   0.23,   0.09,  -0.15 ],
  'W' => [ -0.29,   0.53,   0.63,   0.19,  -0.52,  -1.42,  -2.54 ],
  'Y' => [ -0.76,   0.27,   0.53,   0.43,  -0.22,  -0.95,  -1.87 ],
  'V' => [  0.69,   0.50,   0.19,  -0.13,  -0.46,  -1.13,  -1.76 ]
);

my %HYDROPHOBICITY = (
  'A', 0.02, 'R', -0.42, 'N', -0.77, 'D', -1.04,  'C',  0.77, 'Q', -1.10,  'E', -1.14, 'G', -0.80,  'H', 0.26,  'I', 1.81,
  'L', 1.14, 'K', -0.41, 'M',  1.00, 'F',  1.35,  'P', -0.09, 'S', -0.97,  'T', -0.77, 'W',  1.71,  'Y', 1.11,  'V', 1.13,
);

my $ASA_THRESH = 0.05;
my $CONTTHRESH = 6;     # Angstroms

my %opt;
die "Illegal option\n" unless getopts('psuvc:f:', \%opt);
my $sumfile = shift;
die "PPH summary file not specified\n" unless defined $sumfile;
die "PPH summary file does not exist: $sumfile\n" unless -e $sumfile;
my $setidx = 3; # pointer to default feature set, full 33 features
if (defined $opt{'f'}) {
  if    ($opt{'f'} == 33) { $setidx = 3; }
  elsif ($opt{'f'} == 11) { $setidx = 4; }
  elsif ($opt{'f'} ==  8) { $setidx = 5; }
  else  { die("Unsupported feature set: ".$opt{'f'}."\n"); }
} else {
  $opt{'f'} = '33';
}
my $signed = 1; # defaults to signed deltas
$signed = 0 if $opt{'u'};
my $parse_comments = $opt{'p'} || 0;  # parse comments to extract sids and class
my $sids = $opt{'s'} || 0; # input format with set IDs at array pos -2
my $forced_class = $opt{'c'} || ''; # user-supplied override for the class value
my $verbose = $opt{'v'} || 0;

# Compile a hash of feature indices by feature names
my %FTIDX;
for (my $i=0; $i<@PPH2WEKA; $i++) {
  # generated features use WEKA names for keys
  my $pname = ref($PPH2WEKA[$i]->[0]) eq 'CODE' ? $PPH2WEKA[$i]->[1] : $PPH2WEKA[$i]->[0];
  $FTIDX{$pname} = $i if defined $pname && length $pname;
}

# Compile indexing array for mapping all input fields
# to a (variable) set of output features
my @subset;
foreach my $rec (@PPH2WEKA) {
  my $idx = $rec->[$setidx];
  push @subset, ($idx ? $idx - 1 : undef);
}

open(FH, $sumfile) or die "Can't open summary file: $sumfile\n";
# WEKA features, AoA
my @weka;
# Unique values for each WEKA feature, HoH
my $input_line_no = 0;
my %values;
while (<FH>) {
  $input_line_no++;
  # Skip comments
  next if /^#/;
  # Skip lines with errors & warnings and empty lines
  next if /\tError:\s/;
  next if /^(ERROR|WARNING):\s/;
  next if /^\s*$/;
  chomp;
  # Trim comments if present
  my $comments = $1 if s/\s *#\s*(.+)$//;
  my @pf = split /\t/, $_, 51;  # preserve any trailing empty fields
  # 50 mandatory features
  die('Incorrect number of fields (' . scalar(@pf) . ") at line $input_line_no in summary file: $sumfile\n") unless @pf == 50;
  map { s/^\s+//; s/\s+$//; } @pf;
  # Skip lines with missing alignment data, i.e., the ones with 'unknown' prediction outcome
  next if $pf[$FTIDX{'prediction'}] eq 'unknown' && $pf[$FTIDX{'dScore'}] eq '';
  my ($set_id, $class);
  # Parse comments to extract sids and/or class columns
  if ($parse_comments && defined $comments) {
    my @c = split /\t/, $comments;
    map { s/^\s+//; s/\s+$//; } @c;
    $set_id = shift @c if @c > 1;
    $class  = shift @c;
  }
  my ($aa1, $aa2) = ($pf[$FTIDX{'aa1'}], $pf[$FTIDX{'aa2'}]);
  # Sanity check
  die "AA variant(s) missing or incorrect at line $input_line_no: $_\n" unless length($aa1)==1 && length($aa2)==1;

  # NormASA indicates the presence of 3D structure
  my $with3D = length $pf[$FTIDX{'NormASA'}] ? 1 : 0;
  my @wf; # array of all possible classifier input values (some are dummy)
  # Loop over all PPHv2 output fields as defined in the mappings array
  for ($i=0; $i<@PPH2WEKA; $i++) {

    my ($pname, $wname, $type) = ($PPH2WEKA[$i]->[0], $PPH2WEKA[$i]->[1], $PPH2WEKA[$i]->[2]);
    # skip completely if PPHv2 name is empty (feature processed and inserted elsewhere)
    next unless length $pname;
    # insert dummy values for fields not included in the output
    push(@wf, undef), next unless length $wname;
    # process and insert generated feature
    &$pname($i, \@wf, $aa1, $aa2), next if ref($pname) eq 'CODE';

    # Sanity check
    die "Summary file parsing error at line $input_line_no: $_\n" unless length($pname) && length($wname);
    my $pf = $pf[$i];
    if (defined $pf && length $pf) {

      # Specific features:
      # note: will still be set to '?' if undefined/missing
      if      ($pname eq 'SecStr') {
        if    ($pf eq 'H' || $pf eq 'G')  { $pf = 'HELIX'; }
        elsif ($pf eq 'E')                { $pf = 'SHEET'; }
        else                              { $pf = 'OTHER'; }
      } elsif ($pname eq 'MapReg') {
        if    ($pf eq 'A' || $pf eq 'a')  { $pf = 'ALPHA'; }
        elsif ($pf eq 'B' || $pf eq 'b')  { $pf = 'BETA';  }
        else                              { $pf = 'OTHER'; }
      } elsif ($pname eq 'dVol') {
        # Need to recalculate since PPH dVol is abs() value
        die('Side chain volume undefined for residue: '.$aa1."\n") unless exists $VOLUME{$aa1};
        die('Side chain volume undefined for residue: '.$aa2."\n") unless exists $VOLUME{$aa2};
        $pf = $VOLUME{$aa1} - $VOLUME{$aa2};
        $pf = abs($pf) unless $signed;
      } elsif ($pname eq 'dProp') {
        if ($with3D) {
          my $j = &propens_interval($pf[$FTIDX{'NormASA'}]);
          my $dp = $PROPENSITY{$aa1}->[$j] - $PROPENSITY{$aa2}->[$j];
          $pf = $signed ? $dp : abs($dp);
        } else {
          $pf = '?';
        }
      } elsif ($pname eq 'dScore') {
        # Recalculate if we are in signed mode
        $pf = sprintf('%.3f', $pf[$FTIDX{'Score1'}] - $pf[$FTIDX{'Score2'}]) if $signed;
      } elsif ($pname eq 'CpG') {
        my $cpgmod = $pf % 2;
        $pf = $cpgmod ? 'YES' : 'NO';

      # Generic feature types:
      # numeric: remove extra + sign
      } elsif   ($type eq 'numeric') {
        $pf =~ s/^\+//;
      # nominal: translate any value other than NO to YES
      } elsif ($type eq 'anyyes') {
        $pf = 'YES' unless $pf eq 'NO';
      # binary: translate boolean {0,1} to {YES,NO}
      } elsif ($type eq 'yesno') {
        $pf = $pf ? 'YES' : 'NO';
      }

    # feature is undefined / missing
    } else {

      # Specific features:
      if      ($pname eq 'AveNHet' || $pname eq 'AveNInt') {
        $pf = $with3D ? 0 : '?';
      } elsif ($pname eq 'MinDHet' || $pname eq 'MinDInt') {
        # CONTTHRESH*2 Angstroms denotes arbitrary large distance
        $pf = $with3D ? $CONTTHRESH*2 : '?';

      # Generic feature types:
      # Force missno type to NO if undefined/missing
      } elsif ($type eq 'missno') {
        $pf = 'NO';
      # Force misszero types to 0 if undefined/missing
      } elsif   ($type eq 'misszero') {
        $pf = 0;
      # All other undefined / missing features are denoted by '?'
      } else { $pf = '?'; }

    }

    &values($wname, $pf, $type) if $verbose;
    push @wf, $pf;

  }

  my @wfset; $#wfset = scalar(grep(defined, @subset)) - 1;
  for ($i=0; $i<@wf; $i++) {
    $wfset[$subset[$i]] = $wf[$i] if defined $subset[$i];
  }

  push @wfset, length $set_id ? $set_id : '?' if $sids;
  push @wfset, $forced_class || $class || 'NEUTRAL';

  push @weka, \@wfset;

}
close(FH);

$opt{'f'} += 1 if $sids; # set_id is an extra feature
my $basename = fileparse($sumfile, '.output');
my $testfile = $basename . '.f' . $opt{'f'} . '.arff';
open(FH, ">$testfile") or die "Can't create temporary file: $testfile\n";
# Create an empty .arff if no valid lines were found in the PPH file:
# this can happen especially in a cluster (-q) mode when each PPH file
# holds output of a single query (which may easily fail). Empty arff
# serves as a prediction pipleine failure flag for further processing.
close(FH), exit unless @weka;

my @idxout; $#idxout = scalar(grep(defined, @subset)) - 1;
for ($i=0; $i<@subset; $i++) {
  $idxout[$subset[$i]] = $i if defined $subset[$i];
}

my $relation = '@relation PPHv2.f' . $opt{'f'};
print FH $relation, "\n\n";
foreach my $idx (@idxout) {
  warn("WARNING: pph2arff: Undefined output feature index skipped\n"), next unless defined $idx;
  my $rec = $PPH2WEKA[$idx];
  # Numeric attribute
  if ($rec->[2] eq 'numeric' || $rec->[2] eq 'misszero') {
    print FH '@attribute ', $rec->[1], " numeric\n";
  # Nominal attribute
  } else {
    print FH '@attribute ', $rec->[1], ' {', join(',', @{ $ATTRIBUTE{$rec->[1]} }), "}\n";
  }
}
print FH '@attribute class {', join(',', @{ $ATTRIBUTE{'class'} }), "}\n\n";
print FH '@data', "\n";

foreach my $rec (@weka) {
  print FH join(',', @$rec), "\n";
}
close(FH);

if ($verbose) {
  my @wnames; $#wnames = scalar(grep(defined, @subset)) - 1;
  for ($i=0; $i<@PPH2WEKA; $i++) {
    $wnames[ $subset[$i] ] = $PPH2WEKA[$i]->[1] if defined $subset[$i];
  }
  foreach $wname (@wnames) {
    printf '%22s', "$wname: ";
    my $type = ref $values{$wname};
    if      ($type eq 'HASH') {
      my $total = 0;
      my $missing = $values{$wname}->{'?'} || 0;
      $total += $missing;
      foreach $name (sort keys %{ $values{$wname} }) {
        next if $name eq '?';
        my $count = $values{$wname}->{$name};
        $total += $count;
        print "$name($count) + ";
      }
      print "?($missing) = $total\n";
    } elsif ($type eq 'ARRAY') {
      my ($min, $max, $sum, $count, $missing) = @{ $values{$wname} };
      print
        sprintf('%.4g:%.4g, AVE=', $min, $max),
        sprintf('%.3g', $sum / $count),
        ', COUNT=', $count,
        ', ?=', $missing,
        ', TOTAL=', ($missing + $count),
        "\n";
    } else {
      die "Illegal type encountered in values: $type\n";
    }
  }
}


# Subroutines
#------------------------------------------------------------------------------
sub values {
  my ($wname, $value, $type) = @_;
  # Count unique values for nominal features and check ranges for type numeric
  if ($type eq 'numeric' || $type eq 'misszero') {
    if ($value eq '?') {
      $values{$wname}->[4]++;
    } else {
      my ($min, $max, $sum, $count, $missing) = @{ $values{$wname} };
      if (defined $min) {
        if    ($value < $min) { $min = $value; }
        elsif ($value > $max) { $max = $value; }
        $count++;
        $sum += $value;
      } else {
        $min = $max = $sum = $value;
        $count = 1;
        $missing = 0 unless defined $missing;
      }
      @{ $values{$wname} } = ( $min, $max, $sum, $count, $missing );
    }
  } else {
    $values{$wname}->{$value}++;
  }
}
#------------------------------------------------------------------------------
sub charge {
  my($aa) = shift;

  if ($aa eq 'K' || $aa eq 'R') {
    return 1;
  } elsif ($aa eq 'E' || $aa eq 'D') {
    return -1;
  }

  return 0;
}
#------------------------------------------------------------------------------
sub propens_interval {
  my($x) = shift;
  ## 0 (0,5) [5,15) [15,30) [30,50) [50,75) [75,...)
  if($x==0) { return 0; }
  elsif (  0<$x && $x<5 ) { return 1; }
  elsif ( 5<=$x && $x<15) { return 2; }
  elsif (15<=$x && $x<30) { return 3; }
  elsif (30<=$x && $x<50) { return 4; }
  elsif (50<=$x && $x<75) { return 5; }
  else { return 6; }
}
#------------------------------------------------------------------------------
sub cpg_var1_var2 {
  my ($i, $wf, $aa1, $aa2) = @_;
  my $cpg          = $wf->[$FTIDX{'CpG'}];
  my $transversion = $wf->[$FTIDX{'Transv'}];
  die "cpg_var1_var2: invalid input values\n" unless defined $cpg && defined $transversion;
  my ($cpg_var1_var2, $cpg_transition);
  if ($cpg ne '?' && $transversion ne '?') {
    $cpg_var1_var2 = $cpg . '_' . $aa1 . '_' . $aa2;
    $cpg_var1_var2 = 'NO' unless $cpg eq 'YES' && $transversion eq 'NO';
    $cpg_transition =
      $cpg eq 'YES' ?
        ('YES_' . ($transversion eq 'YES' ? 'TRANSVERSION' : 'TRANSITION')) : 'NO';
  } else {
    $cpg_var1_var2 = $cpg_transition = '?';
  }
  # This relies upon the two features being consecutive
  push @{ $wf }, $cpg_var1_var2;
  push @{ $wf }, $cpg_transition;
  if ($verbose) {
    my ($wname, $type) = ($PPH2WEKA[$i]->[1], $PPH2WEKA[$i]->[2]);
    &values($wname, $cpg_var1_var2, $type);
    ($wname, $type) = ($PPH2WEKA[$i+1]->[1], $PPH2WEKA[$i+1]->[2]);
    &values($wname, $cpg_transition, $type);
  }
}
#------------------------------------------------------------------------------
sub charge_change {
  my ($i, $wf, $aa1, $aa2) = @_;
  my $cc;

  if ($aa1 eq 'H' && charge($aa2) == 1 || $aa2 eq 'H' && charge($aa1) == 1) {
    $cc = 0;
  } else {
    $cc = abs(charge($aa1) - charge($aa2));
  }
  push @{ $wf }, $cc;
  if ($verbose) {
    my ($wname, $type) = ($PPH2WEKA[$i]->[1], $PPH2WEKA[$i]->[2]);
    &values($wname, $cc, $type);
  }
}
#------------------------------------------------------------------------------
sub hydroph_change {
  my ($i, $wf, $aa1, $aa2) = @_;
  my $hyd = abs($HYDROPHOBICITY{$aa1} - $HYDROPHOBICITY{$aa2});
  push @{ $wf }, $hyd;
  if ($verbose) {
    my ($wname, $type) = ($PPH2WEKA[$i]->[1], $PPH2WEKA[$i]->[2]);
    &values($wname, $hyd, $type);
  }
}
#------------------------------------------------------------------------------
sub delta_volume_new {
  my ($i, $wf, $aa1, $aa2) = @_;
  my $asa_normed   = $wf->[$FTIDX{'NormASA'}];
  my $delta_volume = $wf->[$FTIDX{'dVol'}];
  my $delta_prop   = $wf->[$FTIDX{'dProp'}];
  my ($delta_volume_new, $delta_prop_new);
  if($asa_normed ne '?') {
    $delta_volume_new = $asa_normed > $ASA_THRESH ? 0 : $delta_volume;
    $delta_prop_new   = $asa_normed > $ASA_THRESH ? 0 : $delta_prop;
  } else {
    $delta_volume_new = $delta_prop_new = '?';
  }
  # This relies upon the two features being consecutive
  push @{ $wf }, $delta_volume_new;
  push @{ $wf }, $delta_prop_new;
  if ($verbose) {
    my ($wname, $type) = ($PPH2WEKA[$i]->[1], $PPH2WEKA[$i]->[2]);
    &values($wname, $delta_volume_new, $type);
    ($wname, $type) = ($PPH2WEKA[$i+1]->[1], $PPH2WEKA[$i+1]->[2]);
    &values($wname, $delta_prop_new, $type);
  }
}
#------------------------------------------------------------------------------
