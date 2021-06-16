package PPH::Data;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Data

=head1 DESCRIPTION

Contains biological data and other data objects for PolyPhen

=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2011-09-30 16:17:16 -0400 (Fri, 30 Sep 2011) $
 $LastChangedRevision: 370 $
 $LastChangedBy: ivan $

=head1 ROUTINES / VARIABLES

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use Carp qw(confess);
use strict;
use warnings;

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION = 0.1;
    @ISA           = qw(Exporter);
    @EXPORT        = qw();
    @EXPORT_OK     = qw(
      %Effect
      %BONDS %SITES1 %SITES2
      @PhiPsiMap %MaxAccessibility %AA_Volume
      @Non_biol_ligands @Predictions
      %Propensities @GARBAGE_FILES
      %BLOSUM62N
    );
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

#----------------------------------------

=head2 @Predictions

 Function: Possible predictions set by PPH::Rules

=cut

#----------------------------------------
our @Predictions = (
    'unknown',
    'benign',
    'possibly damaging',
    'probably damaging',
);
#----------------------------------------

#----------------------------------------------------------------------

=head2 %Effect

 Function: Predicted effect of a substituion set by PPH::Rules

=cut

#----------------------------------------

our %Effect = (
    '1' =>   {
        Name  => 'structural effect',
    },

    '1.1' => {
        Name  => 'buried site',
    },

    '1.1.1'  => {
        Name  => 'hydrophobicity',
        Comment  => 'Hydrophobicity change at buried site',
    },

    '1.1.2'  => {
        Name  => 'overpacking',
        Comment  => 'Overpacking at buried site',
    },

    '1.1.3'  => {
        Name  => 'cavity',
        Comment  => 'Cavity creation at buried site',
    },

    '1.2' => {
        Name  => 'bond formation',
        Comment  => 'Disruption of annotated bond formation site',
    },

    '1.2.1'  => {
        Name  => 'covalent bond',
    },

    '1.2.1.1'   => {
        Name  => 'disulphide',
    },

    '1.2.1.2'   => {
        Name  => 'thioesther',
    },

    '1.2.1.3'   => {
        Name  => 'thioether',
    },

    '1.2.2'  => {
        Name  => 'non-covalent',
    },

    '1.2.2.1'   => {
        Name  => 'hydrogen',
    },

    '1.2.2.2'   => {
        Name  => 'salt bridge (electrostatic)',
    },



    '2'   => {
        Name  => 'functional effect',
    },

    '2.1' => {
        Name  => 'indirect',
        Comment  => 'Contact with functional site',
    },

    '2.2' => {
        Name  => 'functional site',
        Comment  => 'Disruption of annotated functional site',
    },

    '2.2.1'  => {
        Name  => 'signal peptide',
    },

    '2.2.2'  => {
        Name  => 'transmembrane',
        Comment  => 'Improper substitution in the transmembrane region',
    },

    '2.2.3'  => {
        Name  => 'ligand binding',
        Comment  => 'Disruption of ligand binding site',
    },

    '2.2.4'  => {
        Name  => 'protein interaction',
    },


    '3'   => {
        Name  => 'unknown',
        Comment  => 'Derived from multiple alignment',

    },


);
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 UniProt Annotations

 Function: The given annotations extracted from UniProt will be used to
           calculate e.g. transmembrane PHAT scores or distances to
           sites such as catalytic site or ligand

=head3 BONDS

=cut

our %BONDS;
@BONDS{ qw(DISULFID CROSSLNK) } = ();


=head3 %REGIONS

=cut

our %REGIONS;
@REGIONS{ qw(TRANSMEM INTRAMEM COMPBIAS REPEAT COILED SIGNAL PROPEP) } = ();

=head3 %SITES1

 Function: Critical sites - c.f. PPH::Seq::get_uniprot_annotation
           The sites are not in %REGIONS

=cut

our %SITES1;
@SITES1{ qw(BINDING ACT_SITE LIPID METAL) } = ();

=head3 %SITES2

=cut

our %SITES2;
@SITES2{ qw(SITE DISULFID CROSSLNK MOD_RES CARBOHYD NON_STD) } = ();

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 Constants for PPH::3D

=head3 @PhiPsiMap

 Function: needed for mapping phipsi angles (c.f. PPH::3D)

=cut

#----------------------------------------
our @PhiPsiMap = qw(
  bBBBBBBBBBBbbbbb.....pppppgggggbbbbb
  bBBBBBBBBBBBbbbbb......ggggggggbbbbb
  bBBBBBBBBBBBbbbbbbb....ggggggggbbbbb
  bbBBBBBBBBBBBbbbbbb....ggggggggbbbbb
  bbBBBBBBBBBBBBbbbbb....ggggggggbbbbb
  bbBBBBBBBBBBBbbbbbb....gggggggggbbbb
  bbbBBBBBBBBBbbbbbbb...llllllgggggbbb
  bbbbBBBBBBBbbbbbbbblllllllllgggggbbb
  bbbbbBbbBbbbbbbbbbblllllllllgggggbbb
  bbbbbbbbbbbbbbbbbbblllllllllgggggbbb
  bbbbbbbbbbbbbbbbbbblllllllllgggggbbb
  bbbbbbbbbbbbbbb....llllllllllggggbbb
  bbbbbbbbbbbbbbb....llllllllllggggbbb
  aaaaaaaaaaaaaaa....llllLlllllggggbbb
  aaaaaaaaaaaaaaa....llllLLllllggggbbb
  aaaaaaaAaaaaaaaa...llllllllllggggbbb
  aaaaaaAAAAaaaaaaa...lllllllllggggbbb
  aaaaaaAAAAAaaaaaa...lllllllllggggbbb
  aaaaaAAAAAAAaaaaaa..lllllllllggggbbb
  aaaaaaAAAAAAAaaaaaa.lllllllllggggbbb
  aaaaaaAAAAAAAAaaaaa.lllllllllggggbbb
  aaaaaaaAAAAAAAaaaaaa....gggggggggbbb
  aaaaaaaaAAAAAAAaaaaa...ggggggggggbbb
  aaaaaaaaaAAAAAAaaaaaa..ggggggggggbbb
  aaaaaaaaaaaAAAAaaaaaa..ggggggggggbbb
  aaaaaaaaaaaaaaaaaaaaa..ggggggggggbbb
  aaaaaaaaaaaaaaaaaaaaa.gggggggggggbbb
  .aaaaaaaaaaaaaaaaaaaa.gggggggggggggg
  gggaaaaaaaaaaaaaaaaa..gggggggggggggg
  gggbbbbbbbbbbbbbb.....gggggggggggggg
  bbbbbbbbbbbbbb.......pppppgggggggggg
  bbbbbbbbbbbbbb.......pppppgggggggggg
  bbbbbbbbbbbbbbb......ppPpppggggggggg
  bbbbbbbbbbbbbbb......ppPpppggggggbbb
  bbbbbbbbbbbbbbb......ppPPppggggggbbb
  bbbbbbbbbbbbbbb......ppppppggggggbbb
);

#----------------------------------------

=head3 %MaxAccessibility

 Function: Maximum accessibility of residues

=cut

#----------------------------------------
our %MaxAccessibility = (
  'A' => 103, 'C' => 105, 'D' => 164, 'E' => 187,
  'F' => 169, 'G' =>  87, 'H' => 170, 'I' => 133,
  'K' => 197, 'L' => 137, 'M' => 157, 'N' => 160,
  'P' => 133, 'Q' => 181, 'R' => 234, 'S' => 125,
  'T' => 139, 'V' => 120, 'W' => 195, 'Y' => 181
);
#----------------------------------------

#----------------------------------------

=head3 %AA_Volume

 Function: Volumina of AA

=cut

#----------------------------------------
our %AA_Volume = (
  'A' => 88, 'C' => 108, 'D' => 111, 'E' => 138, 'F' => 190,
  'G' => 60, 'H' => 153, 'I' => 167, 'K' => 168, 'L' => 167,
  'M' => 163,'N' => 114, 'P' => 112, 'Q' => 144, 'R' => 173,
  'S' => 89, 'T' => 116, 'V' => 140, 'W' => 227, 'Y' => 193
);
#----------------------------------------

#----------------------------------------

=head3 @Non_biol_ligands

 Function: Ligands which are not considered for distance measurments

=cut

#----------------------------------------
our @Non_biol_ligands = qw( PO4 PI SO4 SUL CL BR NO3 SCN NH4 K NA LI MG DOD NAG MAN GOL SO4 CL CO3 FS4 );
#----------------------------------------

#----------------------------------------

=head3 %Propensities

=cut

our %Propensities = (
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
#----------------------------------------

#------------------------------------------------------------------------

#----------------------------------------

=head2 blastmatrix

 Usage   : PPH::Data->blastmatrix($filename)
 Function: Simply Cache BLAST matrix
 Args    : filename of BLAST matrix
 Returns : hash of matrix

=cut

#----------------------------------------
sub blastmatrix {
    my $file = shift;
    my (@alph, $i, %matr);

    open my $MAT, '<', $file or confess "Can't read matrix from $file";

    while (<$MAT>) {
        (/^\#/) and next;

        if ( /^\s+(\w.*)$/ ) {
            @alph = (split);    ### leading space is lost only for $_
        } elsif ( /^[A-Z\*] / ) {
            my($let, @row) = (split);
            foreach $i (0..($#row)) {
                $matr{$let}{$alph[$i]} = $row[$i];
            }
        }                       # end elsif
    }
    close $MAT;
    return %matr;
}
#----------------------------------------

# my @BLA = qw( A R N D C Q E G H I L K M F P S T W Y V );
#
# # BLOSUM62 raw frequencies = fij pairs (off-diagonals = 2*fij)
# # Source: http://blocks.fhcrc.org/blocks/uploads/blosum/blosum.tar.gz
# #         (file: blosum62.out)
# my @BL62F = qw(
#  26782.81
#   5848.09  22123.18
#   4857.14   4930.09  17616.72
#   5400.73   3953.21   9269.63  26504.98
#   3962.65    980.88   1092.50    994.42  14864.71
#   4794.70   6194.03   3813.41   4106.64    770.86   9131.13
#   7444.98   6711.13   5505.96  12248.73    955.29   8817.09  20100.50
#  14491.39   4291.00   7124.25   6284.36   1917.42   3409.41   4829.16  47098.64
#   2759.96   3091.55   3563.41   2376.71    572.50   2613.56   3405.63   2387.39  11561.82
#   7943.72   3098.65   2477.68   3076.61   2730.06   2220.18   3037.90   3450.79   1447.42  22975.42
#  11009.91   6028.42   3411.38   3787.90   3907.74   4030.19   4990.96   5198.87   2459.23  28361.76  46272.78
#   8338.88  15533.06   6080.29   6092.99   1248.97   7716.37  10296.37   6327.07   2958.56   3900.98   6138.15  20074.90
#   3341.90   2001.08   1319.20   1156.88    939.80   1843.69   1692.11   1825.96    953.44   6249.52  12282.50   2264.38   5042.88
#   4076.59   2321.82   1868.88   1894.29   1280.49   1351.88   2122.42   2984.24   2019.30   7589.45  13492.91   2364.01   2965.57  22771.30
#   5374.35   2386.63   2143.28   3083.07    899.80   2109.68   3542.26   3398.98   1190.29   2508.69   3524.79   3930.17   1017.22   1308.84  23753.47
#  15579.38   5646.43   7840.40   6985.53   2599.50   4717.08   7360.29   9553.90   2753.84   4291.98   6049.14   7728.29   2132.99   2975.00   4151.82  15680.23
#   9264.29   4435.93   5571.65   4724.65   2318.32   3437.77   5106.27   5446.48   1853.12   6716.05   8281.94   5847.29   2515.49   2896.40   3366.54  11712.22  15591.62
#   1003.70    662.10    402.65    404.11    360.67    566.56    660.02   1015.14    377.94    901.59   1824.09    677.71    495.17   2115.89    352.63    715.91    712.02   8060.58
#   3239.21   2308.18   1745.40   1491.09    862.22   1684.00   2168.88   2079.76   3790.82   3443.82   5505.85   2489.43   1423.89  10562.84   1126.86   2566.42   2346.22   2211.24  12765.13
#  12627.73   3939.43   2993.54   3278.68   3390.37   2905.57   4232.85   4539.45   1616.69  29832.35  23617.94   4824.05   5761.66   6419.39   3102.50   5877.26   9070.39    886.51   3859.61  24458.47
# );
#
# my $BL62S;    # Total sum of all AA pairs (e.g.: 1245852.53 for BLOSUM62)
# my @BL62idx;  # Array of indices into plain BL62F array which step through it by AA column
# foreach my $i (0 .. $#BLA) {
#   my $step;
#   $step += $_ for (0 .. $i);
#   push @BL62idx, $step;
# }
#
# # Load fij frequencies into a hash representing full 20x20 symmetrical matrix
# our %BLOSUM62N;
# for my $i (0 .. $#BLA) {
#   for my $j ($i .. $#BLA) {
#     my $k = $BL62idx[$j] + $i;
#     my $fij = $BL62F[$k];
#     if ($i == $j) {
#       $BLOSUM62N{$BLA[$i]}{$BLA[$i]} = $fij;
#     } else {
#       $BLOSUM62N{$BLA[$i]}{$BLA[$j]} = $BLOSUM62N{$BLA[$j]}{$BLA[$i]} = $fij / 2;
#     }
#     $BL62S += $fij;
#   }
# }
#
# # BLOSUM62N marginal frequencies, i.e. frequencies of each of the 20 AA residues
# my %BL62AA;
# foreach my $aa1 (@BLA) {
#   foreach my $aa2 (@BLA) {
#     if ($aa1 eq $aa2) {
#       $BL62AA{$aa1} += $BLOSUM62N{$aa1}{$aa1};
#     } else {
#       $BL62AA{$aa1} += ($BLOSUM62N{$aa1}{$aa2} + $BLOSUM62N{$aa2}{$aa1}) / 2;
#     }
#   }
# }
#
# # Normalize matrix by dividing each column by its marginal frequency
# foreach my $aa1 (@BLA) {
#   foreach my $aa2 (@BLA) {
#     $BLOSUM62N{$aa1}{$aa2} /= $BL62AA{$aa1};
#   }
# }
#
# # my %BLOSUM62N; # Symmetrical 20x20 matrix of AA pairs probabilities (qij)
# # my $BL62Sx2 = $BL62S * 2;
# # for my $i (0 .. $#BLA) {
# #   for my $j ($i .. $#BLA) {
# #     my $k = $BL62idx[$j] + $i;
# #     my $fij = $BL62F[$k];
# #     my $qij;
# #     if ($i == $j) {
# #       $qij = $fij / $BL62S;
# #     } else {
# #       $qij = $fij / $BL62Sx2;
# #     }
# #     $BLOSUM62N{$BLA[$i]}{$BLA[$j]} = $BLOSUM62N{$BLA[$j]}{$BLA[$i]} = $qij;
# #   }
# # }
#
# print "our %BLOSUM62N = (\n";
# for my $i (0 .. $#BLA) {
#   print "  '$BLA[$i]' => { ";
#   my $s;
#   for my $j (0 .. $#BLA) {
#     $s .= sprintf "'$BLA[$j]'=>%.8f, ", $BLOSUM62N{$BLA[$i]}{$BLA[$j]};
#   }
#   $s =~ s/, $//;
#   print "$s },\n";
# }
# print ");\n";

# Full symmetrical 20x20 matrix of AA pairs frequencies normalized by column marginal AA frequencies
# For details, see: http://blocks.fhcrc.org/blocks/uploads/blosum/
our %BLOSUM62N = (
  'A' => { 'A'=>0.28966145, 'R'=>0.03162413, 'N'=>0.02626547, 'D'=>0.02920499, 'C'=>0.02142843, 'Q'=>0.02592782, 'E'=>0.04025947, 'G'=>0.07836364, 'H'=>0.01492476, 'I'=>0.04295646, 'L'=>0.05953719, 'K'=>0.04509333, 'M'=>0.01807166, 'F'=>0.02204457, 'P'=>0.02906234, 'S'=>0.08424706, 'T'=>0.05009757, 'W'=>0.00542761, 'Y'=>0.01751635, 'V'=>0.06828571 },
  'R' => { 'A'=>0.04547219, 'R'=>0.34404031, 'N'=>0.03833422, 'D'=>0.03073843, 'C'=>0.00762689, 'Q'=>0.04816206, 'E'=>0.05218281, 'G'=>0.03336494, 'H'=>0.02403854, 'I'=>0.02409374, 'L'=>0.04687435, 'K'=>0.12077827, 'M'=>0.01555952, 'F'=>0.01805345, 'P'=>0.01855739, 'S'=>0.04390417, 'T'=>0.03449185, 'W'=>0.00514820, 'Y'=>0.01794740, 'V'=>0.03063128 },
  'N' => { 'A'=>0.04366197, 'R'=>0.04431773, 'N'=>0.31672165, 'D'=>0.08332688, 'C'=>0.00982074, 'Q'=>0.03427964, 'E'=>0.04949436, 'G'=>0.06404155, 'H'=>0.03203233, 'I'=>0.02227245, 'L'=>0.03066569, 'K'=>0.05465715, 'M'=>0.01185860, 'F'=>0.01679980, 'P'=>0.01926645, 'S'=>0.07047919, 'T'=>0.05008487, 'W'=>0.00361952, 'Y'=>0.01568981, 'V'=>0.02690963 },
  'D' => { 'A'=>0.04041852, 'R'=>0.02958542, 'N'=>0.06937297, 'D'=>0.39672118, 'C'=>0.00744214, 'Q'=>0.03073368, 'E'=>0.09166826, 'G'=>0.04703152, 'H'=>0.01778706, 'I'=>0.02302504, 'L'=>0.02834826, 'K'=>0.04559932, 'M'=>0.00865797, 'F'=>0.01417667, 'P'=>0.02307338, 'S'=>0.05227900, 'T'=>0.03535880, 'W'=>0.00302432, 'Y'=>0.01115917, 'V'=>0.02453731 },
  'C' => { 'A'=>0.06441879, 'R'=>0.01594567, 'N'=>0.01776022, 'D'=>0.01616578, 'C'=>0.48329613, 'Q'=>0.01253148, 'E'=>0.01552967, 'G'=>0.03117053, 'H'=>0.00930684, 'I'=>0.04438120, 'L'=>0.06352615, 'K'=>0.02030387, 'M'=>0.01527785, 'F'=>0.02081628, 'P'=>0.01462759, 'S'=>0.04225876, 'T'=>0.03768775, 'W'=>0.00586323, 'Y'=>0.01401667, 'V'=>0.05511553 },
  'Q' => { 'A'=>0.05616709, 'R'=>0.07255942, 'N'=>0.04467186, 'D'=>0.04810687, 'C'=>0.00903017, 'Q'=>0.21393165, 'E'=>0.10328703, 'G'=>0.03993923, 'H'=>0.03061632, 'I'=>0.02600810, 'L'=>0.04721131, 'K'=>0.09039274, 'M'=>0.02159775, 'F'=>0.01583648, 'P'=>0.02471366, 'S'=>0.05525782, 'T'=>0.04027146, 'W'=>0.00663692, 'Y'=>0.01972707, 'V'=>0.03403705 },
  'E' => { 'A'=>0.05501381, 'R'=>0.04959111, 'N'=>0.04068565, 'D'=>0.09051055, 'C'=>0.00705900, 'Q'=>0.06515285, 'E'=>0.29706058, 'G'=>0.03568451, 'H'=>0.02516550, 'I'=>0.02244821, 'L'=>0.03688011, 'K'=>0.07608382, 'M'=>0.01250365, 'F'=>0.01568337, 'P'=>0.02617512, 'S'=>0.05438800, 'T'=>0.03773218, 'W'=>0.00487714, 'Y'=>0.01602668, 'V'=>0.03127815 },
  'G' => { 'A'=>0.07843686, 'R'=>0.02322569, 'N'=>0.03856109, 'D'=>0.03401506, 'C'=>0.01037833, 'Q'=>0.01845395, 'E'=>0.02613856, 'G'=>0.50985714, 'H'=>0.01292211, 'I'=>0.01867793, 'L'=>0.02813968, 'K'=>0.03424623, 'M'=>0.00988329, 'F'=>0.01615265, 'P'=>0.01839750, 'S'=>0.05171194, 'T'=>0.02947990, 'W'=>0.00549460, 'Y'=>0.01125702, 'V'=>0.02457047 },
  'H' => { 'A'=>0.04225614, 'R'=>0.04733293, 'N'=>0.05455730, 'D'=>0.03638843, 'C'=>0.00876521, 'Q'=>0.04001470, 'E'=>0.05214162, 'G'=>0.03655194, 'H'=>0.35403261, 'I'=>0.02216061, 'L'=>0.03765184, 'K'=>0.04529679, 'M'=>0.01459757, 'F'=>0.03091633, 'P'=>0.01822384, 'S'=>0.04216244, 'T'=>0.02837204, 'W'=>0.00578642, 'Y'=>0.05803904, 'V'=>0.02475220 },
  'I' => { 'A'=>0.04694037, 'R'=>0.01831028, 'N'=>0.01464090, 'D'=>0.01818005, 'C'=>0.01613224, 'Q'=>0.01311930, 'E'=>0.01795130, 'G'=>0.02039112, 'H'=>0.00855297, 'I'=>0.27152886, 'L'=>0.16759294, 'K'=>0.02305134, 'M'=>0.03692914, 'F'=>0.04484694, 'P'=>0.01482414, 'S'=>0.02536181, 'T'=>0.03968592, 'W'=>0.00532760, 'Y'=>0.02034993, 'V'=>0.17628283 },
  'L' => { 'A'=>0.04467415, 'R'=>0.02446110, 'N'=>0.01384212, 'D'=>0.01536990, 'C'=>0.01585617, 'Q'=>0.01635302, 'E'=>0.02025147, 'G'=>0.02109510, 'H'=>0.00997865, 'I'=>0.11508155, 'L'=>0.37551572, 'K'=>0.02490635, 'M'=>0.04983785, 'F'=>0.05474925, 'P'=>0.01430230, 'S'=>0.02454518, 'T'=>0.03360506, 'W'=>0.00740148, 'Y'=>0.02234071, 'V'=>0.09583288 },
  'K' => { 'A'=>0.05754650, 'R'=>0.10719344, 'N'=>0.04196000, 'D'=>0.04204764, 'C'=>0.00861913, 'Q'=>0.05325056, 'E'=>0.07105511, 'G'=>0.04366302, 'H'=>0.02041698, 'I'=>0.02692061, 'L'=>0.04235929, 'K'=>0.27707323, 'M'=>0.01562646, 'F'=>0.01631400, 'P'=>0.02712205, 'S'=>0.05333282, 'T'=>0.04035207, 'W'=>0.00467687, 'Y'=>0.01717952, 'V'=>0.03329070 },
  'M' => { 'A'=>0.05366944, 'R'=>0.03213646, 'N'=>0.02118577, 'D'=>0.01857898, 'C'=>0.01509277, 'Q'=>0.02960885, 'E'=>0.02717454, 'G'=>0.02932411, 'H'=>0.01531183, 'I'=>0.10036454, 'L'=>0.19725153, 'K'=>0.03636494, 'M'=>0.16197286, 'F'=>0.04762575, 'P'=>0.01633610, 'S'=>0.03425488, 'T'=>0.04039766, 'W'=>0.00795221, 'Y'=>0.02286705, 'V'=>0.09252972 },
  'F' => { 'A'=>0.03450269, 'R'=>0.01965099, 'N'=>0.01581748, 'D'=>0.01603254, 'C'=>0.01083758, 'Q'=>0.01144179, 'E'=>0.01796335, 'G'=>0.02525746, 'H'=>0.01709058, 'I'=>0.06423419, 'L'=>0.11419881, 'K'=>0.02000807, 'M'=>0.02509945, 'F'=>0.38545507, 'P'=>0.01107752, 'S'=>0.02517926, 'T'=>0.02451402, 'W'=>0.01790808, 'Y'=>0.08939982, 'V'=>0.05433125 },
  'P' => { 'A'=>0.05596804, 'R'=>0.02485417, 'N'=>0.02231994, 'D'=>0.03210684, 'C'=>0.00937044, 'Q'=>0.02197003, 'E'=>0.03688880, 'G'=>0.03539670, 'H'=>0.01239558, 'I'=>0.02612529, 'L'=>0.03670687, 'K'=>0.04092847, 'M'=>0.01059325, 'F'=>0.01363015, 'P'=>0.49473337, 'S'=>0.04323671, 'T'=>0.03505887, 'W'=>0.00367226, 'Y'=>0.01173503, 'V'=>0.03230918 },
  'S' => { 'A'=>0.10925397, 'R'=>0.03959688, 'N'=>0.05498260, 'D'=>0.04898763, 'C'=>0.01822959, 'Q'=>0.03307960, 'E'=>0.05161572, 'G'=>0.06699891, 'H'=>0.01931193, 'I'=>0.03009849, 'L'=>0.04242098, 'K'=>0.05419640, 'M'=>0.01495808, 'F'=>0.02086287, 'P'=>0.02911559, 'S'=>0.21992241, 'T'=>0.08213462, 'W'=>0.00502048, 'Y'=>0.01799761, 'V'=>0.04121563 },
  'T' => { 'A'=>0.07305861, 'R'=>0.03498194, 'N'=>0.04393828, 'D'=>0.03725880, 'C'=>0.01828238, 'Q'=>0.02711041, 'E'=>0.04026827, 'G'=>0.04295119, 'H'=>0.01461379, 'I'=>0.05296307, 'L'=>0.06531175, 'K'=>0.04611199, 'M'=>0.01983727, 'F'=>0.02284114, 'P'=>0.02654869, 'S'=>0.09236309, 'T'=>0.24591243, 'W'=>0.00561502, 'Y'=>0.01850240, 'V'=>0.07152950 },
  'W' => { 'A'=>0.03091465, 'R'=>0.02039313, 'N'=>0.01240190, 'D'=>0.01244686, 'C'=>0.01110888, 'Q'=>0.01745044, 'E'=>0.02032907, 'G'=>0.03126701, 'H'=>0.01164081, 'I'=>0.02776959, 'L'=>0.05618322, 'K'=>0.02087393, 'M'=>0.01525158, 'F'=>0.06517086, 'P'=>0.01086125, 'S'=>0.02205052, 'T'=>0.02193070, 'W'=>0.49654278, 'Y'=>0.06810771, 'V'=>0.02730512 },
  'Y' => { 'A'=>0.04027065, 'R'=>0.02869586, 'N'=>0.02169924, 'D'=>0.01853760, 'C'=>0.01071933, 'Q'=>0.02093590, 'E'=>0.02696405, 'G'=>0.02585608, 'H'=>0.04712840, 'I'=>0.04281441, 'L'=>0.06845007, 'K'=>0.03094920, 'M'=>0.01770215, 'F'=>0.13131981, 'P'=>0.01400940, 'S'=>0.03190636, 'T'=>0.02916878, 'W'=>0.02749068, 'Y'=>0.31739843, 'V'=>0.04798361 },
  'V' => { 'A'=>0.06950040, 'R'=>0.02168180, 'N'=>0.01647582, 'D'=>0.01804517, 'C'=>0.01865989, 'Q'=>0.01599165, 'E'=>0.02329673, 'G'=>0.02498419, 'H'=>0.00889793, 'I'=>0.16419105, 'L'=>0.12998823, 'K'=>0.02655057, 'M'=>0.03171098, 'F'=>0.03533099, 'P'=>0.01707551, 'S'=>0.03234722, 'T'=>0.04992154, 'W'=>0.00487917, 'Y'=>0.02124249, 'V'=>0.26922867 }
);

1;
