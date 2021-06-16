package PPH::Gene;
use strict;
use warnings;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Gene - Map protein sequence to a known gene transcript and retrieve
            nucleotide sequence-based annotations for the SNP position

=head1 DESCRIPTION

Uses database of precomputed sequence mappings for known UniProtKB proteins or
searches for near-exact matches to UniProtKB sequences using UCSC BLAT

=head1 SUBVERSION

 $LastChangedDate: 2009-11-30 00:39:17 -0500 (Mon, 30 Nov 2009) $
 $LastChangedRevision: 261 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use Carp qw(cluck confess);

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

use DBI;
use PPH::Seq;

# 0 - transition, 1 - transversion
my %TRNV = (
  'A' => { 'A' => undef, 'C' => 1, 'G' => 0, 'T' => 1 },
  'C' => { 'A' => 1, 'C' => undef,'G' => 1, 'T' => 0  },
  'G' => { 'A' => 0, 'C' => 1, 'G' => undef, 'T' => 1 },
  'T' => { 'A' => 1, 'C' => 0, 'G' => 1, 'T' => undef }
);

# Codon usage tables.
# fields: [triplet] [frequency: per thousand]
my %CDNUSE = (
  'human' =>
  {
    qw{ TTT 17.6  TCT 15.2  TAT 12.2  TGT 10.6
        TTC 20.3  TCC 17.7  TAC 15.3  TGC 12.6
        TTA  7.7  TCA 12.2  TAA  1.0  TGA  1.6
        TTG 12.9  TCG  4.4  TAG  0.8  TGG 13.2

        CTT 13.2  CCT 17.5  CAT 10.9  CGT  4.5
        CTC 19.6  CCC 19.8  CAC 15.1  CGC 10.4
        CTA  7.2  CCA 16.9  CAA 12.3  CGA  6.2
        CTG 39.6  CCG  6.9  CAG 34.2  CGG 11.4

        ATT 16.0  ACT 13.1  AAT 17.0  AGT 12.1
        ATC 20.8  ACC 18.9  AAC 19.1  AGC 19.5
        ATA  7.5  ACA 15.1  AAA 24.4  AGA 12.2
        ATG 22.0  ACG  6.1  AAG 31.9  AGG 12.0

        GTT 11.0  GCT 18.4  GAT 21.8  GGT 10.8
        GTC 14.5  GCC 27.7  GAC 25.1  GGC 22.2
        GTA  7.1  GCA 15.8  GAA 29.0  GGA 16.5
        GTG 28.1  GCG  7.4  GAG 39.6  GGG 16.5
    }
  },
);

my @TRANS = (
  # Standard genetic code
  {
    'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
    'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
    'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
    'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',

    'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
    'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
    'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
    'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',

    'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
    'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
    'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
    'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',

    'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
    'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
    'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
    'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
  },
  # Vertebrate mitochondrial code
  #
  # Differences from the Standard code:
  #         Code 2        Standard
  #  AGA    Ter  *          Arg  R
  #  AGG    Ter  *          Arg  R
  #  ATA    Met  M          Ile  I
  #  TGA    Trp  W          Ter  *
  {
    'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
    'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
    'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
    'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',

    'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
    'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
    'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
    'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',

    'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
    'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
    'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => '*',
    'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => '*',

    'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
    'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
    'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
    'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
  },
);

# AA pairs with their codons differing in a single nucleotide. Assumes that
# only AA1, CODON1 and AA2 are known a priory, so there is some level of
# synonymous codon uncertanty.
my @AASINGLE = (
  # For standard genetic code
  {
    'A' => { 'GCT' => { 'D' => [ 'GAT' ], 'G' => [ 'GGT' ], 'P' => [ 'CCT' ], 'S' => [ 'TCT' ], 'T' => [ 'ACT' ], 'V' => [ 'GTT' ] },
             'GCC' => { 'D' => [ 'GAC' ], 'G' => [ 'GGC' ], 'P' => [ 'CCC' ], 'S' => [ 'TCC' ], 'T' => [ 'ACC' ], 'V' => [ 'GTC' ] },
             'GCA' => { 'E' => [ 'GAA' ], 'G' => [ 'GGA' ], 'P' => [ 'CCA' ], 'S' => [ 'TCA' ], 'T' => [ 'ACA' ], 'V' => [ 'GTA' ] },
             'GCG' => { 'E' => [ 'GAG' ], 'G' => [ 'GGG' ], 'P' => [ 'CCG' ], 'S' => [ 'TCG' ], 'T' => [ 'ACG' ], 'V' => [ 'GTG' ] } },
    'C' => { 'TGT' => { 'F' => [ 'TTT' ], 'G' => [ 'GGT' ], 'R' => [ 'CGT' ], 'S' => [ 'AGT', 'TCT' ], 'W' => [ 'TGG' ], 'Y' => [ 'TAT' ] },
             'TGC' => { 'F' => [ 'TTC' ], 'G' => [ 'GGC' ], 'R' => [ 'CGC' ], 'S' => [ 'AGC', 'TCC' ], 'W' => [ 'TGG' ], 'Y' => [ 'TAC' ] } },
    'D' => { 'GAT' => { 'A' => [ 'GCT' ], 'E' => [ 'GAA', 'GAG' ], 'G' => [ 'GGT' ], 'H' => [ 'CAT' ], 'N' => [ 'AAT' ], 'V' => [ 'GTT' ], 'Y' => [ 'TAT' ] },
             'GAC' => { 'A' => [ 'GCC' ], 'E' => [ 'GAA', 'GAG' ], 'G' => [ 'GGC' ], 'H' => [ 'CAC' ], 'N' => [ 'AAC' ], 'V' => [ 'GTC' ], 'Y' => [ 'TAC' ] } },
    'E' => { 'GAA' => { 'A' => [ 'GCA' ], 'D' => [ 'GAC', 'GAT' ], 'G' => [ 'GGA' ], 'K' => [ 'AAA' ], 'Q' => [ 'CAA' ], 'V' => [ 'GTA' ] },
             'GAG' => { 'A' => [ 'GCG' ], 'D' => [ 'GAC', 'GAT' ], 'G' => [ 'GGG' ], 'K' => [ 'AAG' ], 'Q' => [ 'CAG' ], 'V' => [ 'GTG' ] } },
    'F' => { 'TTT' => { 'C' => [ 'TGT' ], 'I' => [ 'ATT' ], 'L' => [ 'CTT', 'TTA', 'TTG' ], 'S' => [ 'TCT' ], 'V' => [ 'GTT' ], 'Y' => [ 'TAT' ] },
             'TTC' => { 'C' => [ 'TGC' ], 'I' => [ 'ATC' ], 'L' => [ 'CTC', 'TTA', 'TTG' ], 'S' => [ 'TCC' ], 'V' => [ 'GTC' ], 'Y' => [ 'TAC' ] } },
    'G' => { 'GGT' => { 'A' => [ 'GCT' ], 'C' => [ 'TGT' ], 'D' => [ 'GAT' ], 'R' => [ 'CGT' ], 'S' => [ 'AGT' ], 'V' => [ 'GTT' ] },
             'GGC' => { 'A' => [ 'GCC' ], 'C' => [ 'TGC' ], 'D' => [ 'GAC' ], 'R' => [ 'CGC' ], 'S' => [ 'AGC' ], 'V' => [ 'GTC' ] },
             'GGA' => { 'A' => [ 'GCA' ], 'E' => [ 'GAA' ], 'R' => [ 'AGA', 'CGA' ], 'V' => [ 'GTA' ] },
             'GGG' => { 'A' => [ 'GCG' ], 'E' => [ 'GAG' ], 'R' => [ 'AGG', 'CGG' ], 'V' => [ 'GTG' ], 'W' => [ 'TGG' ] } },
    'H' => { 'CAT' => { 'D' => [ 'GAT' ], 'L' => [ 'CTT' ], 'N' => [ 'AAT' ], 'P' => [ 'CCT' ], 'Q' => [ 'CAA', 'CAG' ], 'R' => [ 'CGT' ], 'Y' => [ 'TAT' ] },
             'CAC' => { 'D' => [ 'GAC' ], 'L' => [ 'CTC' ], 'N' => [ 'AAC' ], 'P' => [ 'CCC' ], 'Q' => [ 'CAA', 'CAG' ], 'R' => [ 'CGC' ], 'Y' => [ 'TAC' ] } },
    'I' => { 'ATT' => { 'F' => [ 'TTT' ], 'L' => [ 'CTT' ], 'M' => [ 'ATG' ], 'N' => [ 'AAT' ], 'S' => [ 'AGT' ], 'T' => [ 'ACT' ], 'V' => [ 'GTT' ] },
             'ATC' => { 'F' => [ 'TTC' ], 'L' => [ 'CTC' ], 'M' => [ 'ATG' ], 'N' => [ 'AAC' ], 'S' => [ 'AGC' ], 'T' => [ 'ACC' ], 'V' => [ 'GTC' ] },
             'ATA' => { 'K' => [ 'AAA' ], 'L' => [ 'CTA', 'TTA' ], 'M' => [ 'ATG' ], 'R' => [ 'AGA' ], 'T' => [ 'ACA' ], 'V' => [ 'GTA' ] } },
    'K' => { 'AAA' => { 'E' => [ 'GAA' ], 'I' => [ 'ATA' ], 'N' => [ 'AAC', 'AAT' ], 'Q' => [ 'CAA' ], 'R' => [ 'AGA' ], 'T' => [ 'ACA' ] },
             'AAG' => { 'E' => [ 'GAG' ], 'M' => [ 'ATG' ], 'N' => [ 'AAC', 'AAT' ], 'Q' => [ 'CAG' ], 'R' => [ 'AGG' ], 'T' => [ 'ACG' ] } },
    'L' => { 'TTA' => { 'F' => [ 'TTC', 'TTT' ], 'I' => [ 'ATA' ], 'S' => [ 'TCA' ], 'V' => [ 'GTA' ] },
             'TTG' => { 'F' => [ 'TTC', 'TTT' ], 'M' => [ 'ATG' ], 'S' => [ 'TCG' ], 'V' => [ 'GTG' ], 'W' => [ 'TGG' ] },
             'CTT' => { 'F' => [ 'TTT' ], 'H' => [ 'CAT' ], 'I' => [ 'ATT' ], 'P' => [ 'CCT' ], 'R' => [ 'CGT' ], 'V' => [ 'GTT' ] },
             'CTC' => { 'F' => [ 'TTC' ], 'H' => [ 'CAC' ], 'I' => [ 'ATC' ], 'P' => [ 'CCC' ], 'R' => [ 'CGC' ], 'V' => [ 'GTC' ] },
             'CTA' => { 'I' => [ 'ATA' ], 'P' => [ 'CCA' ], 'Q' => [ 'CAA' ], 'R' => [ 'CGA' ], 'V' => [ 'GTA' ] },
             'CTG' => { 'M' => [ 'ATG' ], 'P' => [ 'CCG' ], 'Q' => [ 'CAG' ], 'R' => [ 'CGG' ], 'V' => [ 'GTG' ] } },
    'M' => { 'ATG' => { 'I' => [ 'ATA', 'ATC', 'ATT' ], 'K' => [ 'AAG' ], 'L' => [ 'CTG', 'TTG' ], 'R' => [ 'AGG' ], 'T' => [ 'ACG' ], 'V' => [ 'GTG' ] } },
    'N' => { 'AAT' => { 'D' => [ 'GAT' ], 'H' => [ 'CAT' ], 'I' => [ 'ATT' ], 'K' => [ 'AAA', 'AAG' ], 'S' => [ 'AGT' ], 'T' => [ 'ACT' ], 'Y' => [ 'TAT' ] },
             'AAC' => { 'D' => [ 'GAC' ], 'H' => [ 'CAC' ], 'I' => [ 'ATC' ], 'K' => [ 'AAA', 'AAG' ], 'S' => [ 'AGC' ], 'T' => [ 'ACC' ], 'Y' => [ 'TAC' ] } },
    'P' => { 'CCT' => { 'A' => [ 'GCT' ], 'H' => [ 'CAT' ], 'L' => [ 'CTT' ], 'R' => [ 'CGT' ], 'S' => [ 'TCT' ], 'T' => [ 'ACT' ] },
             'CCC' => { 'A' => [ 'GCC' ], 'H' => [ 'CAC' ], 'L' => [ 'CTC' ], 'R' => [ 'CGC' ], 'S' => [ 'TCC' ], 'T' => [ 'ACC' ] },
             'CCA' => { 'A' => [ 'GCA' ], 'L' => [ 'CTA' ], 'Q' => [ 'CAA' ], 'R' => [ 'CGA' ], 'S' => [ 'TCA' ], 'T' => [ 'ACA' ] },
             'CCG' => { 'A' => [ 'GCG' ], 'L' => [ 'CTG' ], 'Q' => [ 'CAG' ], 'R' => [ 'CGG' ], 'S' => [ 'TCG' ], 'T' => [ 'ACG' ] } },
    'Q' => { 'CAA' => { 'E' => [ 'GAA' ], 'H' => [ 'CAC', 'CAT' ], 'K' => [ 'AAA' ], 'L' => [ 'CTA' ], 'P' => [ 'CCA' ], 'R' => [ 'CGA' ] },
             'CAG' => { 'E' => [ 'GAG' ], 'H' => [ 'CAC', 'CAT' ], 'K' => [ 'AAG' ], 'L' => [ 'CTG' ], 'P' => [ 'CCG' ], 'R' => [ 'CGG' ] } },
    'R' => { 'CGT' => { 'C' => [ 'TGT' ], 'G' => [ 'GGT' ], 'H' => [ 'CAT' ], 'L' => [ 'CTT' ], 'P' => [ 'CCT' ], 'S' => [ 'AGT' ] },
             'CGC' => { 'C' => [ 'TGC' ], 'G' => [ 'GGC' ], 'H' => [ 'CAC' ], 'L' => [ 'CTC' ], 'P' => [ 'CCC' ], 'S' => [ 'AGC' ] },
             'CGA' => { 'G' => [ 'GGA' ], 'L' => [ 'CTA' ], 'P' => [ 'CCA' ], 'Q' => [ 'CAA' ] },
             'CGG' => { 'G' => [ 'GGG' ], 'L' => [ 'CTG' ], 'P' => [ 'CCG' ], 'Q' => [ 'CAG' ], 'W' => [ 'TGG' ] },
             'AGA' => { 'G' => [ 'GGA' ], 'I' => [ 'ATA' ], 'K' => [ 'AAA' ], 'S' => [ 'AGC', 'AGT' ], 'T' => [ 'ACA' ] },
             'AGG' => { 'G' => [ 'GGG' ], 'K' => [ 'AAG' ], 'M' => [ 'ATG' ], 'S' => [ 'AGC', 'AGT' ], 'T' => [ 'ACG' ], 'W' => [ 'TGG' ] } },
    'S' => { 'TCT' => { 'A' => [ 'GCT' ], 'C' => [ 'TGT' ], 'F' => [ 'TTT' ], 'P' => [ 'CCT' ], 'T' => [ 'ACT' ], 'Y' => [ 'TAT' ] },
             'TCC' => { 'A' => [ 'GCC' ], 'C' => [ 'TGC' ], 'F' => [ 'TTC' ], 'P' => [ 'CCC' ], 'T' => [ 'ACC' ], 'Y' => [ 'TAC' ] },
             'TCA' => { 'A' => [ 'GCA' ], 'L' => [ 'TTA' ], 'P' => [ 'CCA' ], 'T' => [ 'ACA' ] },
             'TCG' => { 'A' => [ 'GCG' ], 'L' => [ 'TTG' ], 'P' => [ 'CCG' ], 'T' => [ 'ACG' ], 'W' => [ 'TGG' ] },
             'AGT' => { 'C' => [ 'TGT' ], 'G' => [ 'GGT' ], 'I' => [ 'ATT' ], 'N' => [ 'AAT' ], 'R' => [ 'AGA', 'AGG', 'CGT' ], 'T' => [ 'ACT' ] },
             'AGC' => { 'C' => [ 'TGC' ], 'G' => [ 'GGC' ], 'I' => [ 'ATC' ], 'N' => [ 'AAC' ], 'R' => [ 'AGA', 'AGG', 'CGC' ], 'T' => [ 'ACC' ] } },
    'T' => { 'ACT' => { 'A' => [ 'GCT' ], 'I' => [ 'ATT' ], 'N' => [ 'AAT' ], 'P' => [ 'CCT' ], 'S' => [ 'AGT', 'TCT' ] },
             'ACC' => { 'A' => [ 'GCC' ], 'I' => [ 'ATC' ], 'N' => [ 'AAC' ], 'P' => [ 'CCC' ], 'S' => [ 'AGC', 'TCC' ] },
             'ACA' => { 'A' => [ 'GCA' ], 'I' => [ 'ATA' ], 'K' => [ 'AAA' ], 'P' => [ 'CCA' ], 'R' => [ 'AGA' ], 'S' => [ 'TCA' ] },
             'ACG' => { 'A' => [ 'GCG' ], 'K' => [ 'AAG' ], 'M' => [ 'ATG' ], 'P' => [ 'CCG' ], 'R' => [ 'AGG' ], 'S' => [ 'TCG' ] } },
    'V' => { 'GTT' => { 'A' => [ 'GCT' ], 'D' => [ 'GAT' ], 'F' => [ 'TTT' ], 'G' => [ 'GGT' ], 'I' => [ 'ATT' ], 'L' => [ 'CTT' ] },
             'GTC' => { 'A' => [ 'GCC' ], 'D' => [ 'GAC' ], 'F' => [ 'TTC' ], 'G' => [ 'GGC' ], 'I' => [ 'ATC' ], 'L' => [ 'CTC' ] },
             'GTA' => { 'A' => [ 'GCA' ], 'E' => [ 'GAA' ], 'G' => [ 'GGA' ], 'I' => [ 'ATA' ], 'L' => [ 'CTA', 'TTA' ] },
             'GTG' => { 'A' => [ 'GCG' ], 'E' => [ 'GAG' ], 'G' => [ 'GGG' ], 'L' => [ 'CTG', 'TTG' ], 'M' => [ 'ATG' ] } },
    'W' => { 'TGG' => { 'C' => [ 'TGC', 'TGT' ], 'G' => [ 'GGG' ], 'L' => [ 'TTG' ], 'R' => [ 'AGG', 'CGG' ], 'S' => [ 'TCG' ] } },
    'Y' => { 'TAT' => { 'C' => [ 'TGT' ], 'D' => [ 'GAT' ], 'F' => [ 'TTT' ], 'H' => [ 'CAT' ], 'N' => [ 'AAT' ], 'S' => [ 'TCT' ] },
             'TAC' => { 'C' => [ 'TGC' ], 'D' => [ 'GAC' ], 'F' => [ 'TTC' ], 'H' => [ 'CAC' ], 'N' => [ 'AAC' ], 'S' => [ 'TCC' ] } }
  },
  # For vertebrate mitochondrial code
  {
    'A' => { 'GCT' => { 'D' => [ 'GAT' ], 'G' => [ 'GGT' ], 'P' => [ 'CCT' ], 'S' => [ 'TCT' ], 'T' => [ 'ACT' ], 'V' => [ 'GTT' ] },
             'GCC' => { 'D' => [ 'GAC' ], 'G' => [ 'GGC' ], 'P' => [ 'CCC' ], 'S' => [ 'TCC' ], 'T' => [ 'ACC' ], 'V' => [ 'GTC' ] },
             'GCA' => { 'E' => [ 'GAA' ], 'G' => [ 'GGA' ], 'P' => [ 'CCA' ], 'S' => [ 'TCA' ], 'T' => [ 'ACA' ], 'V' => [ 'GTA' ] },
             'GCG' => { 'E' => [ 'GAG' ], 'G' => [ 'GGG' ], 'P' => [ 'CCG' ], 'S' => [ 'TCG' ], 'T' => [ 'ACG' ], 'V' => [ 'GTG' ] } },
    'C' => { 'TGT' => { 'F' => [ 'TTT' ], 'G' => [ 'GGT' ], 'R' => [ 'CGT' ], 'S' => [ 'AGT', 'TCT' ], 'W' => [ 'TGA', 'TGG' ], 'Y' => [ 'TAT' ] },
             'TGC' => { 'F' => [ 'TTC' ], 'G' => [ 'GGC' ], 'R' => [ 'CGC' ], 'S' => [ 'AGC', 'TCC' ], 'W' => [ 'TGA', 'TGG' ], 'Y' => [ 'TAC' ] } },
    'D' => { 'GAT' => { 'A' => [ 'GCT' ], 'E' => [ 'GAA', 'GAG' ], 'G' => [ 'GGT' ], 'H' => [ 'CAT' ], 'N' => [ 'AAT' ], 'V' => [ 'GTT' ], 'Y' => [ 'TAT' ] },
             'GAC' => { 'A' => [ 'GCC' ], 'E' => [ 'GAA', 'GAG' ], 'G' => [ 'GGC' ], 'H' => [ 'CAC' ], 'N' => [ 'AAC' ], 'V' => [ 'GTC' ], 'Y' => [ 'TAC' ] } },
    'E' => { 'GAA' => { 'A' => [ 'GCA' ], 'D' => [ 'GAC', 'GAT' ], 'G' => [ 'GGA' ], 'K' => [ 'AAA' ], 'Q' => [ 'CAA' ], 'V' => [ 'GTA' ] },
             'GAG' => { 'A' => [ 'GCG' ], 'D' => [ 'GAC', 'GAT' ], 'G' => [ 'GGG' ], 'K' => [ 'AAG' ], 'Q' => [ 'CAG' ], 'V' => [ 'GTG' ] } },
    'F' => { 'TTT' => { 'C' => [ 'TGT' ], 'I' => [ 'ATT' ], 'L' => [ 'CTT', 'TTA', 'TTG' ], 'S' => [ 'TCT' ], 'V' => [ 'GTT' ], 'Y' => [ 'TAT' ] },
             'TTC' => { 'C' => [ 'TGC' ], 'I' => [ 'ATC' ], 'L' => [ 'CTC', 'TTA', 'TTG' ], 'S' => [ 'TCC' ], 'V' => [ 'GTC' ], 'Y' => [ 'TAC' ] } },
    'G' => { 'GGT' => { 'A' => [ 'GCT' ], 'C' => [ 'TGT' ], 'D' => [ 'GAT' ], 'R' => [ 'CGT' ], 'S' => [ 'AGT' ], 'V' => [ 'GTT' ] },
             'GGC' => { 'A' => [ 'GCC' ], 'C' => [ 'TGC' ], 'D' => [ 'GAC' ], 'R' => [ 'CGC' ], 'S' => [ 'AGC' ], 'V' => [ 'GTC' ] },
             'GGA' => { 'A' => [ 'GCA' ], 'E' => [ 'GAA' ], 'R' => [ 'CGA' ], 'V' => [ 'GTA' ], 'W' => [ 'TGA' ] },
             'GGG' => { 'A' => [ 'GCG' ], 'E' => [ 'GAG' ], 'R' => [ 'CGG' ], 'V' => [ 'GTG' ], 'W' => [ 'TGG' ] } },
    'H' => { 'CAT' => { 'D' => [ 'GAT' ], 'L' => [ 'CTT' ], 'N' => [ 'AAT' ], 'P' => [ 'CCT' ], 'Q' => [ 'CAA', 'CAG' ], 'R' => [ 'CGT' ], 'Y' => [ 'TAT' ] },
             'CAC' => { 'D' => [ 'GAC' ], 'L' => [ 'CTC' ], 'N' => [ 'AAC' ], 'P' => [ 'CCC' ], 'Q' => [ 'CAA', 'CAG' ], 'R' => [ 'CGC' ], 'Y' => [ 'TAC' ] } },
    'I' => { 'ATT' => { 'F' => [ 'TTT' ], 'L' => [ 'CTT' ], 'M' => [ 'ATA', 'ATG' ], 'N' => [ 'AAT' ], 'S' => [ 'AGT' ], 'T' => [ 'ACT' ], 'V' => [ 'GTT' ] },
             'ATC' => { 'F' => [ 'TTC' ], 'L' => [ 'CTC' ], 'M' => [ 'ATA', 'ATG' ], 'N' => [ 'AAC' ], 'S' => [ 'AGC' ], 'T' => [ 'ACC' ], 'V' => [ 'GTC' ] } },
    'K' => { 'AAA' => { 'E' => [ 'GAA' ], 'M' => [ 'ATA' ], 'N' => [ 'AAC', 'AAT' ], 'Q' => [ 'CAA' ], 'T' => [ 'ACA' ] },
             'AAG' => { 'E' => [ 'GAG' ], 'M' => [ 'ATG' ], 'N' => [ 'AAC', 'AAT' ], 'Q' => [ 'CAG' ], 'T' => [ 'ACG' ] } },
    'L' => { 'TTA' => { 'F' => [ 'TTC', 'TTT' ], 'M' => [ 'ATA' ], 'S' => [ 'TCA' ], 'V' => [ 'GTA' ], 'W' => [ 'TGA' ] },
             'TTG' => { 'F' => [ 'TTC', 'TTT' ], 'M' => [ 'ATG' ], 'S' => [ 'TCG' ], 'V' => [ 'GTG' ], 'W' => [ 'TGG' ] },
             'CTT' => { 'F' => [ 'TTT' ], 'H' => [ 'CAT' ], 'I' => [ 'ATT' ], 'P' => [ 'CCT' ], 'R' => [ 'CGT' ], 'V' => [ 'GTT' ] },
             'CTC' => { 'F' => [ 'TTC' ], 'H' => [ 'CAC' ], 'I' => [ 'ATC' ], 'P' => [ 'CCC' ], 'R' => [ 'CGC' ], 'V' => [ 'GTC' ] },
             'CTA' => { 'M' => [ 'ATA' ], 'P' => [ 'CCA' ], 'Q' => [ 'CAA' ], 'R' => [ 'CGA' ], 'V' => [ 'GTA' ] },
             'CTG' => { 'M' => [ 'ATG' ], 'P' => [ 'CCG' ], 'Q' => [ 'CAG' ], 'R' => [ 'CGG' ], 'V' => [ 'GTG' ] } },
    'M' => { 'ATA' => { 'I' => [ 'ATC', 'ATT' ], 'K' => [ 'AAA' ], 'L' => [ 'CTA', 'TTA' ], 'T' => [ 'ACA' ], 'V' => [ 'GTA' ] },
             'ATG' => { 'I' => [ 'ATC', 'ATT' ], 'K' => [ 'AAG' ], 'L' => [ 'CTG', 'TTG' ], 'T' => [ 'ACG' ], 'V' => [ 'GTG' ] } },
    'N' => { 'AAT' => { 'D' => [ 'GAT' ], 'H' => [ 'CAT' ], 'I' => [ 'ATT' ], 'K' => [ 'AAA', 'AAG' ], 'S' => [ 'AGT' ], 'T' => [ 'ACT' ], 'Y' => [ 'TAT' ] },
             'AAC' => { 'D' => [ 'GAC' ], 'H' => [ 'CAC' ], 'I' => [ 'ATC' ], 'K' => [ 'AAA', 'AAG' ], 'S' => [ 'AGC' ], 'T' => [ 'ACC' ], 'Y' => [ 'TAC' ] } },
    'P' => { 'CCT' => { 'A' => [ 'GCT' ], 'H' => [ 'CAT' ], 'L' => [ 'CTT' ], 'R' => [ 'CGT' ], 'S' => [ 'TCT' ], 'T' => [ 'ACT' ] },
             'CCC' => { 'A' => [ 'GCC' ], 'H' => [ 'CAC' ], 'L' => [ 'CTC' ], 'R' => [ 'CGC' ], 'S' => [ 'TCC' ], 'T' => [ 'ACC' ] },
             'CCA' => { 'A' => [ 'GCA' ], 'L' => [ 'CTA' ], 'Q' => [ 'CAA' ], 'R' => [ 'CGA' ], 'S' => [ 'TCA' ], 'T' => [ 'ACA' ] },
             'CCG' => { 'A' => [ 'GCG' ], 'L' => [ 'CTG' ], 'Q' => [ 'CAG' ], 'R' => [ 'CGG' ], 'S' => [ 'TCG' ], 'T' => [ 'ACG' ] } },
    'Q' => { 'CAA' => { 'E' => [ 'GAA' ], 'H' => [ 'CAC', 'CAT' ], 'K' => [ 'AAA' ], 'L' => [ 'CTA' ], 'P' => [ 'CCA' ], 'R' => [ 'CGA' ] },
             'CAG' => { 'E' => [ 'GAG' ], 'H' => [ 'CAC', 'CAT' ], 'K' => [ 'AAG' ], 'L' => [ 'CTG' ], 'P' => [ 'CCG' ], 'R' => [ 'CGG' ] } },
    'R' => { 'CGT' => { 'C' => [ 'TGT' ], 'G' => [ 'GGT' ], 'H' => [ 'CAT' ], 'L' => [ 'CTT' ], 'P' => [ 'CCT' ], 'S' => [ 'AGT' ] },
             'CGC' => { 'C' => [ 'TGC' ], 'G' => [ 'GGC' ], 'H' => [ 'CAC' ], 'L' => [ 'CTC' ], 'P' => [ 'CCC' ], 'S' => [ 'AGC' ] },
             'CGA' => { 'G' => [ 'GGA' ], 'L' => [ 'CTA' ], 'P' => [ 'CCA' ], 'Q' => [ 'CAA' ], 'W' => [ 'TGA' ] },
             'CGG' => { 'G' => [ 'GGG' ], 'L' => [ 'CTG' ], 'P' => [ 'CCG' ], 'Q' => [ 'CAG' ], 'W' => [ 'TGG' ] } },
    'S' => { 'TCT' => { 'A' => [ 'GCT' ], 'C' => [ 'TGT' ], 'F' => [ 'TTT' ], 'P' => [ 'CCT' ], 'T' => [ 'ACT' ], 'Y' => [ 'TAT' ] },
             'TCC' => { 'A' => [ 'GCC' ], 'C' => [ 'TGC' ], 'F' => [ 'TTC' ], 'P' => [ 'CCC' ], 'T' => [ 'ACC' ], 'Y' => [ 'TAC' ] },
             'TCA' => { 'A' => [ 'GCA' ], 'L' => [ 'TTA' ], 'P' => [ 'CCA' ], 'T' => [ 'ACA' ], 'W' => [ 'TGA' ] },
             'TCG' => { 'A' => [ 'GCG' ], 'L' => [ 'TTG' ], 'P' => [ 'CCG' ], 'T' => [ 'ACG' ], 'W' => [ 'TGG' ] },
             'AGT' => { 'C' => [ 'TGT' ], 'G' => [ 'GGT' ], 'I' => [ 'ATT' ], 'N' => [ 'AAT' ], 'R' => [ 'CGT' ], 'T' => [ 'ACT' ] },
             'AGC' => { 'C' => [ 'TGC' ], 'G' => [ 'GGC' ], 'I' => [ 'ATC' ], 'N' => [ 'AAC' ], 'R' => [ 'CGC' ], 'T' => [ 'ACC' ] } },
    'T' => { 'ACT' => { 'A' => [ 'GCT' ], 'I' => [ 'ATT' ], 'N' => [ 'AAT' ], 'P' => [ 'CCT' ], 'S' => [ 'AGT', 'TCT' ] },
             'ACC' => { 'A' => [ 'GCC' ], 'I' => [ 'ATC' ], 'N' => [ 'AAC' ], 'P' => [ 'CCC' ], 'S' => [ 'AGC', 'TCC' ] },
             'ACA' => { 'A' => [ 'GCA' ], 'K' => [ 'AAA' ], 'M' => [ 'ATA' ], 'P' => [ 'CCA' ], 'S' => [ 'TCA' ] },
             'ACG' => { 'A' => [ 'GCG' ], 'K' => [ 'AAG' ], 'M' => [ 'ATG' ], 'P' => [ 'CCG' ], 'S' => [ 'TCG' ] } },
    'V' => { 'GTT' => { 'A' => [ 'GCT' ], 'D' => [ 'GAT' ], 'F' => [ 'TTT' ], 'G' => [ 'GGT' ], 'I' => [ 'ATT' ], 'L' => [ 'CTT' ] },
             'GTC' => { 'A' => [ 'GCC' ], 'D' => [ 'GAC' ], 'F' => [ 'TTC' ], 'G' => [ 'GGC' ], 'I' => [ 'ATC' ], 'L' => [ 'CTC' ] },
             'GTA' => { 'A' => [ 'GCA' ], 'E' => [ 'GAA' ], 'G' => [ 'GGA' ], 'L' => [ 'CTA', 'TTA' ], 'M' => [ 'ATA' ] },
             'GTG' => { 'A' => [ 'GCG' ], 'E' => [ 'GAG' ], 'G' => [ 'GGG' ], 'L' => [ 'CTG', 'TTG' ], 'M' => [ 'ATG' ] } },
    'W' => { 'TGA' => { 'C' => [ 'TGC', 'TGT' ], 'G' => [ 'GGA' ], 'L' => [ 'TTA' ], 'R' => [ 'CGA' ], 'S' => [ 'TCA' ] },
             'TGG' => { 'C' => [ 'TGC', 'TGT' ], 'G' => [ 'GGG' ], 'L' => [ 'TTG' ], 'R' => [ 'CGG' ], 'S' => [ 'TCG' ] } },
    'Y' => { 'TAT' => { 'C' => [ 'TGT' ], 'D' => [ 'GAT' ], 'F' => [ 'TTT' ], 'H' => [ 'CAT' ], 'N' => [ 'AAT' ], 'S' => [ 'TCT' ] },
             'TAC' => { 'C' => [ 'TGC' ], 'D' => [ 'GAC' ], 'F' => [ 'TTC' ], 'H' => [ 'CAC' ], 'N' => [ 'AAC' ], 'S' => [ 'TCC' ] } }
  },
);

# UniProtKB->knownGene sequence mappings database
my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
my $dbMap  = "$CONFIG{GOLDENPATH}/$CONFIG{GENESET}/genes/upToKg.sqlite";
my $dbGene = "$CONFIG{GOLDENPATH}/$CONFIG{GENESET}/genes/knownGene.sqlite";

my $selectMap = q{
SELECT macc, qfrom, qto, mfrom, mto
FROM map JOIN hits USING(id)
WHERE qacc = ? AND ? BETWEEN qfrom AND qto
ORDER BY mcan DESC, qfrac DESC, hits ASC, ident DESC LIMIT 1
};
my $selectGene = q{SELECT * FROM knownGene WHERE name = ?};
my $selectExChr = q{SELECT exStart, exEnd FROM kgExChr WHERE id = ? ORDER BY exStart ASC};
my $selectExCds = q{SELECT exFrom, exTo FROM kgExCds WHERE id = ? ORDER BY exFrom ASC};

my ($dbhMap, $dbhGene);
my ($sthMap, $sthGene, $sthExChr, $sthExCds);

sub find_gene {
  my ($PROT, $SNP) = @_;
  #### PPH__Gene__find_gene: @_

  # Search for a known UniProtKB protein match in the precomputed database
  my ($txname, $cdnpos);
  if ($PROT->{UniProt}) {
    # Connect to database and prepare select statement if this has not been done yet
    unless (defined $dbhMap) {
      $dbhMap = DBI->connect("dbi:SQLite:dbname=$dbMap", '', '', $dbargs);
      $sthMap = $dbhMap->prepare($selectMap);
    }
    # Search table
    $sthMap->execute($SNP->{Acc}, $SNP->{Pos} - 1);
    # Fetch results (should be a single row)
    my ($macc, $qfrom, $qto, $mfrom, $mto) = $sthMap->fetchrow_array;
    $sthMap->finish;
    return 0 unless defined $macc;
    $txname = $macc;
    $cdnpos = $SNP->{Pos} + $mfrom - $qfrom;
    ##### |     find_gene found (TxName, CdnPos): $txname, $cdnpos
  # Try to map user-specified protein sequence to a translated transcript CDS
  } else {
    my $dbFASTA = "$CONFIG{GOLDENPATH}/$CONFIG{GENESET}/genes/knownGeneAA.seq";
    my $offset;
    my ($macc, $mname, $mdesc, $mseq) = PPH::Seq::grep_in_FASTA($PROT->{Seq}, $dbFASTA);
    if (defined $macc) {
      $offset = index($mseq, $PROT->{Seq});
      $txname = $macc;
      ##### |     find_gene exact grep hit (offset): $offset
    # Try BLAT if grep failed
    } else {
      # Fetch single top hit only
      my $hits = PPH::Seq::map_by_BLAT($PROT->{Acc}, $PROT->{Seq}, $dbFASTA);
      return 0 unless defined $hits;
      ##### |     find_gene BLAT hits: $hits
      my ($qacc, $macc, $qlen, $mlen, $overlap, $blocks, $qstart, $qend, $mstart, $mend,
          $qfrac, $ident, $mcan, $bsizes, $qstarts, $mstarts) = @{ $hits->[0] };
      my $qpos = $SNP->{Pos} - 1;
      my @bsizes  = split /,/, $bsizes;
      my @qstarts = split /,/, $qstarts;
      my @mstarts = split /,/, $mstarts;
      for (my $i=0; $i<@qstarts; $i++) {
        if ($qpos >= $qstarts[$i] && $qpos < $qstarts[$i] + $bsizes[$i]) {
          $offset = $mstarts[$i] - $qstarts[$i];
          last;
        }
      }
      $txname = $macc;
    }
    return 0 unless defined $offset;
    $cdnpos = $SNP->{Pos} + $offset;
    ##### |     find_gene mapped (TxName, CdnPos): $txname, $cdnpos
  }

  # Now fetch the transcript nucleotide sequence and annotations, unless already cached
  unless (exists $PROT->{Gene}{$txname}) {
    my $TX;

    # The empty hashref stub indicating fetching the data was attempted.
    # Replaced with ref to valid data if processing completes successfully.
    $PROT->{Gene}{$txname} = {};

    # Fetch complete set of annotations for the transcript
    unless (defined $dbhGene) {
      $dbhGene = DBI->connect("dbi:SQLite:dbname=$dbGene", '', '', $dbargs);
      $sthGene = $dbhGene->prepare($selectGene);
    }
    $sthGene->execute($txname);
    $TX = $sthGene->fetchrow_hashref;
    $sthGene->finish;
    return 0 unless defined $TX;

    my $chrom = $TX->{chrom};
    unless ($chrom =~ /^chr\d{1,2}$/ || $chrom =~ /^chr[MXY]$/) {
      warn "ERROR: find_gene: Non-canonical chromosome assembly ($chrom) encountered for transcript: $txname\n";
      return 0;
    }

    # Fetch full transcript nucleotide sequence: [txStart..txEnd), introns included
    my $nucseq = PPH::Seq::fetch_txNuc($txname, "$CONFIG{GOLDENPATH}/$CONFIG{GENESET}/genes/knownGeneNuc.2bit");
    return 0 unless defined $nucseq;  # nothing to do if we do not have a sequence
    $TX->{TxSeq} = $nucseq;
    $TX->{TxLen} = length $nucseq;

    $sthExCds = $dbhGene->prepare($selectExCds) unless defined $sthExCds;
    $sthExCds->execute($TX->{id});
    my @cds2tx; # maps nucleotide positions in CDS to nucleotide positions in the full transcript sequence
    while (my $cdsExon = $sthExCds->fetchrow_hashref) {
      push @cds2tx, $cdsExon->{exFrom} .. $cdsExon->{exTo};
    }
    $sthExCds->finish;
    my $cdslen = scalar @cds2tx;
    if ($cdslen % 3) {
      warn "ERROR: find_gene: Incomplete CDS for transcript: $txname\n";
      return 0;
    }
    $TX->{Cds2Tx} = \@cds2tx;

    # Get full list of exons, used later to find exon/intron junctions.
    # Only makes sense if the transcript has at least one intron.
    my @chrExStarts;
    my @chrExEnds;
    if ($TX->{exonCount} > 1) {
      $sthExChr = $dbhGene->prepare($selectExChr) unless defined $sthExChr;
      $sthExChr->execute($TX->{id});
      while (my $exon = $sthExChr->fetchrow_hashref) {
        push @chrExStarts, $exon->{exStart};
        push @chrExEnds, $exon->{exEnd};
      }
      $sthExChr->finish;
    }
    $TX->{chrExStarts} = \@chrExStarts;
    $TX->{chrExEnds}   = \@chrExEnds;

    ##### |     find_gene cached: $txname, $TX
    $PROT->{Gene}{$txname} = $TX;
  }

  # Process nucleotide-based features for the SNP if data for the corresponding
  # transcript was found and successfully loaded.
  my $TX = $PROT->{Gene}{$txname};
  if (keys %$TX && defined $cdnpos) {
    $SNP->{Gene}{TxName} = $txname;       # transcript id
    $SNP->{Gene}{CdnPos} = $cdnpos;       # codon (= AA residue) position in CDS, base one
    my $cdspos = ($cdnpos - 1) * 3;       # first codon nucleotide position in CDS, base zero
    my $txpos  = $TX->{Cds2Tx}[$cdspos];  # first codon nucleotide position in full transcript, base zero
    my $cdn1  =
      substr($TX->{TxSeq}, $TX->{Cds2Tx}[$cdspos],     1) .
      substr($TX->{TxSeq}, $TX->{Cds2Tx}[$cdspos + 1], 1) .
      substr($TX->{TxSeq}, $TX->{Cds2Tx}[$cdspos + 2], 1);
    # Standard genetics code index is 0, vertebrate mitochondrial code is 1
    my $gencode = $TX->{chrom} eq 'chrM' ? 1 : 0;
    my $aa1 = $TRANS[$gencode]{$cdn1};
    unless (defined $aa1) {
      warn "ERROR: find_gene: Illegal codon ($cdn1) found in $txname nucleotide sequence at position: ".($cdspos+1)."\n";
      return 0;
    }
    # Which of the two AA residues / codons we expect to be encoded on the transcript sequence?
    # Depends on the value of REVERSEDIRECTION and whether codon in the transcript sequence
    # matches AA1 or AA2.
    my $txaakey1 = $CONFIG{'REVERSEDIRECTION'} ? 'Aa2' : 'Aa1';
    my $txaakey2 = $txaakey1 eq 'Aa1' ? 'Aa2' : 'Aa1';
    my $swapcodons;
    if ($aa1 eq $SNP->{$txaakey1}) {
      $swapcodons = 0;
    } elsif ($aa1 eq $SNP->{$txaakey2}) {
      ($txaakey1, $txaakey2) = ($txaakey2, $txaakey1);
      $swapcodons = 1;
    }
    unless (exists $AASINGLE[$gencode]{$SNP->{$txaakey1}}{$cdn1}) {
      warn "ERROR: find_gene: Incorrect codon ($cdn1) specified for ".uc($txaakey1)."($SNP->{$txaakey1})\n";
      return 0;
    }
    unless (exists $AASINGLE[$gencode]{$SNP->{$txaakey1}}{$cdn1}{$SNP->{$txaakey2}}) {
      warn "ERROR: find_gene: Not a single-nucleotide AA substitution: $SNP->{Aa1}($cdn1)>$SNP->{Aa2}\n";
      return 0;
    }
    my @codons2 = @{ $AASINGLE[$gencode]{$SNP->{$txaakey1}}{$cdn1}{$SNP->{$txaakey2}} };
    unless (@codons2) {
      confess "ERROR: find_gene: Unable to locate codon(s) for the substitution path: $SNP->{$txaakey1}($cdn1)>$SNP->{$txaakey2}";
    }
    my $cdn2;
    if (@codons2 > 1) {
      ##### |     find_gene codons2: @codons2
      $cdn2 = (sort { $CDNUSE{$CONFIG{REFORGCOMMON}}{$b} <=> $CDNUSE{$CONFIG{REFORGCOMMON}}{$a}} @codons2)[0];
    } else {
      $cdn2 = $codons2[0];
    }
    my $frame;
    for (my $i=0; $i<length($cdn1); $i++) {
      if (substr($cdn1, $i, 1) ne substr($cdn2, $i, 1)) {
        $frame = $i;
        last;
      }
    }
    # Swap codons if codon in the transcript sequence actually encodes AA2
    if ($swapcodons) {
      warn "WARNING: find_gene: Swapped codons ($cdn1>$cdn2) in $txname nucleotide sequence at position: ".($cdspos+1)."\n";
      ($cdn1, $cdn2) = ($cdn2, $cdn1);
    }
    if (defined $frame) {
      $SNP->{Gene}{Frame} = $frame;
    } else {
      warn "ERROR: find_gene: Synonymous substitution: $SNP->{Aa1}($cdn1)>$SNP->{Aa2}($cdn2)\n";
      return 0;
    }
    ###### |     find_gene reconstructed substitution path: "$SNP->{Aa1}($cdn1)>$SNP->{Aa2}($cdn2)"
    $SNP->{Gene}{Cdn1} = $cdn1;
    $SNP->{Gene}{Cdn2} = $cdn2;
    # mutation nucleotide position in full transcript, base zero
    my $mutpos = $txpos + $frame;
    # transcript (reference) allele nucleotide
    my $nt1 = substr($cdn1, $frame, 1);
    # mutant allele nucleotide
    my $nt2 = substr($cdn2, $frame, 1);
    # nucleotides flanking the mutation position (5'3')
    my $f5 = substr($TX->{TxSeq}, $mutpos - 1, 1) || '?';
    my $f3 = substr($TX->{TxSeq}, $mutpos + 1, 1) || '?';
    $SNP->{Gene}{Nt1}    = $nt1;
    $SNP->{Gene}{Nt2}    = $nt2;
    $SNP->{Gene}{Flanks} = $f5 . $f3;
    # mutation nucleotide position in full transcript, base one
    $SNP->{Gene}{TxPos}  = $mutpos + 1;
    # mutation nucleotide position on the chromosome assembly, base zero
    my $chrpos = $TX->{strand} eq '+' ? $TX->{txStart} + $mutpos : $TX->{txEnd} - $mutpos - 1;
    # mutation nucleotide position on the chromosome assembly, base one
    $SNP->{Gene}{ChrPos} = $chrpos + 1;
    $SNP->{Gene}{Transv} = $TRNV{$nt1}{$nt2} if defined $TRNV{$nt1}{$nt2};
    $SNP->{Gene}{CpG}    = CpG($f5, $nt1, $nt2, $f3);
    # Locate nearest exon/intron boundaries
    # Only makes sense if the transcript has at least one intron
    if ($TX->{exonCount} > 1) {
      my ($dstart, $dend) = ($TX->{chrExEnds}[-1], $TX->{chrExStarts}[0]);
      for (my $i=0; $i<@{$TX->{chrExStarts}}; $i++) {
        my $start = $TX->{chrExStarts}[$i];
        my $end   = $TX->{chrExEnds}[$i];
        # Skip exon start if it is at the transcript start -- no junction
        unless ($start == $TX->{txStart}) {
          my $ds  = $start - $chrpos;
          $dstart = $ds if abs($ds) < abs($dstart);
        }
        # Skip exon end if it is at the transcript end -- no junction
        unless ($end == $TX->{txEnd}) {
          my $de  = $end - $chrpos - 1;
          $dend   = $de if abs($de) < abs($dend);
        }
      }
      # dstart and dend are swapped and reversed for minus strand transcripts
      ($dstart, $dend) = (-$dend, -$dstart) if $TX->{strand} eq '-';
      # 5' donor (>GT), 3' acceptor (AG<)
      # Offset dstart and dend to point to intron end (last intron base) and intron start (first intron base)
      $SNP->{Gene}{JXdon} = $dend   + 1;
      $SNP->{Gene}{JXacc} = $dstart - 1;
    }
  }

  ##### |     find_gene returns: $SNP->{Gene}
  return 1;
}

sub CpG {
  my ($nt5, $allele1, $allele2, $nt3) = @_;
  my $t1  = $nt5 . $allele1 . $nt3;
  my $t2  = $nt5 . $allele2 . $nt3;
  warn("WARNING: Incorrect triplet(s): $t1>$t2\n"), return
    unless length($t1) == 3 && length($t2) == 3;
  my $cpg = 0;
  $cpg    = 1 if index($t1, 'CG') >= 0;
  $cpg   += 2 if index($t2, 'CG') >= 0;
  return $cpg;
}

1;
