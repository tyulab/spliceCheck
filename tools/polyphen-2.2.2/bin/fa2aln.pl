#!/usr/bin/env perl
use warnings;
use strict;

=head1 NAME

fa2aln.pl - convert FASTA multiple alignment file(s) to PolyPhen-2 format

=head1 SYNOPSIS

fa2aln.pl [options] <msa1.fa> [msa2.fa]...

where options are:

  -p        name converted .aln file using query's protein accession

  -i        include top (query) sequence in the converted .aln file

=head1 SUBVERSION

 $LastChangedDate: 2009-11-28 00:15:09 -0500 (Sat, 28 Nov 2009) $
 $LastChangedRevision: 254 $
 $LastChangedBy: ivan $

=cut

use Getopt::Std;
use Pod::Usage;
use File::Basename;

BEGIN {
  use FindBin qw($RealBin);
  BEGIN {
    unless ($ENV{'PPH'}) {
      my $binpath = $RealBin;
      $binpath =~ s|/[^/]+$||;
      $ENV{'PPH'} = $binpath;
    }
  }
  use lib "$ENV{PPH}/perl";
}

use PPH::Config;
use PPH;
use SNP;
use PPH::Align;

my %opts;
pod2usage('No arguments specified!') unless @ARGV;
getopts('ip', \%opts) or pod2usage('Wrong arguments!');
pod2usage('No input file specified!') unless @ARGV;

my $SNP = SNP->new_empty();

foreach my $fafile (@ARGV) {

  my ($name, $path, $ext)   = fileparse($fafile, qw{ .msa.fa .fasta .fas .fa });
  my $aln_file = "$name.aln";
  die "ERROR: Same input and output file: $aln_file\n" if $aln_file eq $fafile;

  # Load FASTA MAF into memory
  my ($maf, $acc) = PPH::Align::_read_fasta_maf($fafile);

  # Convert alignments to PSIC/CLUSTAL format
  my $out = '';
  my $query_pid = $acc->[0];
  my $gapped_q_seq   = $maf->{$query_pid}{seq};
  my $ungapped_q_seq = $gapped_q_seq;
  $ungapped_q_seq =~ tr/-//d; # This is a hack!
  die "ERROR: Empty or missing QUERY sequence in MAF: $fafile" unless length $gapped_q_seq;
  my $nseq = scalar @$acc;
  $nseq-- unless $opts{'i'};
  foreach my $accession (@{ $acc }) {
    my $alignment = $maf->{$accession};
    # Skip QUERY
    next if !$opts{'i'} && $accession eq $acc->[0];
    # Ungap aligned sequence
    my (undef, $seq) = PPH::Align::_ungap_HSP_alignment($ungapped_q_seq, $gapped_q_seq, $alignment->{seq});
    # Output PSIC/CLUSTAL-formatted description and sequence
    $out .= PPH::Align::_formatDesc4PSIC_aln($alignment->{tag}, $alignment->{desc}) . $seq . "\n";
  }

  $aln_file = "$query_pid.aln" if $opts{'p'};

  # Create output .aln file
  open my $C, '>', $aln_file or die "ERROR: Can't create output file: $aln_file\n";

  # Print CLUSTAL file header
  print $C "CLUSTAL $query_pid ($nseq)\n\n";
  print $C $out;

  close($C);
}
