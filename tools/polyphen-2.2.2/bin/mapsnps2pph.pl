#!/usr/bin/env perl

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=pod

=head1 NAME

mapsnps2pph.pl - Convert output from MapSNPs annotation pipeline into
                 PolyPhen-2 query input format

=head1 SYNOPSIS

mapsnps2pph.pl [-o <snps.txt>] [-s] file1.snps [file2.snps ...]

where options are:

  -o <filename>   save full annotation report to <filename>;
                  default filename: 'mapsnps.txt'

  -s              substitute knownGene IDs for the missing protein accessions;
                  default is to skip missense SNPs with no protein accessions
                  and do not include them into PolyPhen-2 query output

=head1 DESCRIPTION

Outputs missense SNP annotations in PolyPhen-2 query input format to STDOUT
and writes full annotation report for all SNP functional categories found
in the input file(s) to 'mapsnps.txt' file by default, or to a filename
specified via '-o' option.

Specifying '-s' option assumes that '-s' option was also used when obtaining
MapSNPs annotations and corresponding CDS sequences tagged by UCSC knownGene
transcript IDs were already saved to a file in FASTA format.

=head1 AUTHORS

IAA

=head1 SUBVERSION

 $LastChangedDate: 2010-11-27 20:35:38 -0500 (Sat, 27 Nov 2010) $
 $LastChangedRevision: 338 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use Getopt::Std;
use Pod::Usage;

pod2usage('ERROR: No input file(s) specified') unless @ARGV;
my %opts;
getopts('hso:', \%opts) or pod2usage('ERROR: Illegal argument(s)');
pod2usage(-verbose=>2) if $opts{'h'};
my $subst_acc  = $opts{'s'} || 0;
my $snpfile    = $opts{'o'} || 'mapsnps.txt';

my $pph2format = join("\t", qw{ %-12s %6s %3s %3s %12s %6s %5s %3s %3s %4s %3s %3s %6s %6s });

my %totals;
open(SNPF, ">$snpfile") or die "Can't create output file: $snpfile\n";
foreach my $file (@ARGV) {
  open(FIN, $file) or die "Can't open file: $file\n";
  while (<FIN>) {
    next if /^\s*$/;
    # Aggregate total MapSNPs statistics from all input files
    if (/^##/) {
      next if /^## Totals/;
      my @a = split;
      my $key   = join(' ', @a[1..($#a-1)]);
      my $value = $a[-1];
      next unless $value =~ /^\d+$/;
      $totals{$key} += $value;
      next;
    }
    my @a;
    my $comments;
    if (/^#/) {
      print SNPF;
      next;
    } else {
      chomp;
      # Trim optional user comments
      $comments = $1 if s/\s *#\s*(.*)$//;
      @a = split /\t/, $_, 38;  # preserve any trailing empty fields
      map { s/^\s?(\s*)$/\1?/; } @a;
      print SNPF join("\t", @a);
      print SNPF "\t# $comments" if defined $comments && length $comments;
      print SNPF "\n";
    }
    @a = split /\t/, $_, 38;  # preserve any trailing empty fields
    map { s/^\s+//; s/\s+$//; } @a;
    # Skip all but missenses
    next unless $a[8] eq 'missense';
    # If mapping to UniProtKB was attempted and failed, i.e., spmap (column #29)
    # equals 0 but not empty, then output knownGene transcript name (column #4)
    # as a protein ID (colums #30) instead and cdnpos (column #22) as aapos (column #28),
    # but only if '-s' option was specified.
    if ($a[28] ne '') { # spmap non-empty, mapping was attempted
      unless ($a[28]) { # spmap code is 0, mapping failed
        if ($subst_acc) { $a[29] = $a[3]; $a[27] = $a[21]; }
        else { next; }
      }
    # No UniProtKB mapping was attempted, skip when spacc (column #30) value is absent,
    # unless '-s' option was set in which case substitute columns #30 as above. Also
    # substitute column #28 as above in any case.
    } else {
      $a[27] = $a[21];
      unless ($a[29]) { # spacc is empty, protein accession is not available
        if ($subst_acc) { $a[29] = $a[3]; }
        else { next; }
      }
    }
    $a[7] =~ s|/||;
    my $snpid = "$a[0]|$a[7]|$a[3]$a[1]|$a[2]|$a[31]";
    printf $pph2format,
      $a[29],
      $a[27],
      $a[25],
      $a[26],
      $a[3],
      $a[21],
      $a[22],
      $a[10],
      $a[11],
      $a[12],
      $a[13],
      $a[14],
      $a[15],
      $a[16]
    ;
    my $idcomm = defined $comments && length $comments ? "$snpid\t$comments" : $snpid;
    print "\t# $idcomm\n";
    $ocount++;
  }
  close(FIN);
}
# Print out total aggregated statistics
if (%totals) {
  print  SNPF "## Totals:\n";
  printf SNPF "##   %-18s %8s\n", 'lines input', $totals{'lines input'} || 0;
  printf SNPF "##   %-18s %8s\n", 'lines skipped', $totals{'lines skipped'} || 0;
  printf SNPF "##   %-18s %8s\n", 'alleles annotated', $totals{'alleles annotated'} || 0;
  for my $key (qw{ missense nonsense coding-synon intron utr-3 utr-5 }) {
    printf SNPF "##     %-16s %8s\n", $key, $totals{$key} || 0;
  }
}
close(SNPF);
