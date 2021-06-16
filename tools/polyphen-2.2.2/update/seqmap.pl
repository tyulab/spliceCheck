#!/usr/bin/env perl

use warnings;
use strict;

# Maps one set of protein sequences (e.g., UniProtKB) to another set
# of protein sequences (e.g. UCSC knownGene transcripts).
#
# Usage:
#   seqmap.pl [options] <query.seq> <target.seq>
#
# Both sequence files should be FASTA-formatted with each sequence
# listed on one line, unwrapped.
#

use Getopt::Std;
use Pod::Usage;

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
use PPH::Align;

my $usage = <<'EOT';
Usage: seqmap.pl [options] <seq1.fa> <seq2.fa>
  where options are:
    -i X   sequence identity threshold (default: 0.97)
           specify X=0 to disable identity filtering
    -o Y   fraction of overlap threshold (default: 0.20)
           specify Y=0 to disable overlap filtering
    -h K   maximum number of top scoring BLAT hits per
           each query to retain (default: K=0, include
           all hits that have passed identity and
           overlap filters)
    -r N/M partition input into M equally sized stripes
           and process stripe N
    -v     verbose output
EOT
$usage =~ s/\n$//s;

# Note: Unlike UCSC/blat positional data, all coordinates in output are zero-based, closed end [start,end]
my $header = join("\t", 'qacc      ', 'macc', 'mname      ',
  qw{ qlen mlen overlap hits qstart qend mstart mend qfrac ident map qcan mcan bsizes qstarts mstarts });

my %opts;
getopts('vni:h:o:r:', \%opts) or pod2usage($usage);
my $ident_level     = exists $opts{'i'} ? $opts{'i'} : 0.97;
my $overlap_level   = exists $opts{'o'} ? $opts{'o'} : 0.20;
my $no_of_hsps      = exists $opts{'h'} ? $opts{'h'} : 0;
my $verbose         = $opts{'v'} || 0;
die "$usage\n" unless @ARGV == 2;
my ($query_file, $target_file) = @ARGV;

die "ERROR: Target file not found: $target_file\n" unless -e $target_file;

# Check if we are running under SGE or LSF and if this is an array job
my ($grid_mode, $jobid, $taskid, $maxtask, $grid_msg);
if (exists $ENV{'SGE_ROOT'} && exists $ENV{'JOB_ID'}) {       # SGE
  $grid_mode = 'SGE';
  $jobid  = $ENV{'JOB_ID'} ne 'undefined' ? $ENV{'JOB_ID'} : 0;
  $taskid = exists $ENV{'SGE_TASK_ID'} && $ENV{'SGE_TASK_ID'} ne 'undefined' ? $ENV{'SGE_TASK_ID'} : 0;
} elsif ($ENV{'LSB_QUEUE'} && exists $ENV{'LSB_JOBID'}) {     # LSF
  $grid_mode = 'LSF';
  $jobid  = $ENV{'LSB_JOBID'};
  $taskid = exists $ENV{'LSB_JOBINDEX'} ? $ENV{'LSB_JOBINDEX'} : 0;
}
my $array_mode = $jobid && $taskid ? 1 : 0;

if ($grid_mode) {
  # Adjust PATH for a grid job
  if ($CONFIG{'GRID_PATH'}) {
    my $ARCH = `uname -m`; chomp $ARCH;
    $ENV{'PATH'} = eval "qq{$CONFIG{GRID_PATH}}";
  }
  $grid_msg = "I am $grid_mode job $jobid \@ $ENV{HOSTNAME}\n";
}

if ($array_mode) {
  # GGI_STRIPES overwrites any array size specified initially
  $maxtask = $ENV{'GGI_STRIPES'} || $ENV{'SGE_TASK_LAST'} || $ENV{'LSB_JOBINDEX_END'};
  $grid_msg = "I am $grid_mode job $jobid task $taskid of $maxtask \@ $ENV{HOSTNAME}\n";
}

my ($stripe_no, $stripes);
if ($opts{'r'}) {
  ($stripe_no, $stripes) = $opts{'r'} =~ m|^(\d+)/(\d+)$|;
  die "ERROR: Illegal option's argument: -r $opts{'r'}\n" unless $stripe_no && $stripes;
}

# Debugging
warn $grid_msg if defined $grid_msg;

if ($array_mode) {
  warn "WARNING: Explicit input stripe specification via '-r' option ignored in $grid_mode array mode\n" if $opts{'r'};
  ($stripe_no, $stripes) = ($taskid, $maxtask);
}

# Count how many query sequences we have in input
my $input_query_count = 0;
my ($range_start, $range_end);
if ($stripes) {
  open(FQUERY, $query_file) or die "ERROR: Can't open query file: $query_file\n";
  while (<FQUERY>) {
    $input_query_count++ if m|^>\S+\s|;
  }
  close(FQUERY);
  # Calculate stripe range
  ($range_start, $range_end) = stripe_range($stripe_no, $stripes, $input_query_count);
}

my $matched_count   = 0;
my $unmatched_count = 0;
my $skipped_count   = 0;
my $totals_count    = 0;
my $query_count     = 0;
my $lineno          = 0;
my %map;
open(FQUERY, $query_file) or die "ERROR: Can't open input file: $query_file\n";
print "#$header\n";
while (<FQUERY>) {

  $lineno++;
  chomp;
  next if /^\s*$/;
  $totals_count++;
  if ($stripes) {
    <FQUERY>, next if $totals_count < $range_start;
    <FQUERY>, last if $totals_count > $range_end;
  }
  $query_count++;
  my $qdef = $_;
  die "ERROR: Incorrect FASTA definition format at line $lineno of input file: $query_file\n"
    unless $qdef =~ m|^>\S+\s|;
  my $qseq = <FQUERY>;
  chomp $qseq;
  die "ERROR: Sequence missing at line $lineno of input file: $query_file\n"
    unless defined $qseq && length $qseq;
  my $qfa  = $qdef . "\n" . $qseq . "\n";
  my ($qtag, $qdb, $qacc, $qname, $qdesc) = PPH::Align::_parse_defline($qdef);
  my $qcan = $qdesc =~ /^Canonical;/ ? 1 : 0;
  my $qlen = length $qseq;
  my ($macc, $mid, $mlen, $overlap, $hits, $qstart, $qend, $mstart, $mend, $qfrac, $mapped, $mcan, $bsizes, $qstarts, $mstarts);
  my $ident = 0;

  # Try grep mapping first (returns a single first hit found)
  my ($mdesc, $mseq);
  ($macc, $mid, $mdesc, $mseq) = grep_in_FASTA($qseq, $target_file);

  my @hsps;
  if (defined $macc) {

    $mcan    = $mdesc =~ /^Canonical;/ ? 1 : 0;
    $mlen    = length $mseq;
    $overlap = $qlen;     # always full size query for grep matches
    $hits    = 1;         # always single ungapped hit for grep matches
    $qstart  = 0;         # always 0 for grep matches
    $qend    = $qlen - 1; # zero base, always full size query for grep matches
    my $qoffset = $qlen == $mlen ? 0 : index($mseq, $qseq);
    $mstart  = $qstart + $qoffset;
    $mend    = $qend   + $qoffset;
    $qfrac   = sprintf "%.3g", $qlen / $mlen;
    $ident   = 1;         # always 1.00 for grep matches
    $mapped  = 'g';
    $bsizes  = $mlen;
    $qstarts = $qstart;
    $mstarts = $mstart;
    push @hsps, [ $qacc, $macc, $mid, $qlen, $mlen, $overlap, $hits, $qstart, $qend, $mstart, $mend, $qfrac, $ident, $mapped, $qcan, $mcan,
                  $bsizes, $qstarts, $mstarts ];

  # Try BLAT if grep failed
  } else {

    foreach my $hsp (map_by_BLAT($qfa, $target_file)) {
      ($macc, $mid, $mlen, $ident, $overlap, $hits, $qstart, $qend, $mstart, $mend, $bsizes, $qstarts, $mstarts) = @$hsp;
      if (defined $macc) {
        $ident  = sprintf("%.3g", $ident);
        $qfrac  = sprintf "%.3g", $overlap / $qlen;  # fraction of query sequence aligned
        $mapped = 'b';
        $mdesc  = (grep_FASTA_defline($macc, $target_file))[2];
        $mcan   = $mdesc =~ /^Canonical;/ ? 1 : 0;
        for ($bsizes, $qstarts, $mstarts) { s/,$//; }
        push @hsps, [ $qacc, $macc, $mid, $qlen, $mlen, $overlap, $hits, $qstart, $qend, $mstart, $mend, $qfrac, $ident, $mapped, $qcan, $mcan,
                      $bsizes, $qstarts, $mstarts ];
      }
    }

  }

  my $hsp_count = 0;
  if (defined $mapped) {
    foreach my $hsp (
        # Sort on:
        #   canonical target hits first, mcan (desc)
        #   fraction of query sequence in overlap (desc),
        #   number of alignment blocks (asc)
        #   fraction of identity in overlap (desc),
        sort { $b->[15] <=> $a->[15] || $b->[11] <=> $a->[11] || $a->[6] <=> $b->[6] || $b->[12] <=> $a->[12] } @hsps
      ) {
      last if $no_of_hsps && $hsp_count >= $no_of_hsps;
      if      ($hsp->[12] < $ident_level) {    # ident
        warn "WARNING: ($hsp->[0] > $hsp->[1]) Identity ($hsp->[12]) below threshold ($ident_level), skipped\n" if $verbose;
        $skipped_count++;
      } elsif ($hsp->[11] < $overlap_level) {   # overlap
        warn "WARNING: ($hsp->[0] > $hsp->[1]) Overlap ($hsp->[11]) below threshold ($overlap_level), skipped\n" if $verbose;
        $skipped_count++;
      } else {
        print join("\t", @$hsp), "\n";
        $hsp_count++;
        $matched_count++;
      }
    }
  } else {
    $unmatched_count++;
  }
  #last if $totals_count >= 100;

}
close(FQUERY);

printf "%-16s%6d\n%-16s%6d\n%-16s%6d\n%-16s%6d\n",
  '## Input:',        $query_count,
  '## Unmapped:',     $unmatched_count,
  '## Skipped:',      $skipped_count,
  '## Mapped:',       $matched_count;

warn "TOTAL: $query_count queries processed ($range_start-$range_end/$input_query_count)\n" if $array_mode;

#------------------------------------------------------------------------------

=head2 grep_in_FASTA

 Usage   : my ($acc, $name, $desc, $seq) = grep_in_FASTA($query_seq, $target_file)
 Function: Find a (sub) sequence match in a FASTA formatted target file
           File has to keep each sequence on ONE LINE!
 Args    : query sequence
           target file with FASTA sequences
 Returns : Returns FIRST hit which contains the match to query sequence
           Accession
           Entry name
           Description
           Sequence

=cut

#----------------------------------------
sub grep_in_FASTA {
  my ($query_seq, $target_file) = @_;

  my $pid = open my $F, '-|', "fgrep -s -m 1 -B 1 '$query_seq' $target_file";
  die "ERROR: grep: $!\n" unless defined $pid;
  my $head = <$F>;
  close($F), return () unless defined $head;
  my $seq = <$F>;
  close $F or die "\n";

  chomp $head;
  chomp $seq;

  ($head =~ s/^>//) or die "ERROR: Incorrect format of FASTA defline: $head\n";
  my ($tag, $db, $acc, $name, $desc) = PPH::Align::_parse_defline($head);

  return ($acc, $name, $desc, $seq);
}

sub grep_FASTA_defline {
  my ($qacc, $target_file) = @_;

  my $pid = open my $F, '-|', "fgrep -s -m 1 '|$qacc' $target_file";
  die "ERROR: grep: $!\n" unless defined $pid;
  my $head = <$F>;
  close($F), return () unless defined $head;
  close $F or die "\n";
  chomp $head;
  ($head =~ s/^>//) or die "ERROR: Incorrect format of FASTA defline: $head\n";
  my ($tag, $db, $acc, $name, $desc) = PPH::Align::_parse_defline($head);

  return ($acc, $name, $desc);
}

sub map_by_BLAT {
  my ($query, $db) = @_;
  my $pident  = $ident_level * 100;  # convert to percent
  my $queryfa = File::Temp->new();
  print $queryfa $query;
  $queryfa->flush;
  # usage: blat [options] database query output.psl
  my $cmd = qq(blat -prot -minIdentity=$pident -out=psl -noHead $db $queryfa stdout);
  my @rc = `$cmd`;
  my @hsps;
  for (@rc) {
    chomp;
    next unless length;
    #     0        1       2   3     4        5     6        7      8     9    10     11   12    13    14     15   16     17     18      19      20
    # match mismatch repmatch ns qgaps qgapsize mgaps mgapsize strand qname qsize qstart qend mname msize mstart mend bcount bsizes qstarts tstarts
    my @a = split /\t/;
    my ($mtag, $mdb, $macc, $mname, $mdesc) = PPH::Align::_parse_defline($a[13]);
    my $match = $a[0] + $a[2];  # count all matches
    my $overlap = 0; map { $overlap += $_; } split /,/, $a[18];
    # macc mname mlen ident length gaps qstart qend mstart mend bsizes qstarts mstarts (13)
    push @hsps, [
      $macc,
      $mname,
      $a[14],
      $match / $overlap,
      $overlap,
      $a[17],
      $a[11],
      $a[12] - 1, # qend, convert to zero-based closed end
      $a[15],
      $a[16] - 1, # mend, convert to zero-based closed end
      $a[18],
      $a[19],
      $a[20]
    ];
  }
  return @hsps;
}

# Calculate [start,end] range of a stripe in dataset
# Note: all indices passed and returned are one-based
#
# Arguments:  stripe index, number of stripes, dataset size
# Returns:    data start index, data end index; both indices
#             returned are zero if no more data fall within
#             a stripe specified (e.g. number of stripes is
#             larger than the number of data items)
#
sub stripe_range {
  my ($stripe_no, $stripes, $size) = @_;

  die "ERROR: Stripe index too large: $stripe_no > $stripes\n" if $stripe_no > $stripes;

  my $p = int($size / $stripes);
  my $r = $size % $stripes;

  my $shift_start = ($stripe_no <= $r ? $stripe_no : $r + 1) - 1;
  my $shift_end   = ($stripe_no <= $r ? $stripe_no : $r);

  my $start = ($stripe_no - 1)*$p + 1 + $shift_start;
  my $end   = $stripe_no*$p + $shift_end;

  if ($start > $end) {
    return (0, 0);
  } else {
    return ($start, $end);
  }
}
