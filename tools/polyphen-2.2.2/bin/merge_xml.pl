#!/usr/bin/env perl

use warnings;
use strict;

=head1 NAME

merge_xml.pl - merge results of PolyPhen-2 classifier run into corresponding
               SNP(s) contained in a file in the PPH XML format

=head1 SYNOPSIS

merge_xml.pl pph.xml modelfile prediction pph2_class pph2_prob pph2_FPR pph2_TPR pph2_FDR

--or--

merge_xml.pl -f pph2.tab [-m modelfile] pph.xml

where mandatory argument is:

  -f pph2.tab      PolyPhen-2 report summary file (tab-delimited)

and options are:

  -m modelfile     WEKA model file used by classifier; default
                   is $PPH/models/HumDiv.UniRef100.NBd.f11.model

=head1 SUBVERSION

 $LastChangedDate: 2009-11-28 00:15:09 -0500 (Sat, 28 Nov 2009) $
 $LastChangedRevision: 254 $
 $LastChangedBy: ivan $

=cut

use Getopt::Std;
use Pod::Usage;
use File::Path;
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

my %opts;
pod2usage('No arguments specified!') unless @ARGV;
getopts('f:m:', \%opts) or pod2usage('Wrong arguments!');
if ($opts{'f'}) {
  pod2usage('XML file not specified!') unless @ARGV;
  pod2usage('Multiple XML files not supported!') if @ARGV > 1;
} else {
  pod2usage('Wrong number of arguments!') unless @ARGV == 8;
  pod2usage('Options not allowed in single SNP mode!') if %opts;
}
my $xmlfile = shift @ARGV;
warn "ERROR: File does not exist: $xmlfile\n" unless -e $xmlfile;
warn "ERROR: File empty: $xmlfile\n" unless -s $xmlfile;
open my $FX, '<', $xmlfile or die "ERROR: failed to open XML file: $xmlfile\n";

# Read classifier outcome and scores from a tabulated report file
my $PPH = PPH->new;
if ($opts{'f'}) {
  my $modfile = $opts{'m'} || $CONFIG{'WEKAMODEL'}; # WEKA model filename (not checked!)
  my $modname = fileparse($modfile, qw(.model .mod));
  my $sumfile = $opts{'f'};
  die "File does not exist: $sumfile\n" unless -e $sumfile;
  die "File empty: $sumfile\n" unless -s $sumfile;
  my %summary;
  open(SUMF, $sumfile) or die "Can't open file: $sumfile\n";
  while (<SUMF>) {
    # Skip comments
    next if /^#/;
    # Skip lines with errors & warnings and empty lines
    next if /\tError:\s/;
    next if /^(ERROR|WARNING):\s/;
    next if /^\s*$/;
    chomp;
    my @a = split /\t/;
    map { if (defined) { s/^\s+//; s/\s+$//; } } @a;
    # Verify data
    warn("Invalid input: $_\n"), next
      # column  12: prediction
      # columns 15-19: pph2_class pph2_prob pph2_FPR pph2_TPR pph2_FDR
      unless
        # valid prediction
        (defined $a[14] && ($a[14] eq 'deleterious' || $a[14] eq 'neutral') &&
         defined $a[15] &&  $a[15] =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && $a[15] >= 0 && $a[15] <= 1 &&
         defined $a[16] &&  $a[16] =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && $a[16] >= 0 && $a[16] <= 1 &&
         defined $a[17] &&  $a[17] =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && $a[17] >= 0 && $a[17] <= 1 &&
         defined $a[18] &&  $a[18] =~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/ && $a[18] >= 0 && $a[18] <= 1) ||
        # dummy but properly formatted prediction (aka "unknown")
        (defined $a[11] && $a[11] eq 'unknown' && defined $a[14] && $a[14] eq 'none')
      ;
    # Unique identifier for the variant: acc : pos aa1 aa2
    my $snpid = "$a[5]:$a[6]$a[7]$a[8]";
    push @{ $summary{$snpid} }, @a[11, 14..18];
  }
  close(SUMF);
  $PPH = $PPH->xml2pph($FX);
  close($FX);
  foreach my $SNP (@{ $PPH->{Variant} }) {
    my $snpid = $SNP->{Acc} . ':' . $SNP->{Pos} . $SNP->{Aa1} . $SNP->{Aa2};
    if (exists $summary{$snpid}) {
      my @a = @{ $summary{$snpid} };
      push @{ $SNP->{Prediction} }, {
        Method  => 'NBd',
        Outcome => $a[0],
        Model   => $modname,
        Class   => $a[1],
        Prob    => $a[2],
        FPR     => $a[3],
        TPR     => $a[4],
        FDR     => $a[5]
      };
    } else {
      warn "WARNING: Unmatched SNP ($snpid) skipped in file: $xmlfile\n";
    }
  }
  open $FX, '>', $xmlfile or die "ERROR: failed to create XML file: $xmlfile\n";
  $PPH->pph2xml($FX);
# Read classifier outcome and scores from command line and add them to
# the first variant in XML file, without checking for SNP identifiers.
# Note: This mode is obsolete, do not use.
} else {
  my $modfile = shift; # WEKA model filename (not checked!)
  my $modname = fileparse($modfile, qw(.model .mod));
  my ($prediction, $class, $prob, $fpr, $tpr, $fdr) = @ARGV;
  $PPH = $PPH->xml2pph($FX);
  close($FX);
  my $SNP = $PPH->{Variant}[0];
  # No data verification!
  push @{ $SNP->{Prediction} },
    { Method=>'NBd', Outcome=>$prediction, Model=>$modname, Class=>$class, Prob=>$prob, FPR=>$fpr, TPR=>$tpr, FDR=>$fdr };
  open $FX, '>', $xmlfile or die "ERROR: failed to create XML file: $xmlfile\n";
  $PPH->pph2xml($FX);
}
