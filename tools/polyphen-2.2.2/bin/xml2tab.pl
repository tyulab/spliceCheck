#!/usr/bin/env perl

use warnings;
use strict;

use Getopt::Std;
use Pod::Usage;
use XML::Simple;

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

use PPH;
use SNP;
$| = 1;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

xml2tab.pl - print contents of PPH XML file(s) as a tab-separetd table

=head1 SYNOPSIS

xml2tab.pl xmlfile1 [xmlfile2 ...]

=head1 DESCRIPTION


=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2011-09-20 21:30:07 -0400 (Tue, 20 Sep 2011) $
 $LastChangedRevision: 368 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

@ARGV or pod2usage();

my $PPH = PPH->new();

print '#', join("\t", $PPH->return_header(-formatted=>1)), "\n";
foreach my $xmlf (@ARGV) {
  -e $xmlf or next;
  open my $F, '<', $xmlf or die "ERROR: Can't open file: $xmlf\n";
  $PPH = $PPH->xml2pph($F);
  close($F);
  foreach my $SNP (@{ $PPH->{Variant} }) {
    print join("\t", $PPH->format_result($SNP)), "\n";
  }
}
