#!/usr/bin/env perl
use strict;

while (<>) {
  if (/^>/) {
    if (/^>UniRef(\d+)_(\S+)(\s+.+)/) {
      $_ = ">gnl|uniref$1|$2$3\n";
    } else {
      die "Unsupported defline format:\n$_";
    }
  } else {
    # Replace Selenocysteine (U) with Cysteine (C) and Pyrrolysine (O) with
    # Lysine (K) in the sequence since MAFFT as well as many other tools
    # do not recognize these non-standard residue codes.
    tr/UOuo/CKck/;
  }
  print;
}
