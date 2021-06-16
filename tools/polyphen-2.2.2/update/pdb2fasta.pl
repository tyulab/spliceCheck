#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Copy;
use DBI;

use FindBin qw($RealBin);
BEGIN {
  unless ($ENV{'PPH'}) {
    my $binpath = $RealBin;
    $binpath =~ s|/[^/]+$||;
    $ENV{'PPH'} = $binpath;
  }
}

#------------------------------------------------------------------------------------
my %aa3to1 = (
  'ALA' => 'A',   'CYS' => 'C',   'ASP' => 'D',   'GLU' => 'E',
  'PHE' => 'F',   'GLY' => 'G',   'HIS' => 'H',   'ILE' => 'I',
  'LYS' => 'K',   'LEU' => 'L',   'MET' => 'M',   'ASN' => 'N',
  'PRO' => 'P',   'GLN' => 'Q',   'ARG' => 'R',   'SER' => 'S',
  'THR' => 'T',   'VAL' => 'V',   'TRP' => 'W',   'TYR' => 'Y',
  'UNK' => 'X',
);

my $ONE_CHAIN = '_';
my $MIN_SEQ_LEN = 10;
my $UNKNOWN_RES = 'X';

#------------------------------------------------------------------------------------

my $USAGE = <<'EOT';
Usage: pdb2fasta.pl [-l|wwpdb_path] [>logfile]
EOT

my %opts;
getopts('l', \%opts) or die $USAGE;
die $USAGE unless $opts{'l'} || @ARGV == 1;
my $path = shift @ARGV;

# no PQS anymore
my $dbase             = 'pdb';
my $seqfile           = "$dbase.seq";
my $fragmfile         = "$dbase.fragm";
my $formatdb_logfile  = "$dbase.makeblastdb.log";
my $dbfile            = "$dbase.sqlite";

$| = 1;

# Load precompiled fragm data into SQLite database and quit
createdb($fragmfile, $dbfile), exit(0) if $opts{'l'};

print "Scanning PDB directory ...\n";
opendir(PDBDIR, $path) || die "Can't readdir $path: $!\n";
my @files = sort grep { /\.(pdb|ent|brk)(\.gz|\.Z)$/i } readdir(PDBDIR);
closedir PDBDIR;
my $pdbcount = scalar(@files) || 0;
print "Found $pdbcount PDB files to process.\n\n";
#--------------------------------------------------------------

open(FRAGM, ">$fragmfile") || die "Can't write to $fragmfile: $!\n";
open(SEQ,   ">$seqfile")   || die "Can't write to $seqfile: $!\n";

my $count = 0;
foreach my $file (@files) {

  $file=~/(\d[0-9A-Z]{3})\./i;
  my $id = $1;

  my $pdbfile = "$path/$file";
  my $r = &read_pdb_atoms($pdbfile);
  my $desc = &get_pdb_desc($pdbfile);

 CHAIN: foreach my $chain (sort keys %{$r}) {

    my $c = $r->{$chain};

    print "$id $chain\t";

    my $seq = &extract_seq_from_chain($c);

    if ($seq=~/[acgtu]/) {
      print "DNA, rejected\n";
      next CHAIN;
    }				# end if

    my $len = length $seq;

    if ($len<$MIN_SEQ_LEN) {
      print "short, rejected\n";
      next CHAIN;
    }				# end if

    if ($seq=~/^$UNKNOWN_RES+$/o) {
      print "all unknown residues, rejected\n";
      next CHAIN;
    }				# end if

    my $rnums = [ map { $_->{Num} } @{$c} ];
    my $rfragms = &extract_fragments($rnums);
    my $f = $#{$rfragms}+1;	## number of fragments

    foreach my $x (@{$rfragms}) {
      printf FRAGM "%-6s %s  %5d %5d  %5s %5s\n",
        $id, $chain,
        $x->{BegIdx}+1, $x->{EndIdx}+1, ## +1 since seq numbering is 1-based
        $x->{Beg}, $x->{End}
      ;
    }				# end foreach

    print SEQ &format_fasta_seq("$dbase|$id|$chain $desc", $seq);

    print "OK\t", $len, " residues in ", scalar @{$rfragms}, " fragments\n";

    undef $rfragms;
    undef $r->{$chain};

  }				# end foreach chain

  undef $r;
  $count++;

}				# end foreach pdb

close SEQ;
close FRAGM;
close LOG;

my $MAKEBLASTDB_NAME = 'makeblastdb';
my $MAKEBLASTDB_ARGS = "-in $seqfile -dbtype prot -out $dbase -logfile $formatdb_logfile";
my $MAKEBLASTDB_PATH;
print "\nRunning $MAKEBLASTDB_NAME... ";
# Attempt to locate makeblastdb executable
my $mkmsg = "Please run the following command manually:\n\t$MAKEBLASTDB_NAME $MAKEBLASTDB_ARGS\n";
if      (-x "$ENV{PPH}/blast/bin/$MAKEBLASTDB_NAME") {
  $MAKEBLASTDB_PATH = "$ENV{PPH}/blast/bin/$MAKEBLASTDB_NAME";
} elsif (-x "$ENV{PPH}/bin/$MAKEBLASTDB_NAME") {
  $MAKEBLASTDB_PATH = "$ENV{PPH}/bin/$MAKEBLASTDB_NAME";
} else {
  my $rc = `which $MAKEBLASTDB_NAME 2>&1`; chomp $rc;
  $MAKEBLASTDB_PATH = $rc unless $rc =~ /no $MAKEBLASTDB_NAME/o;
}
die "failed to locate $MAKEBLASTDB_NAME executable!\n$mkmsg" unless $MAKEBLASTDB_PATH;

my $cmd = "$MAKEBLASTDB_PATH $MAKEBLASTDB_ARGS";
my $rc  = `$cmd 2>&1`;
if ($?) {
  die "failed to execute makeblastdb: $!\n$mkmsg";
} elsif ($rc) {
  die "failed: $rc\n$mkmsg";
} else {
  local $/ = undef;
  open(LOG, "$formatdb_logfile") or
    die "failed: Can't open $formatdb_logfile\n$mkmsg";
  my $log = <LOG>;
  close(LOG);
  $log =~ s/^\s+//s;
  if ($log =~ /^Error:/m) {
    die "failed:\n$log$mkmsg";
  } else {
    print "done.\n";
  }
}

# Load compiled fragm file into SQLite database
createdb($fragmfile, $dbfile);
print "All done.\n";

#---------------------------------------------------------------------------------------

# Create and populate SQLite database
#
#   fstart / fend correspond to sequence ranges or single positions in FASTA sequence
#   astart / aend correspond to ranges or single residue codes in matching PDB ATOM records
#   scale is 1-based for ranges / positions
#
#  pdbid  chain   fstart    fend   astart    aend
#   102m      _        1     154        0     153
#   103l      _        1      34        1      34
#   103l      _       35      35      40A     40A
#   103l      _       36      36      40B     40B
#   103l      _       37      37      40C     40C
#   103l      _       38     159       41     162
#
sub createdb {
  my ($datafile, $dbfile) = @_;
  die "ERROR: createdb: Missing data file name\n" unless length $datafile;
  die "ERROR: createdb: Missing database file name\n" unless length $dbfile;
  die "ERROR: createdb: Missing or empty data file: $datafile\n" unless -s $datafile;
  my $table = q{
  CREATE TABLE fragm (
    id     INTEGER PRIMARY KEY AUTOINCREMENT,
    pdbid  CHAR(4) NOT NULL,
    chain  CHAR(1),
    fstart INTEGER NOT NULL,
    fend   INTEGER NOT NULL,
    astart VARCHAR(6) NOT NULL,
    aend   VARCHAR(6) NOT NULL
  )
  };
  my $insert = q{INSERT INTO fragm VALUES (NULL, } . join(', ', ('?')x6) . ')';
  my @indices = (
    q{CREATE INDEX x_id     ON fragm (pdbid)},
    q{CREATE INDEX x_idch   ON fragm (pdbid, chain)},
    q{CREATE INDEX x_idchst ON fragm (pdbid, chain, fstart)},
    q{CREATE INDEX x_idchen ON fragm (pdbid, chain, fend)}
  );
  print 'Creating SQLite database ...';
  unlink $dbfile or die "ERROR: Failed to remove old copy of SQLite file: $dbfile\n" if -e $dbfile;
  my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
  my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs);
  # Set some pragmas
  $dbh->do('PRAGMA legacy_file_format = OFF');
  $dbh->do('PRAGMA foreign_keys = ON');
  # Create table
  $dbh->do($table);
  print " done.\n";
  # Prepare insert statement
  my $sth = $dbh->prepare($insert);
  # Polupate table
  print 'Populating database    ';
  open(FIN, $datafile) or die "\nERROR: Failed to open file: $datafile\n";
  my $batch_size = 5000;
  my $count      = 0;
  while (<FIN>) {
    next if /^#/ || /^\s*$/;
    chomp;
    unless ($count % $batch_size) {
      if ($count) { print '.'; $dbh->commit; }
      $dbh->begin_work;
    }
    my @r = split;
    eval { $sth->execute(@r); };
    if ($@) {
      warn $@;
      $dbh->rollback;
      print "\n" if $count;
      exit 1;
    }
    $count++
  }
  close(FIN);
  print "\n";
  $dbh->commit;
  # Index table
  print 'Indexing database      ';
  foreach my $index (@indices) {
    print '.';
    $dbh->do($index);
  }
  print "\n";
  print 'Optimizing database    ...';
  $dbh->do('ANALYZE fragm');
  print "\n";
  $dbh->disconnect;
  print "Finished.\n";
}

#----------------------------------------------------------------------------------------------------------------
# 5 6 10 11 11A 12 13 13A 13B 20 21 -> 0-1, 2-3, 4-4, 5-6, 7-7, 8-8, 9-10
sub extract_fragments
  {
    my $r = shift;

    my(@ret);

    $ret[0]{BegIdx}=0;
    $ret[0]{Beg}=$r->[0];

    my $j=0;

    for ( my $i=1; $i<=$#{$r}; $i++) {
      if ( $r->[$i-1]!~/^-?\d+$/ || $r->[$i]!~/^-?\d+$/ || $r->[$i]-$r->[$i-1]!=1 ) {

	$ret[$j]{EndIdx} = $i-1;
	$ret[$j]{End} = $r->[$i-1];

	$j++;

	$ret[$j]{BegIdx} = $i;
	$ret[$j]{Beg} = $r->[$i];


      }				# end if
    }				# end for

    $ret[$j]{EndIdx} = $#{$r};
    $ret[$j]{End} = $r->[$#{$r}];

    return [@ret];

  }				# end sub
#----------------------------------------------------------------------------------------------------------------
sub format_fasta_seq
  {
    my($desc, $seq) = @_;

    my $FASTA_WIDTH = 60;

    my $out = ">$desc\n";

    for (my $i=0; $i<length($seq); $i+=$FASTA_WIDTH ) {
      $out .= substr($seq,$i,$FASTA_WIDTH);
      $out .= "\n";
    }				# end for

    return $out;

  }				# end sub
#--------------------------------------------------------------------------
sub read_pdb_atoms
  {

    my $pdbfile = shift;
    if ($pdbfile =~ /(\.Z|\.gz)$/) {
      $pdbfile = "gunzip < $pdbfile |";
    }

    my $REMOVE_ACE_NH2 = 0;
    my $VERBOSE = 0;

    local($_);

    my(%aaress, @hetress);



    unless(open (PDB, $pdbfile)) {
      print "can't open $pdbfile";
      return undef;
    }				# end unless

    my($curres, $chain);
    my $ter = 0;

    while (<PDB>) {

      if ( /^ATOM/ || (/^HETATM/ && !/HOH/) ) {

        my $a = parse_atom($_);

        if ( defined $curres && ($a->{ResNum} ne $curres->{Num}) ) { # new res
          push @{ $aaress{$chain} }, $curres;
          $curres = undef;
        }			# end if new res

        push @{ $curres->{Atoms} }, $a;

        $curres->{Num}   = $a->{ResNum};
        $curres->{Name}  = $a->{ResName};
        $curres->{Chain} = $chain = $a->{Chain};
        $curres->{Het}   = /^HETATM/ ? 1 : 0;

      }				# end if ATOM
      elsif ( /^TER/ ) {
        $curres->{Ter} = 1;
        $ter = 1;
      }				# end elsif
      elsif ( /^ENDMDL/ ) {
        last;
      }				# end elsif

    }				# end while

    close PDB;

    # process last res:
    push @{ $aaress{$chain} }, $curres;

    #warn "\nNo TER record in $pdbfile" if($VERBOSE && !$ter);

    ####### Additional processing:
    #######
    foreach $chain (keys %aaress) {

      my $i;
      my $r = $aaress{$chain};

      if ($REMOVE_ACE_NH2 && $r->[0]{Name} eq 'ACE') {
        shift @{$r};
        print "\nLeading ACE removed from  $chain chain in $pdbfile" if $VERBOSE;
      }

      if ($REMOVE_ACE_NH2 && $#{$r}!=-1 && $r->[$#{$r}]{Name} eq 'NH2') {
        pop @{$r};
        print "\nTrailing NH2 removed from  $chain chain in $pdbfile" if $VERBOSE;
      }


      for ($i=$#{$r}; $i>=0; $i--) {
        last if( $r->[$i]{Het}==0 || (defined $r->[$i]{Ter} && $r->[$i]{Ter}==1) );
        ###   push @hetress, $r->[$i];
      }

      if ( $i+1 <= $#{$r} ) {
        splice @{ $r }, $i+1;
        print "\nAmbiguous terminal residue removed from $chain chain in $pdbfile" if $VERBOSE && !$ter;
      }				# end if

      delete $aaress{$chain} unless @{$aaress{$chain}};

    }				# end foreach

    ### return { Atoms => \%aaress, Hetatoms => \@hetress };

    return \%aaress;

  }				# end sub
#----------------------------------------------------------------------------------------------------------------
sub extract_seq_from_chain
  {

    my($c) = @_;		# ref to array with chain

    my($seq);
    foreach my $res (@{$c}) {
      $seq .= ($aa3to1{$res->{Name}} || $UNKNOWN_RES);
    }				# end foreach

    return $seq;

  }				# end sub
#----------------------------------------------------------------------------------------------------------------
sub read_seqres
  {

    my($pdbfile) = @_;

    if ($pdbfile =~ /(\.Z|\.gz)$/) {
      $pdbfile = "gunzip < $pdbfile |";
    }

    local($_);

    my(%seqres, $chain, $aa, $gotit);

    (print "Can't open $pdbfile"), return undef unless open (PDB, $pdbfile);

    while (<PDB>) {

      if ( /^SEQRES/ ) {

        $gotit = 1;

        $chain = substr($_,11,1);
        $chain = $ONE_CHAIN if $chain eq ' ';

        for (my $i=19; $i<=67; $i+=4) {

          $aa = substr($_,$i,3);

          $aa =~ s/\s+//g;

          $seqres{$chain} .= ($aa3to1{$aa} || $UNKNOWN_RES ) if $aa =~/\S/;

        }			# end for

      }				# end if SEQRES
      elsif ( $gotit ) {
        last;
      }

    }				# end while <PDB>

    close PDB;

    return { %seqres };

  }				# end sub
#----------------------------------------------------------------------------------------------------------------
sub parse_atom
  {

    #1         2         3         4         5         6         7         8
    #12345678901234567890123456789012345678901234567890123456789012345678901234567890
    #ATOM    145  N   VAL A  25      32.433  16.336  57.540  1.00 11.92      A1   N
    #ATOM    146  CA  VAL A  25      31.132  16.439  58.160  1.00 11.85      A1   C
    #ATOM    147  C   VAL A  25      30.447  15.105  58.363  1.00 12.34      A1   C
    #ATOM    148  O   VAL A  25      29.520  15.059  59.174  1.00 15.65      A1   O
    #ATOM    149  CB AVAL A  25      30.385  17.437  57.230  0.28 13.88      A1   C
    #ATOM    150  CB BVAL A  25      30.166  17.399  57.373  0.72 15.41      A1   C
    #ATOM    151  CG1AVAL    25      28.870  17.401  57.336  0.28 12.64      A1   C
    #ATOM    152  CG1BVAL A  25      30.805  18.788  57.449  0.72 15.11      A1   C
    #ATOM    153  CG2AVAL A  25      30.835  18.826  57.661  0.28 13.58      A1   C
    #ATOM    154  CG2BVAL A  25      29.909  16.996  55.922  0.72 13.25      A1   C

    # from http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
    #COLUMNS        DATA TYPE       FIELD         DEFINITION
    #---------------------------------------------------------------------------------
    # 1 -  6        Record name     "ATOM  "
    # 7 - 11        Integer         serial        Atom serial number.
    #13 - 16        Atom            name          Atom name.
    #17             Character       altLoc        Alternate location indicator.
    #18 - 20        Residue name    resName       Residue name.
    #22             Character       chainID       Chain identifier.
    #23 - 26        Integer         resSeq        Residue sequence number.
    #27             AChar           iCode         Code for insertion of residues.
    #31 - 38        Real(8.3)       x             Orthogonal coordinates for X
    #39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
    #47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
    #55 - 60        Real(6.2)       occupancy     Occupancy.
    #61 - 66        Real(6.2)       tempFactor    Temperature factor.
    #73 - 76        LString(4)      segID         Segment identifier, left-justified.
    #77 - 78        LString(2)      element       Element symbol, right-justified.
    #79 - 80        LString(2)      charge        Charge on the atom.

    my($s) = @_;

    my @vals;
    push @vals,
      substr($s, 6,5), substr($s,12,4), substr($s,16,1), substr($s,17,3), substr($s,21,1),
      substr($s,22,4), substr($s,26,1), substr($s,30,8), substr($s,38,8), substr($s,46,8)
      ####  substr($s,54,6), substr($s,60,6), substr($s,72,4), substr($s,76,2), substr($s,78,2)
    ;

    foreach my $i (0..$#vals) {
      $vals[$i] =~s/\s+//g;
    }

    # print "Empty res name: \n$s" if $vals[3] eq '';

    return {
      AtomNum  => $vals[0],
      AtomName => $vals[1],
      #  AltInd  => $vals[2],
      ResName  => $vals[3],
      Chain  => (length $vals[4])?$vals[4]:$ONE_CHAIN, ## $vals[4] || $ONE_CHAIN doesn't work for '0'
      ResNum  => $vals[5].$vals[6],
      X   => $vals[7],
      Y   => $vals[8],
      Z   => $vals[9],
      #  Occupancy => $vals[10],
      #  Bfactor  => $vals[11],
      #  SegmId  => $vals[12],
      #  ElemSymb => $vals[13],
      #  Charge  => $vals[14],
    };

  }				# end sub
#----------------------------------------------------------------------------------------------------------------
sub get_pdb_desc {
  my $file = shift;

  if ($file =~ /\.(Z|gz)$/) {
    $file = "gunzip -qc $file |";
  }

  open(PDB, $file) or die "Can't open file: $file\n";

  my $got = 0;
  my $title;

  local $_;

  while (<PDB>) {
    if (/^TITLE/) {
      $title .= substr($_, 10, 60);
      $got++;
    } elsif ($got) {
      last;
    }
  }

  close PDB;

  if (defined $title) {
    $title =~ s/\s+/ /g;
    $title =~ s/^\s+//;
    $title =~ s/\s+$//;
  }

  return (defined $title && length $title) ? $title : 'N/A';
}
#----------------------------------------------------------------------------------------------------------------
