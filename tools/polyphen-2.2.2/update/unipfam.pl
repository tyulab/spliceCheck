#!/usr/bin/env perl
use warnings;
use strict;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

unipfam.pl - a simple script to extract Pfam annotations for UniProt proteins

=head1 SYNOPSIS

unipfam.pl [options] [<seq_file>]

Where options are:

 -n <organism>    specify common organism name instead of the default human;
                  supported names are: chimp, orangutan, macaque (crab-eating),
                  mouse, rat, dog, cow, zebrafish and fruitfly; some taxonomic
                  groups of species are also supported, based on UniProtKB set
                  of taxonomic divisions. For a full list, please see:
                    ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/taxonomic_divisions/
 -h               print this help

=head1 DESCRIPTION

This script uses UniProtKB sequence file (default: human.seq) in
FASTA format prepared by uniprot.pl script to get a list of protein
accessions and then downloads and parses swisspfam.gz file to extract
Pfam annotations for all proteins found in the sequence file.

=head1 SUBVERSION

 $LastChangedDate: 2014-05-23 21:09:39 -0400 (Fri, 23 May 2014) $
 $LastChangedRevision: 470 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use LWP::Simple;
use Getopt::Std;
use Pod::Usage;
use DBI;

my %opts;
getopts('hln:', \%opts) or pod2usage();
exists $opts{'h'} and pod2usage({ -verbose=>99 });

# Organism common name: human (default), orangutan, macaque (crab-eating), mouse, rat, dog, cow, zebrafish
# For UniProt organism codes (spname), see: http://www.uniprot.org/docs/speclist
my $CNAME = 'human';

# SPECIES hash indexed by CNAME
my %SPECIES = (
  # Human
  'human'    => {
    'os'     => 'Homo sapiens',
    'spname' => 'HUMAN',
    'file'   => 'human',
  },
  # Chimpanzee
  'chimp'    => {
    'os'     => 'Pan troglodytes',
    'spname' => 'PANTR',
    'file'   => 'mammals',
  },
  # Orangutan
  'orangutan'    => {
    'os'     => 'Pongo pygmaeus|Pongo abelii',
    'spname' => 'PONPY|PONAB',
    'file'   => 'mammals',
  },
  # Crab-eating macaque
  'macaque'    => {
    'os'     => 'Macaca fascicularis',
    'spname' => 'MACFA',
    'file'   => 'mammals',
  },
  # Mouse
  'mouse'    => {
    'os'     => 'Mus musculus',
    'spname' => 'MOUSE',
    'file'   => 'rodents',
  },
  # Rat
  'rat'      => {
    'os'     => 'Rattus norvegicus|Rattus rattus',
    'spname' => 'RAT|RATRT',
    'file'   => 'rodents',
  },
  # Dog
  'dog'      => {
    'os'     => 'Canis lupus familiaris',
    'spname' => 'CANLF',
    'file'   => 'mammals',
  },
  # Cow
  'cow'      => {
    'os'     => 'Bos taurus',
    'spname' => 'BOVIN',
    'file'   => 'mammals',
  },
  # Zebrafish
  'zebrafish' => {
    'os'      => 'Danio rerio',
    'spname'  => 'DANRE',
    'file'    => 'vertebrates',
  },
  # Fruit fly
  'fruitfly' => {
    'os'      => 'Drosophila melanogaster',
    'spname'  => 'DROME',
    'file'    => 'invertebrates',
  },
  #
  # Taxonomic divisions
  #
  # bacteria
  'bacteria' => {
    'os'      => '',
    'spname'  => '',
    'file'    => 'bacteria',
  },
);

if (my $cname = $opts{'n'}) {
  die "Illegal or unsupported organism name: $cname\n"
    unless exists $SPECIES{$cname};
  $CNAME = $cname;
}

my $uniseq_name = "$CNAME.seq";
$uniseq_name    = shift if @ARGV;
my $outfile     = "$CNAME.pfam";
my $dbfile      = "$CNAME.sqlite";

$| = 1;

# Load precompiled pfam data into SQLite database and quit
createdb($outfile, $dbfile), exit(0) if $opts{'l'};

open(F, $uniseq_name) or die "Can't open input file: $uniseq_name\n";

my %HsAcc;
print "Reading protein accessions ...\n";
while (<F>) {
    /^>/ or next;
    my (undef, $acc) = split /\|/;
    next if $acc =~ /-\d+$/;  # skip isoforms
    $HsAcc{$acc} = 1;
}
close F;
print "Done.\n";

my $fetched = 0;
# Reuse existing swisspfam.gz file
if (-e 'swisspfam.gz') {
  print "Found swisspfam.gz in the current directory, reusing.\n";
# Download recent swisspfam.gz file
} else {
  my $PFAM_FTP = 'ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/swisspfam.gz';
  print "Fetching swisspfam.gz ...\n";
  die "ERROR: Can't fetch file: swisspfam.gz\n" unless mirror($PFAM_FTP, 'swisspfam.gz') == 200;
  print "Done.\n";
  $fetched = 1;
}

my $pcount = 0;
my $dcount = 0;
print "Extracting Pfam annotations ...\n";
open(F, 'gunzip -c swisspfam.gz |') or die "ERROR: Can't open input file: swisspfam.gz\n";
open(FOUT, ">$outfile") or die "ERROR: Can't create output file: $outfile\n";
while (<F>) {
  /^>/ or next;
  my (undef, undef, $acc) = split;
  $acc =~ s/\.\d+$//; # remove sequence version number
  exists $HsAcc{$acc} or next;
  while (<F>) {
    /^\s*$/   and last;
    /^Pfam-B/ and next;
    chomp;
    my @ranges;
    while (s/\s(\d+-\d+)$//) { unshift @ranges, $1; }
    s/\s+$//;
    my ($pfname, $desc) = split ' ', $_, 2;
    my $pfid = $1 if $desc =~ s/^.+\(\d+\)\s+(PF\d+\.\d+)\s+//;
    warn("ERROR: Failed to locate Pfam identifier for $acc entry: $_\n"), next unless defined $pfid;
    $desc =~ s/\s+$//; $desc =~ s/^\s+//;
    # acc ranges pfid pfname desc
    print FOUT join("\t", $acc, join(',', @ranges), $pfid, $pfname, $desc), "\n";
    $dcount++;
  }
  $pcount++;
}
close FOUT;
close F;
unlink 'swisspfam.gz' if -e 'swisspfam.gz' && $fetched;
print "Found $dcount Pfam-A domains in $pcount $CNAME proteins.\n";
# Load compiled pfam file into SQLite database
createdb($outfile, $dbfile);
# Remove database dump file
unlink $outfile or die "ERROR: Can't delete file: $outfile\n";
print "All done.\n";

#---------------------------------------------------------------------------------------

# Create and populate SQLite database
#
#   scale for start-end intervals is 1-based
#
#   acc     start   end   pfid         pfname  desc
#   P63104      3   236   PF00244.11   14-3-3  14-3-3 protein 1
#   P63104    450   487   PF00244.11   14-3-3  14-3-3 protein 1
#
sub createdb {
  my ($datafile, $dbfile) = @_;
  die "ERROR: createdb: Missing data file name\n" unless length $datafile;
  die "ERROR: createdb: Missing database file name\n" unless length $dbfile;
  die "ERROR: createdb: Missing or empty data file: $datafile\n" unless -s $datafile;
  my $table = q{
  CREATE TABLE pfam (
    acc    VARCHAR(16) NOT NULL,
    start  INTEGER NOT NULL,
    end    INTEGER NOT NULL,
    pfid   VARCHAR(16) NOT NULL,
    pfname VARCHAR(20),
    desc   VARCHAR(60)
  )
  };
  my $insert = q{INSERT INTO pfam VALUES (} . join(', ', ('?')x6) . ')';
  my @indices = (
    q{CREATE INDEX x_acc  ON pfam (acc)},
    q{CREATE INDEX x_pfid ON pfam (pfid)},
    q{CREATE INDEX x_ast  ON pfam (acc,start)},
    q{CREATE INDEX x_aen  ON pfam (acc,end)},
  );
  print 'Creating pfam table ...';
  my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
  my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs);
  # Set some pragmas
  $dbh->do('PRAGMA legacy_file_format = OFF');
  $dbh->do('PRAGMA foreign_keys = ON');
  # Create table
  $dbh->do('DROP TABLE IF EXISTS pfam');
  $dbh->do($table);
  print " done.\n";
  # Prepare insert statement
  my $sth = $dbh->prepare($insert);
  # Polupate table
  print 'Populating pfam table ';
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
    my ($acc, $ranges, $pfid, $pfname, $desc) = split /\t/;
    my @r = split /,/, $ranges; # ranges
    foreach my $range (@r) {
      my ($start, $end) = split /-/, $range;
      eval { $sth->execute($acc, $start, $end, $pfid, $pfname, $desc) };
      if ($@) {
        warn $@;
        $dbh->rollback;
        print "\n" if $count;
        exit 1;
      }
      $count++;
    }
  }
  close(FIN);
  print "\n";
  $dbh->commit if $count;
  # Index table
  print 'Indexing pfam table   ';
  foreach my $index (@indices) {
    print '.';
    $dbh->do($index);
  }
  print "\n";
  print 'Optimizing pfam table ...';
  $dbh->do('ANALYZE pfam');
  print "\n";
  $dbh->disconnect;
  print "Finished.\n";
}
