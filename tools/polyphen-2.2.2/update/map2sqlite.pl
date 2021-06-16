#!/usr/bin/env perl

use warnings;
use strict;

use DBI;

my $usage = <<'EOT';
Usage: map2sqlt.pl <seqmap_output> <sqlite_database>
EOT

my @tables = qw( map hits );
my @schema = (
q{
CREATE TABLE map (
  id      INTEGER PRIMARY KEY AUTOINCREMENT,
  qacc    VARCHAR(16) NOT NULL,
  macc    VARCHAR(16) NOT NULL,
  mname   VARCHAR(20),
  qlen    INTEGER NOT NULL,
  mlen    INTEGER NOT NULL,
  overlap INTEGER NOT NULL,
  hits    INTEGER NOT NULL,
  qstart  INTEGER NOT NULL,
  qend    INTEGER NOT NULL,
  mstart  INTEGER NOT NULL,
  mend    INTEGER NOT NULL,
  qfrac   REAL NOT NULL,
  ident   REAL NOT NULL,
  method  CHAR(1) NOT NULL,
  qcan    BOOLEAN NOT NULL,
  mcan    BOOLEAN NOT NULL
)
},
q{
CREATE TABLE hits (
  id      INTEGER REFERENCES map,
  qfrom   INTEGER NOT NULL,
  qto     INTEGER NOT NULL,
  mfrom   INTEGER NOT NULL,
  mto     INTEGER NOT NULL
)
}
);

my $insert_map  = 'INSERT INTO map VALUES (NULL, ' . join(', ', ('?')x16) . ')';
my $insert_hits = 'INSERT INTO hits VALUES ('      . join(', ', ('?')x5) . ')';

my @indices = (
  'CREATE INDEX x_qacc    ON map (qacc)',
  'CREATE INDEX x_macc    ON map (macc)',
  'CREATE INDEX x_mname   ON map (mname)',
  'CREATE INDEX x_overlap ON map (overlap DESC)',
  'CREATE INDEX x_hits    ON map (hits ASC)',
  'CREATE INDEX x_qfrac   ON map (qfrac DESC)',
  'CREATE INDEX x_ident   ON map (ident DESC)',
  'CREATE INDEX x_iqfrom  ON hits (id, qfrom)',
  'CREATE INDEX x_iqto    ON hits (id, qto)',
  'CREATE INDEX x_imfrom  ON hits (id, mfrom)',
  'CREATE INDEX x_imto    ON hits (id, mto)'
);

die $usage unless @ARGV == 2;

my $datafile = shift @ARGV;
my $dbfile   = shift @ARGV;
open(FIN, $datafile) or die "ERROR: Can't open input file $datafile: $!\n";

my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
my $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs);

# set some pragmas, SQLite-specific!
$dbh->do('PRAGMA legacy_file_format = OFF');
$dbh->do('PRAGMA foreign_keys = ON');

# (re-)create table
$dbh->do('DROP TABLE IF EXISTS hits');
$dbh->do('DROP TABLE IF EXISTS map');
for my $table (@schema) {
  $dbh->do($table);
}

my $sth_map  = $dbh->prepare($insert_map);
my $sth_hits = $dbh->prepare($insert_hits);

$| = 1;

my $batch_size = 5000;
my $count = 0;
print 'Loading    ';
while (<FIN>) {
  next if /^#/ || /^\s*$/;
  chomp;
  unless ($count % $batch_size) {
    if ($count) { print '.'; $dbh->commit; }
    $dbh->begin_work;
  }
  my @r = split /\t/;
  my @a = @r[0..$#r-3];
  my ($bsizes, $qstarts, $mstarts) = @r[-3..-1];
  my @bsizes  = split /,/, $bsizes;
  my @qstarts = split /,/, $qstarts;
  my @mstarts = split /,/, $mstarts;
  my @qends;  # all ends are zero-based closed [from,to]
  my @mends;  # all ends are zero-based closed [from,to]
  for (my $i=0; $i<@bsizes; $i++) {
    # end = start + length - 1 (zero-based closed end scale)
    push @qends, $qstarts[$i] + $bsizes[$i] - 1;
    push @mends, $mstarts[$i] + $bsizes[$i] - 1;
  }
  eval { $sth_map->execute(@a); };
  if ($@) {
    warn $@;
    $dbh->rollback;
    print "\n" if $count;
    exit 1;
  }
  my $hit = $dbh->last_insert_id('', '', '', '');
  eval {
    for (my $i=0; $i<@bsizes; $i++) {
      $sth_hits->execute($hit, $qstarts[$i], $qends[$i], $mstarts[$i], $mends[$i]);
    }
  };
  if ($@) {
    warn $@;
    $dbh->rollback;
    print "\n" if $count;
    exit 1;
  }
  $count++;
}

print "\n";
$dbh->commit;

print 'Indexing   ';
foreach my $index (@indices) {
  print '.';
  $dbh->do($index);
}
print "\n";

print 'Optimizing ';
foreach my $table (@tables) {
  $dbh->do("ANALYZE $table");
}
print "\n";

$dbh->disconnect;

close(FIN);
