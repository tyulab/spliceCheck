package PPH::Seq;
use strict;
use warnings;
use Carp qw(cluck confess);

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Seq

=head1 DESCRIPTION

Routines to process / convert sequence / annotation database files
(former pph_seq)

=head1 SUBVERSION

 $LastChangedDate: 2012-05-23 14:30:10 -0400 (Wed, 23 May 2012) $
 $LastChangedRevision: 398 $
 $LastChangedBy: ivan $



=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

BEGIN {
    use Exporter ();
    our ($VERSION, @ISA, @EXPORT, @EXPORT_OK, %EXPORT_TAGS);
    $VERSION = 0.1;
    @ISA           = qw(Exporter);
    @EXPORT        = qw();
    @EXPORT_OK     = qw();
    %EXPORT_TAGS = ();
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
our @EXPORT_OK;

use File::Temp;
use Storable;
use DBI;

use PPH::Align;

# UniProtKB annotations database
my $dbfile = "$CONFIG{UNIPROT}.sqlite";
my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
my $dbh;

#IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII


#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

=head1 ROUTINES

#------------------------------------------------------------------------------

=head2 id2acc

 Usage   : PPH::Seq->id2acc($id)
 Function: Converts a UniProtKB identifier / entry name into primary accession
 Args    : sequence identifier
 Returns : accession

=cut

#----------------------------------------
{
my $select = q{SELECT acc FROM id2acc WHERE name = ? LIMIT 1};
my $sth;
sub id2acc {
  my $id = shift;
  #### PPH__Seq__id2acc: $id
  die "ERROR: id2acc: Empty / missing ID\n" unless length $id;
  # Connect to database and prepare select statement if this has not been done yet
  $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  $sth = $dbh->prepare($select) unless defined $sth;
  # Search table
  $sth->execute($id);
  # Fetch results (should be a single row)
  my $acc = ($sth->fetchrow_array)[0];
  $sth->finish;
  ###### |     id2acc returns: $acc
  return $acc || '';
}
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

=head2 xdb2acc

 Usage   : PPH::Seq->xdb2acc($xacc)
 Function: maps a number of different protein identifiers from
           various databases to a unique UniProtKB accession number
 Args    : supported protein sequence identifier
 Returns : UniProtKB accession number, empty string if not matched

=cut

#------------------------------------------------------------------------------
{
my $select = q{SELECT acc FROM xdb2acc WHERE name = ?1 OR name LIKE ?2 LIMIT 1};
my $sth;
sub xdb2acc {
  my $xacc = shift;
  #### PPH__Seq__xdb2acc: $xacc
  die "ERROR: xdb2acc: Empty / missing ID\n" unless length $xacc;
  # Connect to database and prepare select statement if this has not been done yet
  $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  $sth = $dbh->prepare($select) unless defined $sth;
  # Search table
  # First try exact accession match; if it fails try appending accession
  # version number glob and then look for the first version matched
  $sth->execute($xacc, $xacc.'.%');
  # Fetch results (should be a single row)
  my $macc = ($sth->fetchrow_array)[0];
  $sth->finish;
  ###### |    xdb2acc returns: $macc
  return $macc || '';
}
}

#------------------------------------------------------------------------------

=head2 snp2acc

 Usage   : PPH::Seq->snp2acc($rsid)
 Function: map a dbSNP rsID to a sequence variant annotated in Swiss-Prot
 Args    : rsID SNP identifier
 Returns : UniProtKB accession number, variant position and the two AA residues

=cut

#------------------------------------------------------------------------------
{
# We only fetch the first matching variant (can be several)
my $select = q{SELECT acc, pos, aa1, aa2 FROM snp2acc WHERE rsid = ? LIMIT 1};
my $sth;
sub snp2acc {
  my $rsid = shift;
  #### PPH__Seq__snp2acc: $rsid
  die "ERROR: snp2acc: Empty / missing rsID\n" unless length $rsid;
  # Connect to database and prepare select statement if this has not been done yet
  $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  $sth = $dbh->prepare($select) unless defined $sth;
  # Search table
  $sth->execute($rsid);
  # Fetch results (should be a single row)
  my @variant = $sth->fetchrow_array;
  $sth->finish;
  ###### |    snp2acc returns: @variant
  return @variant;
}
}

#------------------------------------------------------------------------------

=head2 acc2snp

 Usage   : PPH::Seq->acc2snp($acc, $pos, $aapair)
 Function: look up a dbSNP rsID for a substitution in UniProtKB annotations
 Args    : UniProtKB protein accession, substitution position, sorted pair of AA residues
 Returns : dbSNP rsID

=cut

#------------------------------------------------------------------------------
{
# We only fetch the first matching rsID (can be several)
my $select = q{SELECT rsid FROM snp2acc WHERE acc = ? AND pos = ? AND aapair = ? LIMIT 1};
my $sth;
sub acc2snp {
  my ($acc, $pos, $aapair) = @_;
  #### PPH__Seq__acc2snp: $acc, $pos, $aapair
  # Connect to database and prepare select statement if this has not been done yet
  $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  $sth = $dbh->prepare($select) unless defined $sth;
  # Search table
  $sth->execute($acc, $pos, $aapair);
  # Fetch results (should be a single row)
  my $rsid = ($sth->fetchrow_array)[0];
  $sth->finish;
  ###### |    acc2snp returns: $rsid
  return $rsid;
}
}

#------------------------------------------------------------------------------

=head2 fastacmd

  Usage   : PPH::Seq->fastacmd($id, $db)
  Function: Fetch a sequence from a database prepared with formatdb (NCBI)
  Args    : Identifier
            Database
  Returns : Unprocessed header
            Accession no
            Identifier
            Description
            The sequence (without white spaces)

=cut

#----------------------------------------
sub fastacmd {
  my ($id, $db) = @_;
  #### PPH__Seq__fastacmd: @_
  my ($head, $seq);
  eval {
    my $pid;
    my $tmpf = File::Temp->new();
    # IAA: switched to blastdbcmd (NCBI BLAST+) because of fastacmd random crashes
    system( "$CONFIG{FASTACMD}",
      '-db', $db, '-dbtype', 'prot', '-entry', $id, '-out', $tmpf, '-outfmt', '%f' ) == 0 or
      die "Failed to execute $CONFIG{FASTACMD}";
    open my $F, '<', $tmpf or die "Can't open file: $tmpf";
    $head = <$F>;
    $seq  = join '', <$F>;
    close $F;
    defined $head or die 'Empty sequence returned';
    chomp $head;
    $head =~ /^>/ or die 'Extracted FASTA entry definition line malformed or missing';
    $seq  =~ s/\s+//g;
    defined $seq && length $seq or die 'Empty sequence returned';
  };
  if ($@) {
    warn "ERROR: fastacmd: $@\n";
    return;
  }

  my ($def, $dbid, $acc, $name, $desc) = PPH::Align::_parse_defline($head);

  ##### |     fastacmd returns: $head, $acc, $name, $desc, $seq
  return ($head, $acc, $name, $desc, $seq);
}
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

=head2 fetch_FASTA

 Usage   : PPH::Seq->fetch_FASTA($id, $db)
 Function: Fetch a sequence from a database
 Args    : Sequence identifier
           FASTA database file name
 Returns : Accession no
           Identifier
           Description
           The sequence (without white spaces)

=cut

#----------------------------------------
sub fetch_FASTA {
    my ($id, $db) = @_;
    #### PPH__Seq__fetch_FASTA: $id, $db

    ##### |    $CONFIG{FETCHSEQ}: $CONFIG{FETCHSEQ}
    open my $F, '-|', qq($CONFIG{FETCHSEQ} '$db' $id)
      or confess "Can't execute fetchseq $db $id\n";

    my $output = do { local $/; <$F> }; close $F;

    return undef unless defined $output && length $output;

    my ($head, $seq) = split /\n/m, $output, 2;
    my ($def, $dbid, $acc, $name, $desc) = PPH::Align::_parse_defline($head);

    # Do some sequence cleanup
    $seq =~ s/\s+//g;
    # Users often submit sequences with translated terminator codon(s)
    $seq =~ s/\*+$//;
    $seq =~ s/X+$//i;
    # Terminator codons inside query sequence look dubious
    warn "WARNING: Terminator codon(s) encountered in query protein sequence: $acc\n" if $seq =~ tr/*/X/;

    ##### |     fetch_FASTA returns: $acc, $name, $desc, $seq
    return ($acc, $name, $desc, $seq);
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 fetch_txNuc

 Usage   : PPH::Seq->fetch_txNuc($txname, $nuc2bit)
 Function: Fetch full nucleotide sequence for a transcript from 2bit database
 Args    : Transcript identifier (name)
 Returns : Nucleotide sequence (without whitespace)

=cut

#----------------------------------------
sub fetch_txNuc {
  my ($txname, $nuc2bit) = @_;
  #### PPH__Gene__fetch_txNuc: $txname, $nuc2bit
  my $nucseq  = `$CONFIG{TWOBIT} -seq='$txname' $nuc2bit stdout 2>&1`;
  my $defline = $1 if $nucseq =~ s/^>(\S+).*?\n//s;
  if (defined $defline) {
    # correct transcript defline
    if ($defline eq $txname) {
      $nucseq =~ s/\s+//g;
      # empty sequence??
      unless (length $nucseq) {
        warn "ERROR: fetch_txNuc: Empty nucleotide sequence for transcript: $txname\n";
        return;
      }
    # wrong transcript??
    } else {
      warn "ERROR: fetch_txNuc: Incorrect transcript fetched from the sequence database ($defline), should be: $txname\n";
      return;
    }
  # something went wrong, error message returned instead of the sequence
  } else {
    $nucseq =~ s/\s+/ /g;
    warn "ERROR: fetch_txNuc: Failed to extract nucleotide sequence for transcript $txname: $nucseq\n";
    return;
  }
  #### PPH__Gene__fetch_txNuc returns: $nucseq
  return $nucseq;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 get_uniprot_annotation

 Usage   : my $status = $PPH::Seq::get_uniprot_annotation($PROT, $SNP)
 Function: Find UniProt annotations for a given protein / SNP
 Args    : Objects of type PROT and SNP
 Returns : $SNP->{Features} = { CritSites, Sites, Regions, PHAT }
           1 if ok

=cut

#----------------------------------------
sub get_uniprot_annotation {
    my ($PROT, $SNP) = @_;
    ### PPH__Seq__get_uniprot_annotation: $PROT, $SNP

    # Load annotations for critical sites, used later to calculate SNP 3D contacts
    unless (exists $PROT->{Features}{CritSites}) {
      my $csites = _search_features($PROT->{Acc});
      if (@$csites) {
        my %critsites;
        foreach my $csite (@$csites) {
          my ($name, $from, $to) = @$csite;
          # BONDS range consists of the two connection points only, otherwise
          # it is a list of all consecutive position between "from" and "to"
          push @{ $critsites{$name} },
            exists $PPH::Data::BONDS{$name} ? ($from, $to) : ($from..$to);
        }
        $PROT->{Features}{CritSites} = \%critsites;
      } else {
        $PROT->{Features}{CritSites} = {};
      }
      ###### |     get_uniprot_annotation CritSites: $PROT->{Features}{CritSites}
    }

    # Rest of features are mapped to SNP position
    my $features = _search_features($PROT->{Acc}, $SNP->{Pos});
    return 0 unless @$features;

    my $transmem_found;
    foreach my $ftname (@$features) {
      # Is feature a region...
      if (exists $PPH::Data::REGIONS{$ftname}) {
        push @{ $SNP->{Features}{Regions} }, $ftname;
        # Transmembrane region
        $transmem_found++ if $ftname eq 'TRANSMEM';
      }
      # ...or a site?
      elsif (exists $PPH::Data::SITES1{$ftname} or
             exists $PPH::Data::SITES2{$ftname}) {
        push @{ $SNP->{Features}{Sites} }, $ftname;
      }
    }
    if ($transmem_found) {
      my %PHAT = PPH::Data::blastmatrix("$CONFIG{MATRIX}/PHAT_T75_B73.bla");
      $SNP->{Features}{PHAT} = $PHAT{$SNP->{Aa1}}{$SNP->{Aa2}};
    }

    ###### |     get_uniprot_annotation returns (Features): $SNP->{Features}
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head3 _search_features

 Usage   : _search_features($SNP->{Acc}[, $SNP->{Pos}]);
 Function: Search for annotated protein features in the database
 Args    : Accession, SNP position (optional)
 Returns : Ref to list of features found

=cut

#----------------------------------------
{
my $bonds  = '"' . join('","', keys %PPH::Data::BONDS)  . '"';
my $csites = '"' . join('","', keys %PPH::Data::SITES1) . '"';
my @select = (
# Select all features at the SNP position
qq{SELECT name, CASE WHEN name IN($bonds) THEN 1 ELSE 0 END AS isbond FROM features
WHERE acc = ?1 AND (isbond AND (start = ?2 OR end = ?2) OR NOT isbond AND ?2 BETWEEN start AND end)},
# Select all critical sites in the protein
qq{SELECT name, start, end FROM features WHERE acc = ? AND name IN($csites)}
);
my @sth;
sub _search_features {
  my ($acc, $pos) = @_;
  #### PPH__Seq__search_features: $acc, $pos
  die "ERROR: _search_features: Empty / missing accession\n" unless length $acc;
  my $sx = defined $pos ? 0 : 1;
  # Connect to database and prepare select statement if this has not been done yet
  $dbh      = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  ####### |     SQL statement: $select[$sx]
  $sth[$sx] = $dbh->prepare($select[$sx]) unless defined $sth[$sx];
  # Search table and fetch results
  my @features;
  if ($sx) {  # fetches FT name and positional range (all critical sites in the protein)
    $sth[$sx]->execute($acc);
    while (my @feature = $sth[$sx]->fetchrow_array) {
      push @features, \@feature;
    }
  } else {    # fetches FT name only (names of bonds, sites, and regions at the SNP position)
    $sth[$sx]->execute($acc, $pos);
    while (my @ft = $sth[$sx]->fetchrow_array) {
      push @features, $ft[0];
    }
  }
  $sth[$sx]->finish;
  ###### |     _search_features returns: @features
  return \@features;
}
}


#------------------------------------------------------------------------------

=head2 grep_in_FASTA

 Usage   : my ($acc, $name, $desc, $seq) = grep_in_FASTA($query_seq, $target_db)
 Function: Find a (sub) sequence in a FASTA-formatted database file
           File has to keep each sequence on ONE LINE!
 Args    : string with the query protein sequence
           target FASTA database filename
 Returns : Returns FIRST hit which contains the match to query sequence
           Accession
           Entry name
           Description
           Sequence

=cut

#------------------------------------------------------------------------------
sub grep_in_FASTA {
  my ($query_seq, $target_db) = @_;
  #### PPH__Seq__grep_in_FASTA: $query_seq, $target_db

  my $pid = open my $F, '-|', "fgrep -s -m 1 -B 1 '$query_seq' $target_db";
  die "ERROR: grep: $!\n" unless defined $pid;
  my $head = <$F>;
  close($F), return () unless defined $head;
  my $seq = <$F>;
  close $F or die "\n";

  chomp $head;
  chomp $seq;

  $head =~ s/^>// or die "ERROR: Incorrect FASTA defline format: $head\n";
  my ($definition, $db, $acc, $name, $desc) = PPH::Align::_parse_defline($head);

  ##### |     grep_in_FASTA returns: $acc, $name, $desc, $seq
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

#------------------------------------------------------------------------------

=head2 map_by_BLAT

 Usage   : @hits = map_by_BLAT($qacc, $qseq, $target_db)
 Function: Find closely matching sequences in a FASTA-formatted database file
 Args    : query sequence in FASTA format (string)
           target FASTA database filename
 Returns : Returns ref to array of hits with matches to query sequence
           Accession
           Entry name
           Description
           Sequence

=cut

#------------------------------------------------------------------------------
{
my $ident_level   = 0.97;
my $overlap_level = 0.20;
my $no_of_hsps    = 0;    # returns all hits
sub map_by_BLAT {
  my ($qacc, $qseq, $db) = @_;
  #### PPH__Gene__map_by_BLAT: $qacc, $qseq, $db

  chomp $qacc;
  die "ERROR: map_by_BLAT: Missing or empty query accession\n" unless defined $qacc && length $qacc;
  chomp $qseq;
  die "ERROR: map_by_BLAT: Missing or empty query sequence\n" unless defined $qseq && length $qseq;
  my $qfa  = '>' . $qacc . "\n" . $qseq . "\n";
  my $qlen = length $qseq;
  my ($macc, $mid, $mlen, $overlap, $hits, $qstart, $qend, $mstart, $mend, $qfrac, $bsizes, $qstarts, $mstarts);
  my ($mdesc, $mcan, $mseq);
  my $ident = 0;

  my @hsps;
  my $mapped = 0;
  foreach my $hsp (run_BLAT($qfa, $db)) {
    ($macc, $mlen, $ident, $overlap, $hits, $qstart, $qend, $mstart, $mend, $bsizes, $qstarts, $mstarts) = @$hsp;
    if (defined $macc) {
      $ident = sprintf("%.3g", $ident);
      $qfrac = sprintf "%.3g", $overlap / $qlen;  # fraction of query sequence aligned
      $mdesc = (grep_FASTA_defline($macc, $db))[2];
      $mcan  = $mdesc =~ /^Canonical;/ ? 1 : 0;
      for ($bsizes, $qstarts, $mstarts) { s/,$//; }
      push @hsps, [ $qacc, $macc, $qlen, $mlen, $overlap, $hits, $qstart, $qend, $mstart, $mend, $qfrac, $ident, $mcan,
                    $bsizes, $qstarts, $mstarts ];
      $mapped  = 1;
    }
  }

  my @matched;
  my $hsp_count = 0;
  if ($mapped) {
    foreach my $hsp (
        # Sort on:
        #   canonical target hits first, mcan (desc)
        #   fraction of query sequence in overlap (desc),
        #   number of alignment blocks (asc)
        #   fraction of identity in overlap (desc),
        sort { $b->[12] <=> $a->[12] || $b->[10] <=> $a->[10] || $a->[5] <=> $b->[5] || $b->[11] <=> $a->[11] } @hsps
      ) {
      last if $no_of_hsps && $hsp_count >= $no_of_hsps;
      if ($hsp->[11] >= $ident_level && $hsp->[10] >= $overlap_level) {
        push @matched, $hsp;
        $hsp_count++;
      }
    }
    ##### |     map_by_BLAT returns: @matched
    return @matched ? \@matched : undef;
  } else {
    ##### |     map_by_BLAT found no matches
    return;
  }
}
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
{
my $pident  = 97;   # minimum percentage of sequence identity
sub run_BLAT {
  my ($query, $db) = @_;
  #### PPH__Gene__run_BLAT: $query, $db
  return unless $CONFIG{'BLAT'};
  my $queryfa = File::Temp->new();
  print $queryfa $query;
  $queryfa->flush;
  # usage: blat [options] database query output.psl
  my $cmd = qq($CONFIG{BLAT} -prot -minIdentity=$pident -out=psl -noHead $db $queryfa stdout);
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
    # macc mlen ident length hits qstart qend mstart mend bsizes qstarts mstarts (12)
    push @hsps, [
      $macc,
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
}
#------------------------------------------------------------------------------


#MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

1;
