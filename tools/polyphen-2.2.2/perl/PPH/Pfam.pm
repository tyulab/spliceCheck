package PPH::Pfam;

use strict;
use warnings;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

Pfam - Using the swisspfam file to find out it a SNP lies within a domain

=head1 DESCRIPTION



=head1 AUTHOR

Steffen Schmidt

=head1 SUBVERSION

 $LastChangedDate: 2011-11-21 13:28:37 -0500 (Mon, 21 Nov 2011) $
 $LastChangedRevision: 376 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use DBI;

BEGIN {
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

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

=head1 ROUTINES

=cut

# Pfam annotations database
my $dbfile = "$CONFIG{UNIPROT}.sqlite";
my $dbargs = { AutoCommit=>1, RaiseError=>1, PrintError=>0 };
my $dbh;

#------------------------------------------------------------------------------

#   Table: pfam
#   scale for start-end intervals is 1-based
#
#   acc     start   end   pfid         pfname  desc
#   P63104      3   236   PF00244.11   14-3-3  14-3-3 protein 1
#   P63104    450   487   PF00244.11   14-3-3  14-3-3 protein 1
#
{
my $select = q{SELECT pfid FROM pfam WHERE acc = ?1 AND ?2 BETWEEN start AND end};
my $sth;
sub _search_pfam {
  my ($acc, $pos) = @_;
  ##### PPH__Pfam__search_pfam: $acc, $pos
  $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  $sth = $dbh->prepare($select) unless defined $sth;
  # Search table
  $sth->execute($acc, $pos);
  # Fetch results (can be several rows or empty)
  my @pfids;
  while (my @id = $sth->fetchrow_array) {
    push @pfids, @id;
  }
  $sth->finish;
  ###### |    _search_pfam returns: @pfids
  return \@pfids;
}
}

#------------------------------------------------------------------------------

=head2 find_pfam

  Usage    : $PPH::Pfam::find_pfam($PROT, $SNP)
  Function : Find matching Pfam annotations for a SNP
  Args     : Objects of types PRTO and SNP
  Sets     : $SNP->{Features}{Domains}[ from, to, pfid ]
=cut

#----------------------------------------

sub find_pfam {
    my ($PROT, $SNP) = @_;
    ### PPH__Pfam_find_pfam: $PROT->{Acc}, $SNP->{Pos}
    $PROT->{Pfam} = _in_pfam($PROT->{Acc}) unless exists $PROT->{Pfam};
    if ($PROT->{Pfam}) {
      my $pfam = _search_pfam($PROT->{Acc}, $SNP->{Pos});
      foreach my $pfid (@$pfam) {
        push @{ $SNP->{Features}{Domains} }, $pfid;
      }
    }
    ##### |     find_pfam returns: $SNP->{Features}{Domains}
}

#------------------------------------------------------------------------------

=head2 _in_pfam

  Usage    : $PPH::Pfam::_in_pfam($PROT->{Acc})
  Function : Check if the protein is present in Pfam database
  Args     : UniProtKB protein accession
  Returns  : 0 | 1
=cut

#----------------------------------------
{
my $select = q{SELECT acc FROM pfam WHERE acc = ? LIMIT 1};
my $sth;
sub _in_pfam {
  my $qacc = shift;
  ##### PPH__Pfam__in_pfam: $qacc
  $dbh = DBI->connect("dbi:SQLite:dbname=$dbfile", '', '', $dbargs) unless defined $dbh;
  $sth = $dbh->prepare($select) unless defined $sth;
  # Search table
  $sth->execute($qacc);
  # Fetch results (should be a single row)
  my $macc = ($sth->fetchrow_array)[0];
  $sth->finish;
  ###### |    _in_pfam returns: $macc
  return defined $macc ? 1 : 0;
}
}

1;
