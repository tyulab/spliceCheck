package SNP;
use warnings;
use strict;
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

SNP

=head1 DESCRIPTION

Object containing all information about a SNP

=head1 AUTHOR

StS

=head1 SUBVERSION

 $LastChangedDate: 2009-11-28 00:15:09 -0500 (Sat, 28 Nov 2009) $
 $LastChangedRevision: 254 $
 $LastChangedBy: ivan $

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

use Carp qw(cluck confess);
use XML::Simple;
use Data::Dumper;
use Storable;

#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
# Begin
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
#BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB

#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# Constructor
sub new {
  #### SNP object constructor: @_
  my ($class, $params, $comments) = @_;
  my ($acc, $pos, $aa1, $aa2, $txname, $cdspos, $frame, $nt1, $nt2, $flanks, $transv, $cpg, $jxdon, $jxacc) = @$params;
  my $self;
  # Do we solely have a dbSNP rsID in input?
  if ($acc =~ /^rs/ && ! grep { defined } $pos, $aa1, $aa2) {
    # rsID only
    $self = {
      rsID => $acc,
    };
  } else {
    # First 4 input params are mandatory, rest are optional
    $self = {
      Acc  => $acc,
      Pos  => $pos,
      Aa1  => $aa1,
      Aa2  => $aa2
    };
  }
  # Store a copy of the original query.
  # This simple copy method only works for single level hash!
  $self->{Query} = { %$self };
  # Add optional nucleotide sequence based parameters.
  # Note: The values are NOT validated since source program (e.g. mapsnps.pl) is
  #       supposed to had them validated already.
  # Note: Nucleotide input parameters are not stored in the Query section but
  #       instead recorded to Gene section of the SNP object because if present
  #       in input, they are never re-mapped or altered in any other way by
  #       the program.
  my @nucvalues = ($txname, $cdspos, $frame, $nt1, $nt2, $flanks, $transv, $cpg, $jxdon, $jxacc);
  my @nucpars = qw{ TxName CdnPos Frame Nt1 Nt2 Flanks Transv CpG JXdon JXacc };
  foreach my $value (@nucvalues) {
    my $name = shift @nucpars;
    $self->{Gene}{$name} = $value if defined $value && length $value;
  }
  # Add optional comments
  $self->{Query}{Comments} = $comments if defined $comments && length $comments;
  bless($self, $class);
  my $msg = $self->_test_init();
  die "$msg\n" if $msg;
  ###### |     SNP object constructor returns: $self
  return $self;
}

sub new_empty {
  my ($class, $self) = @_;
  $self = {};
  bless($self, $class);
  return $self;
}
#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
# Routines

#----------------------------------------------------------------------

=head1 ROUTINES

=head3 _test_init

 Usage   : _test_init()
 Function: Simply test if initial SNP is ok
 Args    : None (Object of type SNP)
 Returns : TRUE if successful

=cut
#----------------------------------------
sub _test_init {
  my $self = shift;

  # If there is ONLY a dbSNP rsID present in input then there is nothing to test
  return '' if exists $self->{rsID} && $self->{rsID};

  # Check for mandatory parameters
  defined $self->{Acc} and length $self->{Acc} or return 'Protein accession missing';
  defined $self->{Pos} and length $self->{Pos} or return 'Substitution position missing';
  defined $self->{Aa1} and length $self->{Aa1} or return 'Reference residue (AA1) missing';
  defined $self->{Aa2} and length $self->{Aa2} or return 'Variant residue (AA2) missing';

  $self->{Pos} =~ /^\d+$/                    or return "Illegal position specified: $self->{Pos}";
  $self->{Aa1} =~ /^[GPAVLIMCFYWHKRQNEDST]$/ or return "Illegal AA1 residue code: $self->{Aa1}";
  $self->{Aa2} =~ /^[GPAVLIMCFYWHKRQNEDST]$/ or return "Illegal AA2 residue code: $self->{Aa2}";

  $self->{Aa1} ne $self->{Aa2}               or return "Synonymous substitution: $self->{Aa1} -> $self->{Aa2}";

  return '';
}
#----------------------------------------------------------------------

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
1;
