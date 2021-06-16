package PPH::Misc;

use Carp qw(cluck);
use strict;
use warnings;
use File::Basename;

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

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Misc

=head1 DESCRIPTION

Module containing useful routines

=head1 SUBVERSION

 $LastChangedDate: 2011-09-20 21:30:07 -0400 (Tue, 20 Sep 2011) $
 $LastChangedRevision: 368 $
 $LastChangedBy: ivan $

=head1 ROUTINES

=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD


#----------------------------------------------------------------------

=head2 locate_file

 Usage   : PPH::Misc->locate_file($filer, @paths)
 Function: Locate a (non-empty) file in several possible locations
           and set $filer to a full pathname of the first (non-empty)
           file found
 Args    : Filename (converted to a ref due to prototype)
           Array of path names (no trailing slashes)
 Returns : -1   if file found in the current dir ($$filer unaltered)
            0   if file not found ($$filer unaltered)
            N   1 <= N <= scalar(@paths) if file found in @paths
                ($$filer is set to a full pathname)

=cut

#----------------------------------------------------------------------

sub locate_file (\$@) {
  #----------------------------------------
  my ($filer, @paths) = @_;
  ##### |     PPH__Misc_locate_file: $$filer, @paths
  ####### check: -e $$filer and -s $$filer
  -e $$filer and -s $$filer and return -1;
  my $base = basename($$filer);
  my $i = 0;
  foreach (@paths) {
    $i++;
    my $pathname = $_ . '/' . $base;
    -e $pathname and -s $pathname or next;
    $$filer = $pathname;
    ####### |     returns: $i, $pathname
    return $i;
  }
  ####### |     returns: 0
  return 0;
}
#----------------------------------------------------------------------

1;
