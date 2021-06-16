package PPH::BLAST;
#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::BLAST

=head1 DESCRIPTION

Routines to BLAST and Parse
(former pph_misc.pl)


=head1 SUBVERSION

 $LastChangedDate: 2011-11-21 13:28:37 -0500 (Mon, 21 Nov 2011) $
 $LastChangedRevision: 376 $
 $LastChangedBy: ivan $


=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
use Carp qw(carp cluck confess);
use strict;
use warnings;

use File::Basename;
use File::Copy;
use XML::Simple;

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
    # Make sure old BLAST environment does not get in the way
    delete @ENV{'BLAST_PATH', 'BLASTDB', 'BLASTMAT'};
}

use PPH::LEON;

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

=head1 ROUTINES

=head2 run_BLAST (pph_psic runblast)

 Usage   : PPH::BLAST->runblast($seqfile, $blastfile, $blast_arg)
 Function: Simple wrapper to execute a BLAST command
 Args    : Sequence filename
           BLAST    filename
           BLAST    arguments
 Returns : TRUE if successful

=cut

#----------------------------------------
sub run_BLAST {
    my ($seq_file, $blast_file, $blast_arg) = @_;
    #### PPH__BLAST__run_BLAST: $seq_file, $blast_file, $blast_arg

    #----------------------------------------
    # Execute blast
    #----------------------------------------
    my $exec = "$CONFIG{BLAST} $blast_arg -query $seq_file -out $blast_file 2>&1";
    ##### |    Open: $exec
    open my $BLAST, '-|', $exec or confess "BLAST error: Can't execute $exec";
    my $output = do { local $/; <$BLAST> }; # join '', (<BLAST>);
    close $BLAST;
    #----------------------------------------

    #----------------------------------------
    # There shouldn't be any kind of output
    #----------------------------------------
    if (length($output)) {
        if ($CONFIG{DIEONBLASTERROR} and $output =~ /ERROR/i) {
            confess  "BLAST error: BLAST '$exec' returned $output";
        } else {
            warn "BLAST warn: BLAST '$exec' returned $output";
            ($output =~ /ERROR/i) and return 0;
        }
    }
    #----------------------------------------

    #----------------------------------------
    # Test if BLAST is complete
    # Check last line of $blast_file if it
    # ends correctly!
    my $line = `tail -n 2 $blast_file`;

    unless ($line =~ /^S2/m || $line =~ m|^</BlastOutput>$|m) {
        if ($CONFIG{DIEONBLASTERROR}) {
            confess  "BLAST error: BLAST '$exec' file ends with '$line'";
        } else {
            warn "BLAST error: BLAST '$exec' file ends with '$line'";
            return 0;
        }
    }
    #----------------------------------------

    return 1;
}
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------

=head2 parse_XML_BLAST

 Usage   : PPH::BLAST::parse_XML_BLAST($file, $HSP)
 Function: Parse a BLAST XML (m7 output) file to retrieve HSPs
           The problem is that we need not only the hsp but also the
           info from the hit (such as hit id) - therefore 2 arrays
 Args    : Filename
           ArrayRef (INOUT) of HSPs
 Returns : empty string if successful, otherwise error message

=cut

#------------------------------------------------------------------------------

sub parse_XML_BLAST {
    my ($file, $HSPr) = @_;
    #### PPH__BLAST__parse_XML_BLAST: $file

    scalar @$HSPr and warn "PPH::BLAST::parse_XML_BLAST: Non-empty INOUT parameter\n";

    return 'BLAST XML file missing' unless -e $file;
    return 'BLAST XML file empty' if -z $file;

    #----------------------------------------
    # Slowly go into Data structure
    # - sometimes it's either a HASH or ARRAY
    #----------------------------------------
    my $hits = XMLin($file);
    $hits    = $hits->{'BlastOutput_iterations'};
    $hits    = (ref $hits eq 'HASH') ? $hits : $hits->[0];
    $hits    = $hits->{'Iteration'};
    # IAA: Return immediately if no BLAST hits found
    return 'no hits found in BLAST XML'
      if exists $hits->{'Iteration_message'} && $hits->{'Iteration_message'} eq 'No hits found';
    $hits    = (ref $hits eq 'HASH') ? $hits : $hits->[0];
    $hits    = $hits->{'Iteration_hits'}->{'Hit'};
    # IAA: If we have just a single hit it will be a hash reference, convert it to array ref
    $hits    = [ $hits ] if ref $hits eq 'HASH';

    foreach my $hit (@$hits) {
      my $hsps = $hit->{'Hit_hsps'}->{'Hsp'};
      $hsps = (ref $hsps eq 'HASH') ? [ $hsps ] : $hsps;
      foreach my $hsp (@$hsps) {
        $hsp->{'Hsp_gaps'} =
            (exists $hsp->{'Hsp_gaps'}) ? $hsp->{'Hsp_gaps'} : 0;
        # IAA: Replace Selenocysteine (U) with Cysteine (C) and Pyrrolysine (O)
        # with Lysine (K) in the sequence since MAFFT as well as many other tools
        # do not recognize these non-standard residue codes.
        $hsp->{'Hsp_hseq'} =~ tr/UOuo/CKck/;
        # IAA: Modified to use PPH::Align::_calc_ident() which properly treats
        # ambiguous and unknown (or masked) aa codes both in query and hit sequences.
        $hsp->{'Identity'} = PPH::Align::_calc_ident($hsp->{'Hsp_qseq'}, $hsp->{'Hsp_hseq'});
        # IAA: This version does not cope well with ambiguous (B, J, Z) aa codes as well
        # as with X's in database hit sequences and hence should be avoided.
        my $segged = ($hsp->{'Hsp_qseq'} =~ tr/X//);
        $hsp->{'IdentitySEG'} =
          $hsp->{'Hsp_identity'} / ($hsp->{'Hsp_align-len'} - $hsp->{'Hsp_gaps'} - $segged);
        $hsp->{'Hit_num'}       = $hit->{'Hit_num'};
        $hsp->{'Hit_id'}        = $hit->{'Hit_id'};
        $hsp->{'Hit_accession'} = $hit->{'Hit_accession'};
        $hsp->{'Hit_def'}       = $hit->{'Hit_def'};
        push @$HSPr, $hsp;
      }
    }
    return scalar @$HSPr ? '' : 'no hits in BLAST XML';
}
#------------------------------------------------------------------------------

1;
