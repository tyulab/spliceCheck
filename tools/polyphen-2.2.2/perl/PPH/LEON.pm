package PPH::LEON;

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD

=head1 NAME

PPH::Leon

=head1 DESCRIPTION

Wrapper for LEON - multiple aLignment Evaluation Of Neighbours

I experienced sometimes programs (query2relacs) which went rogue - I hope that
this implementation controls this issue

 Julie D. Thompson, Véronique Prigent and Olivier Poch*
 Nucleic Acids Research, 2004, Vol. 32, No. 4 1298-1307

=head1 SUBVERSION

 $LastChangedDate: 2011-11-21 13:28:37 -0500 (Mon, 21 Nov 2011) $
 $LastChangedRevision: 376 $
 $LastChangedBy: ivan $


=cut

#DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
use Carp qw(cluck);
use strict;
use warnings;
use File::Basename;
use File::Copy;

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

use PPH::Data;

#RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR

#------------------------------------------------------------------------------

=head1 ROUTINES

=head2 run_LEON

  Usage   : run_LEON(fas_infile, query_id, fas_outfile)
  Function: evaluates a FASTA formatted alignment using LEON
  Args    : infile  FASTA formatted infile, e.g.: /tmp/PPH_AiaY/P04798_msa_ou
            qname   Query sequence name, e.g.: QUERY
            outfile FASTA formatted outfile, e.g.: /tmp/PPH_AiaY/P04798_leon___ou2.tfa
  Returns : "OK" if successful

=cut

#----------------------------------------
sub run_LEON {
    my ($infile, $qname, $outfile) = @_;

    (-e $infile)
        or do {
            warn "PPH::LEON::run_LEON $infile does'n exist";
            return 'LEON INFILE wrong'
        };

    my $input_basename = fileparse($infile);

    #### PPH__LEON__run_LEON: $infile, $qname, $outfile
    # '_' is  important since cluspack is chopping off the postfix push
    # e.g.: P04798_msa_ou_leon
    my $out = "${input_basename}_leon";
    $out =~ s/\./_/;

    #----------------------------------------
    # LEON has some problems with long ids..
    #----------------------------------------
    open my $F, '<', $infile  or die "Can't open $infile to read";
    $infile = "$out.infile";
    open my $W, '>', $infile  or die "Can't write to $infile";
    my $tmp = do { local $/; <$F> };
    close $F;
    $tmp =~ s/(^>\S*).*/$1/mg;
    print $W $tmp;
    close $W;

    #----------------------------------------

    my @PIDS;
    my $pid;
    eval {
        local $SIG{'ALRM'} = sub { die 'TIMEOUT'; } if $CONFIG{'LEON_T_LIMIT'};
        alarm $CONFIG{'LEON_T_LIMIT'} if $CONFIG{'LEON_T_LIMIT'};

        #----------------------------------------
        # Convert input alignment
        #----------------------------------------
        ##### |    convseq: $infile
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}convseq $infile $out.msf gcg 2>>$CONFIG{LIMBO}")
                or die "CONVSEQ (1) failed $?\n";
        }
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}convseq $infile $out.rsf rsf 2>>$CONFIG{LIMBO}")
                or die "CONVSEQ (2) failed $?\n";
        }
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}convseq $infile $out.tfa fasta 2>>$CONFIG{LIMBO}")
                or die "CONVSEQ (3) failed $?\n";

        }
        ##### |    convseq done
        #----------------------------------------

        #----------------------------------------
        # Run Cluspack (SECATOR)
        #----------------------------------------
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            my $cmd =
                sprintf("$CONFIG{CLUSPACK} $CONFIG{CLUSPACK_ARG}", "$out.msf").
                    " >>$CONFIG{LIMBO} 2>&1";
            ##### |    cluspack: $cmd
            exec("$cmd") or die "CLUSPACK failed $?\n";
        }
        ##### |    cluspack done
        #----------------------------------------

        #----------------------------------------
        # Run COILS
        #----------------------------------------
        ##### |    coils (.tfa): $out
        $ENV{COILSDIR}=$CONFIG{MATRIX};
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{COILS} < $out.tfa > $out.coils 2>>$CONFIG{LIMBO}")
                or die "COILS failed: $?\n";
        }
        ##### |    coils done
        #----------------------------------------

        #----------------------------------------
        # TM / LOW COMPLEXITY REGIONS
        #----------------------------------------
        ##### |    TM & LCR regions (.tfa): $out
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
              exec("$CONFIG{LEON}resbias $out.tfa -v > $out.bias 2>>$CONFIG{LIMBO}")
                or die "TM & LCR (1) failed $?\n";
        }
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}bias2relacs $out.tfa $out.bias $out.clu $out.int1.rsf 2>>$CONFIG{LIMBO}")
                or die "TM & LCR (2) failed $?\n";
        }
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}coils2relacs $out.int1.rsf $out.coils $out.clu $out.int2.rsf 2>>$CONFIG{LIMBO}")
                or die "TM & LCR (3) failed $?\n";
        }
        ##### |    resbias, bias2relacs, coils2relacs done
        #----------------------------------------

        #----------------------------------------
        # Calculate core blocks with RASCAL
        #----------------------------------------
        ##### |    RASCAL (.tfa, .clu): $out
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}comp_query $out.tfa $out.clu $qname > $out.cmp 2>>$CONFIG{LIMBO}")
                or die "RASCAL (1) failed $?\n";
        }
        (-e "$out.cmp") or die "RASCAL (1) failed: $out.cmp missing\n";

        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}query_seqerrs $out.tfa $out.clu > $out.seqerrs 2>>$CONFIG{LIMBO}")
                or die "RASCAL (2) failed $?\n";
        }
        (-e "$out.seqerrs") or die "RASCAL (2) failed: $out.seqerrs missing\n";

        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}query2relacs $out.int2.rsf $out.cmp $out.seqerrs $qname $out.int3.rsf 2>>$CONFIG{LIMBO}")
                or die "RASCAL (3) failed $?\n";
        }
        (-e "$out.int3.rsf") or die "RASCAL (3) failed: $out.int3.rsf missing\n";
        #----------------------------------------

        #----------------------------------------
        # chain core blocks into conserved regions
        # set threshold parameters x=5, d=40, T=280, L=21
        #----------------------------------------
        ##### |    CORE BLOCKS (.int3.rsf, .int4.rfs): $out
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}chain_blocks $out.int3.rsf $qname $out.int4.rsf 5 40 280 21 >$out.log 2>&1")
                or die "CORE BLOCKS failed $?\n";
        }
        #----------------------------------------

        #----------------------------------------
        # output in FASTA format
        #----------------------------------------
        ##### |    CONVSEQ2 (.int4.rsf): $out
        # e.g.: P04798_leon___ou2.tfa
        if ($pid = fork) {
            push @PIDS, $pid; waitpid($pid, 0);
        } else {
            exec("$CONFIG{LEON}convseq $out.int4.rsf $outfile fasta 2>>$CONFIG{LIMBO}")
                or die "CONVSEQ (4) failed $?\n";
        }
        #----------------------------------------
        alarm 0 if $CONFIG{'LEON_T_LIMIT'};
    };

    if ($@ =~ /^TIMEOUT/) {
        system("kill -9 @PIDS >>$CONFIG{LIMBO} 2>&1");
        return 'LEON TIMEOUT';
    } elsif ($@) {
        system("kill -9 @PIDS >>$CONFIG{LIMBO} 2>&1");
        return "LEON ERROR: $@";
    }
    #----------------------------------------

    return 'OK';
}
#----------------------------------------

#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
1;
