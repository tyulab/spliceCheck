#!/bin/csh
#               SIFT.csh
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software
#

# Argument 1: the protein sequence file in fasta format
# Argument 2: the pathname to the protein sequence database
# Argument 3: the substitution file (file containing amino acid substitutions to be predicted on.

### Set these for your installation
#       Location of blastpgp
setenv NCBI /home/ubuntu/workspace/final/tools/ncbi-blast-2.7.1+/bin


#       Location of SIFT
setenv SIFT_DIR /home/ubuntu/workspace/final/tools/sift6.2.1

#       SIFT's output files are written here
setenv tmpdir $SIFT_DIR/tmp/

#       Location of BLIMPS.This should not need any editing.
setenv BLIMPS_DIR $SIFT_DIR/blimps

### Shouldn't need to make any more changes, look for output in $tmpsift
set bindir = $SIFT_DIR/bin
set root_file = $1:r
set tail_of_root = $root_file:t
set tmpfasta = $tmpdir/$tail_of_root.alignedfasta
set tmpsift = $tmpdir/$tail_of_root.SIFTprediction

$bindir/seqs_chosen_via_median_info.csh $1 $2 2.75
$bindir/info_on_seqs $tmpfasta $3 $tmpsift
echo "Output in $tmpsift"

exit
