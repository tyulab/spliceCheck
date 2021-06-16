#!/bin/csh
#		SIFT.csh
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

# Argument 1: NCBI protein gi #. Just the number, nothing else.
# Argument 2: The substitution file (file containing amino acid substitutions to be predicted on.
# Argument 3: Type of sequences to get from NCBI BLink. "BEST" or "ALL"
# BEST is recommended if high-confidence scores can be obtained

### Set these for your installation
#	Location of psiblast (absolute path)
setenv NCBI /home/ubuntu/workspace/final/tools/ncbi-blast-2.7.1+/bin
#setenv NCBI /usr/local/blast/

#	Location of SIFT (absolute path)
setenv SIFT_DIR /home/ubuntu/workspace/final/tools/sift6.2.1

#	SIFT's output files are written here
set tmpdir = $SIFT_DIR/tmp

#       Location of BLIMPS. This should not need any editing
setenv BLIMPS_DIR $SIFT_DIR/blimps

### Shouldn't need to make any more changes, look for output in $tmpdir
set bindir = $SIFT_DIR/bin
set gi = $1
set polymorphism_file = $2
set outfile = $tmpdir/$1.SIFTprediction
set bestORall = $3 # if set to BEST, takes best reciprocal hits from BLink,
			# otherwise it will take all hits from BLINK

# get unaligned sequences from NCBI
set unalignedseqs = $tmpdir/$gi.unaligned
perl $bindir/perlscripts/get_BLINK_seq.pl $gi $tmpdir/$gi.unaligned $bestORall

# getting the alignment from PSI-BLAST, quick and dirty
 cat $unalignedseqs | perl -pe 's/\|//' | perl -pe 's/\|/\t/' | perl -pe 's/^\n//' | awk '{if ($1 ~ /^>/) { print $1} else { print;}}' > $unalignedseqs.2
        set alignedseqs = $unalignedseqs.aligned
        perl $bindir/perlscripts/separate_query_from_rest_of_seqs.pl $unalignedseqs.2  $unalignedseqs.queryseq  $unalignedseqs.database
        $NCBI/makeblastdb -in $unalignedseqs.database -dbtype prot
# extremely large evalues and multipass threshold because want to make sure get
#all the sequences the user submits
        $NCBI/psiblast -db $unalignedseqs.database -query $unalignedseqs.queryseq -out $unalignedseqs.psiblastout -outfmt 0 -evalue 10 -num_descriptions 399 -num_alignments 399
        echo QUERY > $unalignedseqs.listseq
        grep ">" $unalignedseqs.database | cut -d" " -f1 | cut -c2- >> $unalignedseqs.listseq
        $bindir/seqs_from_psiblast_res $unalignedseqs.psiblastout $unalignedseqs.listseq 1 $unalignedseqs.queryseq $alignedseqs $$


# final scoring
$bindir/info_on_seqs $alignedseqs $polymorphism_file $outfile

exit

