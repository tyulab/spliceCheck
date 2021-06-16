/* COPYRIGHT 2000 , Fred Hutchinson Cancer Research Center */ 

#ifndef _CLUMPING_H_
#define _CLUMPING_H_

#include <sequences.h>

#define INDEX(n, col, row) (col*n - (col*(col+3))/2 - 1 + row)
#define MAXSEQ 400

struct ClusterOfSequences {
	Sequence* seqs[MAXSEQ];
	int nseqs;
}; 


void free_ClusterOfSequences (struct ClusterOfSequences* cluster, 
				int num_clusters);



struct ClusterOfSequences*
cluster_seqs (double clus, Sequence* seqs[MAXSEQ], int nseq, int* numclusters);

Block**
make_blocks_from_clumps (Sequence* seqs[MAXSEQ], 
			int nseqs, 		
			double clumping_threshold, 
			 int* num_blocks);

void
clump_into_consensus_seqs (Sequence* seqs[MAXSEQ],
                        int nseqs,              
                        double clumping_threshold, 
                         int* num_consensus_seqs,
			Sequence* consensus_seqs[MAXSEQ]);


int
print_consensus_seqs_aligned  (double clus, Sequence* seqs[MAXSEQ], 
				int nseq, FILE* outfp,
			        FILE* keyoutfp, int option);

#endif
