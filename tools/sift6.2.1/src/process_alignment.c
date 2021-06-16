/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */  

/* process_alignment.c 
07-12-00 reads in clustal gapped alignment, removes columns that
correspond to gaps from the query sequence, prints X's at beg and
end of sequence

so an alignment is still printed out (in fasta format) 


*/

#define EXTERN
#include <assert.h>
#include "blocksprogs.h"
#include "blocklist.c"
#include "Matrix_Info.c"
#include "Alignment.c"
#include "PN_blocks.c"
#include "PN_convert.c"
#define MAXSEQ 400 /* Maximum number of sequences */

#define LINE_LEN 800
#define FALSE 0
#define TRUE 1
#define MAXWIDTH 55 /* maximum block width */

/* Local routines */
FILE* errorfp;

void getargs (int argc, char* argv[], FILE** psiblastfp, FILE** outfp);

/* MAIN */

int main
(int argc, char* argv[])
{

	FILE* alignmentfp; FILE* outfp;
	int max_iterations;
	char desc[SMALL_BUFF_LENGTH];
	int converged;
	Sequence *seqs[MAXSEQ];
	Sequence *newseqs[MAXSEQ];
	int nseqs, aa_length, i, last_loc;
	int db_type, seq_type;
	char query_name[SMALL_BUFF_LENGTH];
	char line[LARGE_BUFF_LENGTH]; int done;
	long fp_pos;

	getargs (argc, argv, &alignmentfp, &outfp);


	nseqs = 0;

	db_type = type_dbs (alignmentfp, DbInfo);
	seq_type = UNKNOWN_SEQ;
	seq_type = seq_type_dbs (alignmentfp, DbInfo, db_type, seq_type);
	if (seq_type == NA_SEQ) {
		fprintf (stderr, "WARNING: Sequences appear to be DNA but will be treated as protein\n");
		seq_type = AA_SEQ;
	}
	rewind (alignmentfp);
	if (db_type >= 0  ) {
		while ( nseqs < MAXSEQ && (seqs[nseqs] =
		 read_a_sequence (alignmentfp, db_type, seq_type)) != NULL) 
		{
	printf ("Trying fasta format %d\n", nseqs);
			nseqs++;
		}
	}
	else 
	{	
		/* skip extra lines coming from html stuff */
		done = 0; fp_pos = 0;
		while ( !feof (alignmentfp) && 
			fgets (line, LARGE_BUFF_LENGTH, alignmentfp) != NULL 
			&& !done) {
			if (strncmp (line, "//", 2) == 0 || 
				strncmp (line, "CLUSTAL", 7) == 0) {
				done = 1;
			} else {
				fp_pos = ftell (alignmentfp);
			}
		} /* end skip newlines and meme format for html stuff */
		fseek (alignmentfp, fp_pos, 0);

		nseqs = try_clustal (alignmentfp, seqs, desc);
		printf ("tried clustal %d\n", nseqs);
		if (nseqs <= 0) 
		{
			nseqs = try_msf (alignmentfp, seqs, desc);
			printf ("tried msf %d\n", nseqs);
       		 	change_periods_to_dashes (nseqs, seqs);
		}
	} /* end try clustal */
	fix_names (nseqs, seqs);
	if (nseqs == MAXSEQ) 
	{  printf ("WARNING: Maximum number of sequences = %d\n", nseqs); }
/*	make_query_first (seqs, nseqs); */


	/* have to rename first sequence QUERY so remove-gaps can operate*/
	strcpy (query_name, seqs[0]->name);
	strcpy (seqs[0]->name, "QUERY"); 

	aa_length = get_length (seqs[0]);
	remove_gaps (newseqs, seqs, nseqs, aa_length);  

	
	strcpy ( newseqs[0]->name, query_name);
	for (i=0; i<nseqs; i++) {
		convert_sequence_with_beg_and_end_gaps_to_X (newseqs[i]);
		output_sequence(newseqs[i], outfp); 
	/*	output_sequence(seqs[i], outfp); */
	}
	free_seqs (newseqs, nseqs);  
	fclose (alignmentfp);
	fclose (outfp);
	return 0;
} /* end main */

/*=====================================================================*/
/* sequences over length threshold */
/* return sequences (in newseqs ) that have at least seqs[0]->length* threshold
in alignment .  will provide seed for starting blocks. returns the number
of sequences in newseq*/

int
sequences_over_length_threshold
(Sequence* newseqs[MAXSEQ], Sequence* seqs[MAXSEQ], int nseqs, double threshold)
{
	int newseq_index, i, length, query_length;

	query_length = get_length (seqs[0]);
	newseq_index = 0;
	for (i = 0; i < nseqs; i++) {
		length = get_length (seqs[i]);
		if (length >= threshold * query_length) {
			/*newseqs[newseq_index] = copy_sequence (seqs[i]); */
			newseq_index++;
		}
	}
	return newseq_index;
} /* end of sequences_over_length_threshold */


void getargs (int argc, char* argv[], FILE** psiblastfp, FILE** outfp)
{
	char psiblastfilename[LARGE_BUFF_LENGTH];
	char outfilename[LARGE_BUFF_LENGTH];
	
	if (argc < 2) 
	{
		printf ("process alignment \n");
		
	}

	if (argc > 1) strcpy (psiblastfilename, argv[1]);
	else
	{
		printf ("Enter name of psiblast outfile with alignment:\n");
		fgets (psiblastfilename, LARGE_BUFF_LENGTH, stdin);
	}

	if ((*psiblastfp = fopen (psiblastfilename, "r")) == NULL)
	{
		printf ("cannot open file %s \n", psiblastfilename);
		exit (-1);
	}

	if (argc > 2) strcpy (outfilename, argv[2]);
	else 
	{
		printf ("Enter name of outfile\n");
		fgets (outfilename, LARGE_BUFF_LENGTH, stdin);
	}

	if ( (*outfp = fopen (outfilename, "w")) == NULL)
	{
		printf ("Cannot open file %s \n", outfilename);
		exit (-1);
	}

} /* end of getargs */

