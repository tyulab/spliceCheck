/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */  
/* seqs_to_matrixweb.c
converts sequences (fasta or clustal format) to a matrix for HTML format
given outfilename, generates 2 files, outfilename.matrix
				      outfilename.predictions
01/02/01 took out error annotation if sequence is not allowed because SIFT_alts
call already takes care of it

07/04/01 version 2 to print information of sequence
	change output_predictions so it prints out median information

07/06/01 replaced old seqs_to_matrixweb -- 
	(original version in seqs_to_matrixweb0301 

11/01/01 can pass in argument for sequences with $8 % identity  
	will be filtered out

May 7, 2001 returns predictions for multiple substituitons at a postion. doesn't work -- too slow. array implementation better . try debugging.
 8/25/02 generating too many errors, added ErrorLevelReport.
returned to returning prediction for only one subtstitution per position

*/

#define EXTERN
#include <assert.h>
#include "blocksprogs.h"
#include "Alignment.c"
#include "PN_blocks.c"
#include "Matrix_Info.c"
#include "Psiblast.c"

#define MAXSEQ 400 /* Maximum number of sequences */
#define LINE_LEN 800
#define FALSE 0
#define TRUE 1
#define ADEQUATE_SEQ_INFO 3.25

/* Local routines */

void getargs (int argc, char* argv[], FILE** seqfp, FILE** polymorphfp, 
					char outfilename[LARGE_BUFF_LENGTH], 
				int* scale, double* threshold,
			int* gap_option,int* exp_option, int* seq_identity);
void output_matrix_sPN_web (Matrix* matrix, FILE* omfp, double threshold, Sequence* query_seq, double* fraction_stored);

void
output_predictions (Residue* polymorph_data, FILE* outfp, Matrix* pssm,
                   double threshold);

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

/* MAIN */

int main
(int argc, char* argv[])
{

	FILE* seqfp; FILE* outfp; FILE* polymorphismfp; 
	char outfilename[LARGE_BUFF_LENGTH], currentoutfile[LARGE_BUFF_LENGTH];
	char tempname[LARGE_BUFF_LENGTH];
	char desc[SMALL_BUFF_LENGTH]; int desc_length; char* strptr;
	Sequence *seqs[MAXSEQ];
	int nseqs, aa_length, i;
	Block* block;
	int db_type;
	int seq_type;
	Matrix* pssm;
	double threshold;
	Residue* polymorph_data;
	int subtract_option, diri_option;
	int gap_option, exp_option;
	int scale;
	int seq_identity;

	double* fraction_stored; /* basic aa/total number of seq at each pos */

	ErrorLevelReport = 5;

	init_frq_qij();
	printf ( "tell me i've entered\n");
	getargs (argc, argv, &seqfp, &polymorphismfp, outfilename, 
		&scale, &threshold, &gap_option, &exp_option, &seq_identity);
	nseqs = 0;
	subtract_option = FALSE; diri_option = TRUE;
     /*-----------------------------------------------------------------*/
      /*   Check next for input file of sequences & assume are aligned */
      db_type = type_dbs(seqfp, DbInfo);
      /*   could set db_type = FLAT if it comes back negative   */
      seq_type = UNKNOWN_SEQ;
      seq_type = seq_type_dbs(seqfp, DbInfo, db_type, seq_type);
      if (seq_type == NA_SEQ)
      {
         fprintf(stderr, "WARNING: Sequences appear to be DNA but will be treated as protein\n");
         seq_type = AA_SEQ;
      }
      rewind(seqfp);
      /*-----------------------------------------------------------------*/
      /*   read fasta sequences into memory                    */
      if (db_type >= 0)
      {
         while ( nseqs < MAXSEQ &&
             (seqs[nseqs] = read_a_sequence(seqfp, db_type, seq_type)) != NULL)
         {
            nseqs++;
         }
      }
      else /*  CLUSTAL or MSF? */
      {
         /* get tail of outfilename, will use this as a description to make
	a mablock file */
	strcpy (tempname, outfilename);
	strptr = strtok (tempname, "/");
	while ((strptr = strtok ((char*) NULL, "/ \t\n\r\0")) != NULL) {
		desc_length = strlen(strptr);
		if (desc_length > SMALL_BUFF_LENGTH) { 
			desc_length = SMALL_BUFF_LENGTH;
		}
		strncpy (desc, strptr, desc_length);
		desc[desc_length] = '\0';
	}	 
	 nseqs = try_clustal(seqfp, seqs, desc);
         if (nseqs <= 0)
         { nseqs = try_msf(seqfp, seqs, desc); 
	   change_periods_to_dashes(nseqs, seqs); 
          }
          for (i = 0; i < nseqs; i++) {
                convert_sequence_with_beg_and_end_gaps_to_X (seqs[i]);
          } 
	}
      fix_names(nseqs, seqs);
	if (nseqs == 0) {
		fprintf (errorfp, "Cannot read sequences.  Check format.\n" );
		exit (-1);
	}
	if (nseqs == MAXSEQ) 
	{  fprintf (stderr, "WARNING: Maximum number of sequences = %d\n", nseqs); }
        /* reduce redundancy with query sequence 01/03/00 
           so sequence weighting doesn't include more than 1 sequence
           representing query */
        remove_seqs_percent_identical_to_query (seqs, &nseqs, 
						(double) seq_identity);
	
	aa_length = get_length (seqs[0]);
	
	/* this is an alignment, all sequences should have same length */

	/* READ POLYMORPHISMS */
printf ("about to read polymrophsims\n");
	polymorph_data = read_polymorphism_data (polymorphismfp, seqs[0]);
	fclose (polymorphismfp);

	block = make_block (aa_length, 0, nseqs, seqs, FALSE);
printf ("made the block\n");
        fraction_stored = calculate_basic_aa_fraction (block);
	pssm = SIFT_prediction (block, diri_option, gap_option, exp_option, subtract_option);
	printf ("ok, how's it going");

	strcpy (currentoutfile, outfilename);
	strcat (currentoutfile, ".matrix");
        if ( (outfp = fopen (currentoutfile, "w")) == NULL)
        {
                fprintf (errorfp, "Cannot open matrix outfile %s \n", outfilename);
                exit (-1);
        }
	output_matrix_sPN_web (pssm, outfp, threshold, seqs[0], fraction_stored);
	fclose (outfp);
	
        strcpy (currentoutfile, outfilename);
        strcat (currentoutfile, ".predictions");
        if ( (outfp = fopen (currentoutfile, "w")) == NULL)
        {
                fprintf (errorfp, "Cannot open predictions file %s \n", outfilename);
                exit (-1);
        }
	output_predictions (polymorph_data, outfp, pssm, threshold); 
	fclose (outfp);
	
	free (polymorph_data);
	free_block (block);
        free_matrix (pssm);
	fclose (seqfp); 
	free_seqs (seqs, nseqs);
	fclose (errorfp);
	rm_file (errorfilename);
	exit (0);

} /* end main */

void 
output_predictions ( Residue* polymorph_data, FILE* outfp, Matrix* pssm,
		   double threshold)
{
	Residue substitution, original_aa;
	int pos;
	int substitution_exists;
	double median;
	Block* block_with_seqs_at_pos;
	double* info_array;
	AAnode current;

	substitution_exists = FALSE;

	for (pos = 0; pos < pssm->width; pos++) {
		original_aa = pssm->block->residues[0][pos];
		if (pssm->weights[original_aa][pos] < threshold) {
			fprintf (stderr, "<font color=red>WARNING</font>:  Original amino acid %c at position %d is not allowed by the prediction. <BR><BR>\n", aa_btoa[original_aa], pos + 1);	
		
		}  
		substitution = polymorph_data[pos];
		if (substitution != aa_atob['-'])   { 
		/* at least one substitution exists */
			substitution_exists = TRUE;
/* calculate median information */
                block_with_seqs_at_pos = subblock_of_seqs_with_aa_at_pos
                                                        (pssm->block, pos);
                info_array =  calculate_info_for_each_pos
                                (block_with_seqs_at_pos, FALSE);
                median = median_of_array (info_array, pssm->block->width);


			if (pssm->weights[substitution][pos] >= threshold) {
				fprintf (outfp, "Substitution at pos %d from %c to %c is predicted to be TOLERATED with a score of %.2f.\n", pos + 1, aa_btoa[original_aa], aa_btoa[substitution], pssm->weights[substitution][pos]);
			} else { 
				fprintf (outfp, "Substitution at pos %d from %c to %c is predicted to <font color=red>AFFECT PROTEIN FUNCTION</font> with a score of %.2f.\n", pos + 1, aa_btoa[original_aa], aa_btoa[substitution], pssm->weights[substitution][pos]);
			}
		fprintf (outfp, "    Median sequence conservation: %.2f\n", 
								median);
		fprintf (outfp, "    Sequences represented at this position:%d\n", 
		block_with_seqs_at_pos->num_sequences);	

		if (median > ADEQUATE_SEQ_INFO
			 && (pssm->weights[substitution][pos] < threshold)) {
			fprintf (outfp, "<font color=red>WARNING!!</font>  This substitution may have been predicted to affect function just because \n the sequences used were not diverse enough.  <b>There is LOW CONFIDENCE in this prediction.</b>\n");
		} /* end if median > ADEQUATE_SEQ_INFO & was predicted 
					to be deleterious */
		fprintf (outfp, "\n");

		free (info_array); free_block (block_with_seqs_at_pos);
		} /* end of if there is a substitution to predict on */
	} /* end of for pos */


	if (substitution_exists == FALSE) {
		fprintf (outfp, "No substitutions were submitted.<BR>\n");
	}

} /* end of output_predictions procedure */

int
number_of_digits (int number)
{
	if (number < 10) {
		return 1;
	} else if (number < 100) {
		return 2;
	} else if (number < 1000) {
		return 3;
	} else if (number < 10000) {
		return 4;
	}
	return 5;
}

void
output_matrix_sPN_web (Matrix* matrix, FILE* omfp, double threshold, 
			Sequence* query_seq, double* fraction_stored)
{

	char c; int l;
	int dig_number;
	int pos;
    fprintf (omfp, "Each row corresponds to a position in the reference protein.  ");
    fprintf (omfp, "Below each position is the fraction of sequences that contain one of the basic amino acids.  A low fraction indicates the position is either severely gapped or unalignable and has little information.  Expect poor prediction at these positions.<BR>");
	fprintf (omfp, "Each column corresponds to one of the twenty amino acids. <BR>");
    fprintf (omfp, "Each entry contains the score at a particular position (row) for ");
    fprintf (omfp, "an amino acid substitution (column).  Substitutions predicted to"
);
    fprintf (omfp, " be intolerant are highlighted in red.<BR><BR>\n");
    fprintf (omfp, "<BR>\n");
    fprintf (omfp, "<table cellspacing=0 border=0 width=0 cols=21>\n");

	for (l=0; l<matrix->width; l++) {

	/* every 25 positions, print out amino acids for easier reading*/
	fprintf (omfp, "<tr>");
	if (fmod( (double) l, 25.0) == 0) {
	    fprintf (omfp, "<th>pos</th>\n");
	    for (c='A'; c < 'Z'; c++) {
      		if ((c != 'J') && (c != 'O') && (c != 'U') && 
					(c != 'B') && (c != 'X') ) {
		        fprintf(omfp, "<th>%c</th>\n", c);
      		}
             }
	fprintf (omfp, "</tr>\n");
	}
	pos = l + matrix->block->sequences[0].position;
	fprintf (omfp, "<tr><th>%2d%c %.2f</th>\n", pos, 
					aa_btoa[query_seq->sequence[l]], 
					fraction_stored[l]);
      for (c='A'; c < 'Z'; c++) {
        if ((c != 'J') && (c != 'O') && (c != 'U') && (c != 'B') && (c != 'X') )
	 {
          fprintf (omfp, "<td>\n");
	  if (matrix->weights[aa_atob[(int)c]][l] < threshold) {
		fprintf (omfp, "<font color=red>%.2f</font>\n", matrix->weights[aa_atob[(int)c]][l]);
	   } else {
		fprintf (omfp, "%.2f", matrix->weights[aa_atob[(int)c]][l]);
	   }
	fprintf (omfp, "</td>\n");
	} /* end of if not amino acid characters*/
      } /* end of for alphabet*/
      fprintf (omfp, "</tr>\n");
    }
   fprintf (omfp, "</table>\n");

} /* end of output_matrix */


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


void getargs (int argc, char* argv[], FILE** seqfp, FILE** polymorphfp,
			char outfilename[LARGE_BUFF_LENGTH], 
	int* scale,double* threshold, int* gap_option, int* exp_option,
	int *seq_identity)
{
	char seqfilename[LINE_LEN];
	char seq_outfilename[LINE_LEN];
	char substfilename[LARGE_BUFF_LENGTH];
	
	if (argc < 5) 
	{
		printf ("seqs_to_matrixweb.c :  Converts alignment of sequences to matrix for web\n");
		
	}

	if (argc > 1) strcpy (seqfilename, argv[1]);
	else
	{
		printf ("Enter filename with sequences:\n");
		fgets (seqfilename, LINE_LEN, stdin);
	}
printf ("fawegwa\n");
	if ((*seqfp = fopen (seqfilename, "r")) == NULL)
	{
		printf ("cannot open file %s \n", seqfilename);
		exit (-1);
	}
printf ("eaegrtjkl\n");
	if (argc > 2) strcpy (substfilename, argv[2]);
	else {
		printf ("Enter file with substitutions.\n");
		printf ("format:\n");
		printf ("M1Y\n");
		printf ("S2K\n");
		printf ("...\n\n");
		fgets (substfilename, LARGE_BUFF_LENGTH, stdin);
	}
        if (substfilename[0] == '-') { *polymorphfp = NULL; }
	else if ((*polymorphfp = fopen (substfilename, "r")) == NULL)
        {
                printf ("cannot open file %s \n", substfilename);
                exit (-1);
        }


	if (argc > 3) strcpy (outfilename, argv[3]);
	else 
	{
		printf ("Enter name of outfile\n");
		fgets (outfilename, LARGE_BUFF_LENGTH, stdin);
	}

        strcpy (errorfilename, outfilename);
        strcat (errorfilename, ".error");
        if ((errorfp = fopen (errorfilename, "w")) == NULL) {
                printf ("couldn't open file %s\n", errorfilename);
                exit (-1);
        }


	if (argc > 4) *scale=atoi(argv[4]);
        else
        {
                printf ("Enter conversion scale for matrix\n");
                scanf ("%d", scale);
        }


	if (argc > 5) *threshold=atof(argv[5]);
	else
	{
		printf ("Enter threshold for intolerance\n");
		scanf ("%lf", threshold);
	}

	if (argc > 6) *gap_option = atoi (argv[6]);
	else {
		printf ("Enter gap option option 1 for allow everything\n");
                printf ("                        0 for ignore gaps\n");
                scanf ("%d", gap_option);
        }
	if (argc > 7) *exp_option = atoi (argv[7]);
	else {
		printf ("Enter exp option 0 - # diffaase\n");
                printf ("                1 - similarity scale\n");
                scanf ("%d", exp_option);
        }

	if (argc > 8) *seq_identity = atoi (argv[8]);
	else {
		printf ("Enter integer 0-100 of sequence identity that will be filtered out\n");
		*seq_identity = 100; /* default of 100% */
		/*scanf ("%d", seq_identity); */
	}
 
} /* end of getargs */

