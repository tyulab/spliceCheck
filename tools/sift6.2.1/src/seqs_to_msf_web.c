/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */  
/* seqs_to_msf_web.c  
prints out for web style with references
reads in aligned sequences and prints out in msf format
names must be unique and < 20 characters already
*/

#define EXTERN
#include <assert.h>
#include "blocksprogs.h"
#include "blocklist.c"
#include "Alignment.c"
#include "List_Number.c"
#include "Protdist.c"

#define MAXSEQ 400 /* Maximum number of sequences */
#define LINE_LEN 800
#define FALSE 0
#define TRUE 1
#define MAXWIDTH 55 /* maximum block width */

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

/* Local routines */
void getargs (int argc, char* argv[], FILE** seqfp, 
		FILE** outfp);
void output_msf_web(FILE* fout, Block* block);

/* ################	MAIN 	###########################  */

int main
(int argc, char* argv[])
{

	FILE* seqfp; FILE* outfp; 
	Sequence *seqs[MAXSEQ];
	int nseqs, nnewseqs, aa_length, i;
	Block* block;
	int db_type; int seq_type;
	
	getargs (argc, argv, &seqfp, &outfp);

/*****READ SEQUENCES ************************************/
     /*-----------------------------------------------------------------*/
      /*   Check next for input file of sequences & assume are aligned */
      db_type = type_dbs(seqfp, DbInfo);
      /*   could set db_type = FLAT if it comes back negative   */
      seq_type = UNKNOWN_SEQ;
      seq_type = seq_type_dbs(seqfp, DbInfo, db_type, seq_type);
      if (seq_type == NA_SEQ)
      {
         fprintf(stderr, "WARNING: Sequences appear to be DNA but will be treated as protein.\n");
         seq_type = AA_SEQ;
      }
      rewind(seqfp);
      /*-----------------------------------------------------------------*/
      /*   read fasta sequences into memory                    */
      nseqs = 0;
      if (db_type >= 0)
      {
         while ( nseqs < MAXSEQ &&
             (seqs[nseqs] = read_a_sequence(seqfp, db_type, seq_type)) != NULL)
         {
            nseqs++;
         }
      }
	fclose (seqfp); 
	fix_names (nseqs, seqs); 
	change_Xes_to_dashes (nseqs, seqs);
	for (i = 0; i < nseqs; i++) {
		convert_sequence_with_beg_and_end_gaps_to_X (seqs[i]);
	}
	if (nseqs == MAXSEQ) 
	{  fprintf (stderr, "WARNING: Maximum number of sequences = %d\n", nseqs); }
	
		
	aa_length = get_length (seqs[0]);
	/* this is an alignment, all sequences should have same length */

        block = make_block (aa_length, 0, nseqs, seqs, FALSE);
	pb_weights (block);
	output_msf_web(outfp, block);
	
	fclose (outfp);
	free_seqs (seqs, nseqs);
	exit (0);

} /* end main */


void getargs (int argc, char* argv[], FILE** seqfp, 
                FILE** outfp)
{
	char seqfilename[LARGE_BUFF_LENGTH];
	char outfilename[LARGE_BUFF_LENGTH];
	
	if (argc < 2) 
	{
		printf ("seqs_to_msf \n");
		printf ("takes aligned sequences (fasta format) ");
		printf ("and prints it out in msf format");
	}

	if (argc > 1) strcpy (seqfilename, argv[1]);
	else
	{
		printf ("Enter name of file with  alignmens:\n");
		fgets (seqfilename, LARGE_BUFF_LENGTH, stdin);
	}

	if ((*seqfp = fopen (seqfilename, "r")) == NULL)
	{
		fprintf (stderr, "cannot open file %s \n", seqfilename);
		exit (-1);
	}


        if (argc > 2)  strcpy (outfilename, argv[2]);
        else
        {
                printf ("Enter name of out file\n");
                scanf ("%s", outfilename);
        }

        if ( (*outfp = fopen (outfilename, "w")) == NULL)
        {
                fprintf (stderr, "Cannot open file %s \n", outfilename);
                exit (-1);
        }
 
} /* end of getargs */


/*==========================================================================
Prints the alignment in MSF style
==========================================================================*/
void output_msf_web(FILE* fout, Block* block)
{
   int  nseq, seq, pos, curr_pos ;
   char outfile[80];
   char* error_string;
   char shortened_name[14];
   char* pstr;
   char* ptr;
   int space_index, space_length;
   int ncbi_ref;
 
   assert (block != NULL);
   nseq = block->num_sequences;

/*   fprintf(fout, "MSF of: %s\n", block->sequences[0].name);
   fprintf(fout, "MSF: %d Type: P Check: 5859 ..\n",
                                        block->width);
        for (seq = 0; seq < nseq; seq++)
        {
        	fprintf(fout, "Name: %20s Len: %5d Check: %5d Weight: %5d\n",
                block->sequences[seq].name, block->width,
                 750, (int) (block->sequences[seq].weight * 100));
         }
         fprintf(fout,"//\n");
 */        curr_pos = 0; pos = 0;
         while (curr_pos < block->width) {
            fprintf (fout, "\n");
            fprintf (fout, "                      %d  %50d\n", curr_pos+1, curr_pos + 50);
            for (seq=0; seq < nseq; seq++)
            {
		ptr = strstr (block->sequences[seq].name, "\n");
		if (ptr != NULL) {*ptr = '\0';}
                /* web stuff */
	        space_length = 20;
                space_length -= strlen(block->sequences[seq].name);
                for (space_index = 0; space_index < space_length; space_index++)
                {
                        fprintf (fout, " ");
                }
	ncbi_ref = FALSE;
	if (strstr (block->sequences[seq].name,"gi") != NULL) {
		ncbi_ref = TRUE;
	} /*else if (strstr (block->sequences[seq].name,"gb") != NULL) {
		ncbi_ref = TRUE;
	} else if (strstr (block->sequences[seq].name, "gln") != NULL) {
		ncbi_ref = TRUE;
	} else if (strstr (block->sequences[seq].name, "emb") != NULL) {
		ncbi_ref = TRUE;
	} */

	
	if (seq != 0 && ncbi_ref == FALSE ) {
	/* SWISS-PROT ref */
		fprintf (fout, "<A HREF=\"http://www.expasy.ch/cgi-bin/get-sprot-entry?%s\">", block->sequences[seq].name);
	} else if (seq != 0 && ncbi_ref == TRUE) {
		fprintf (fout, "<A HREF=\"http://www.ncbi.nlm.nih.gov:80/entrez/query.fcgi?cmd=Retrieve&db=Protein&list_uids=%s&dopt=GenPept\">",
					&(block->sequences[seq].name[2]) );

	}
		fprintf(fout, "%s", block->sequences[seq].name);
                if (seq != 0) {
			fprintf (fout, "</A>");
		}
		fprintf (fout, " ");
		for (pos = curr_pos; (pos < block->width)
                                && (pos < curr_pos + 50 ); pos++)
                {
                  if ( pos % 10 == 0) {
                        fprintf (fout, " ");
                  }
                  fprintf(fout, "%c",
                                aa_btoa[block->sequences[seq].sequence[pos]]);
                } /* end of for pos */
                fprintf (fout, "\n");
            } /* end of for seq*/
            curr_pos += 50;
        } /* end of while */

} /* end of output_msf */


