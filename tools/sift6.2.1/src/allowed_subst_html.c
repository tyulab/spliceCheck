/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

/* allowed_subst_ps.c

prints html file of amino acids not tolerated and tolerated in a table

10/23/00 threshold parameter doesn't work, just keeping it for future,
	uses INTOLERANCE_PROB_THRESHOLD
12/08/04  added html before and after 
*/
#define EXTERN
#include "blocksprogs.h"
#include <string.h>
#include "Matrix_Info.c"
#include "PN_convert.c"
#define INTOLERANCE_PROB_THRESHOLD .05
#include "Psiblast.c"
#include "PN_blocks.c"

struct html_table {
	char* table[MATRIX_AA_WIDTH]; /* table of allowed subsitutions */
	char* intolerant_table[MATRIX_AA_WIDTH]; /* table of intolerant subst*/
	int* allowed_subst; /* # of allowed substitutions at allowed_subt[pos]*/
		     /* used to determine spacing when printing out table */
	int length;
	Sequence* ref_seq; /* original sequence */
};

typedef struct html_table HtmlTable;

HtmlTable* initialize_table (int size);
void getargs (int argc, char* argv[], char outfilename[LARGE_BUFF_LENGTH],
		FILE** seqfp, double* threshold,int* gap_option,int* exp_option,
		int* seq_identity);
void copy_aa_scores (struct aa_score list[20],
                Matrix* no_pseudo_pssm, Matrix* pseudo_pssm, int pos);

void make_coordinates (int pos_shift,
                        Matrix* no_pseudo_pssm, Matrix* pseudo_pssm, 
			HtmlTable* html_table);
void make_coordinate (FILE* outfp, int pos, struct aa_score list[20]); 

void
print_html_table (FILE* fout, HtmlTable* html_table, int start_pos,
		int max_length);

void
print_vertical_htmltable (FILE* fout, HtmlTable* html_table, int start_pos,
			int max_length, double* fraction_stored);

FILE* errorfp;
char errorfilename[LARGE_BUFF_LENGTH];

void free_htmltable (HtmlTable* htmltable);

void color_char (FILE* fout, char c);

int 
main (int argc, char* argv[])
{
	FILE* ps_outfp; FILE* seqfp; 
	Sequence* seqs[MAXSEQ];
	int nseqs;
	Block* block;
	Matrix* no_pseudo_pssm;
	Matrix* pseudo_pssm_score;
	char tempname[LARGE_BUFF_LENGTH];
	char* strptr; int desc_length; char desc[SMALL_BUFF_LENGTH];
	int i, db_type, seq_type;	
	char outfilename[LARGE_BUFF_LENGTH];
	char currentoutfile[LARGE_BUFF_LENGTH];
	int aa_length;
	HtmlTable* html_table;
	int file_counter,  aa_printed;
	double threshold; 
	int gap_option, exp_option, diri_option, subtract_option;
	int table_size;
	double* fraction_stored;
	int seq_identity;

	table_size = 100;
	subtract_option = FALSE; diri_option =TRUE;
	init_frq_qij();
	getargs (argc, argv, outfilename, &seqfp, &threshold, 
			&gap_option, &exp_option, &seq_identity);


/* read in sequences */
        nseqs = 0;
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
        printf ("des is desc %s\n", desc);
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
        if (nseqs == MAXSEQ)
        {  fprintf (stderr, "WARNING: Maximum number of sequences = %d\n", nseqs);}
        /* reduce redundancy with query sequence 01/03/00
           so sequence weighting doesn't include more than 1 sequence
           representing query */
        remove_seqs_percent_identical_to_query (seqs, &nseqs,
					 (double) seq_identity);

 
	aa_length = get_length (seqs[0]);

        block = make_block (aa_length, 0, nseqs, seqs, FALSE);
      	fraction_stored = calculate_basic_aa_fraction (block);
	no_pseudo_pssm = PN_block_to_matrix (block, 2);
        pseudo_pssm_score = SIFT_prediction (block, diri_option, 
					gap_option, exp_option, subtract_option);

	file_counter = 1; aa_printed = 0;
	while ( aa_printed < block->width) {
	        sprintf (currentoutfile, "%s%d", outfilename, file_counter);
	        if ( (ps_outfp = fopen (currentoutfile, "w")) == NULL)
       		{
       		         fprintf (errorfp, "Cannot write to %s \n", 
							currentoutfile);
               		 exit (-1);
        	}
		html_table = initialize_table (table_size);
		html_table->ref_seq = seqs[0];
		make_coordinates ( aa_printed,
                         no_pseudo_pssm,  pseudo_pssm_score, html_table);
		print_vertical_htmltable (ps_outfp, html_table, aa_printed, aa_length, fraction_stored);
		fclose (ps_outfp);
		aa_printed += table_size; file_counter++;
	}

	free_block (block);
	free_matrix (pseudo_pssm_score);
	free_matrix (no_pseudo_pssm); 
	free_htmltable (html_table);
	fclose (errorfp);
	rm_file (errorfilename);

	exit (0);

} /* end of MAIN */

HtmlTable*
initialize_table (int size)
{
	int aa, pos;
	HtmlTable* htmltable;
	
	htmltable = (HtmlTable*) calloc (1, sizeof (HtmlTable));
	htmltable->length = size;
	htmltable->allowed_subst = (int*) calloc (size, sizeof (int));
	
	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		htmltable->table[aa] = (char *) calloc (size, sizeof(char));
		htmltable->intolerant_table[aa] = (char *) calloc (size, sizeof(char));
	}

	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		for (pos = 0; pos < size; pos++) {
			htmltable->table[aa][pos] = ' ' ;
			htmltable->intolerant_table[aa][pos] = ' ';
		}
	}
	return htmltable;
}

void
free_htmltable (HtmlTable* htmltable)
{
	int aa;

	for (aa = 0; aa < MATRIX_AA_WIDTH; aa++) {
		free (htmltable->table[aa]);
		free (htmltable->intolerant_table[aa]);
	}
	free (htmltable->allowed_subst);
	free(htmltable);
}

void
copy_aa_scores (struct aa_score list[20], 
		Matrix* no_pseudo_pssm, Matrix* pseudo_pssm, int pos) 
{
	int i;
	int aa;

	copy_values (list, pseudo_pssm, pos);
	for ( i = 0; i < 20; i++) {
		aa = list[i].aa;
		if (no_pseudo_pssm->weights[aa][pos] > 0.0) {
			list[i].bool_variable =  TRUE;
		}
	}
}

void
print_html_table (FILE* fout, HtmlTable* html_table, int start_pos, int max_length)
{
	int i, row, pos;

	fprintf (fout, "<html>");
	fprintf (fout, "<body bgcolor=white>");
	fprintf (fout, "<H1><center> Substitutions allowed at positions %d through %d</center><H1><BR>\n", start_pos +1, start_pos + html_table->length);
	fprintf (fout, "<BR><table cellspacing=0 border=0 width=0 cols=21>");
	for (row = 19; row >= 0 ; row--) {
		fprintf (fout, "<tr>");
		for (pos = 0; pos < html_table->length; pos++) {
			fprintf (fout, "<td>%c</td>", html_table->table[row][pos]);
		}
		fprintf (fout, "</tr>");
	}
	fprintf (fout, "<tr>");
	for (pos = start_pos; pos < start_pos + html_table->length && pos < max_length; pos++) {
		fprintf (fout, "<th>%d</th>", pos + 1);
	}
	fprintf (fout, "</tr>");
	fprintf (fout, "</table>");
	fprintf (fout, "</html>"); 
}


void
print_vertical_htmltable (FILE* fout, HtmlTable* html_table, int start_pos, int max_length, double* fraction_stored)
{
        int row, pos;
	fprintf (fout, "<HTML>"); 
        fprintf (fout, "<body bgcolor=white>");
        fprintf (fout, "<H1><center> Predictions for positions %d through %d</center></H1><BR>\n", start_pos +1, start_pos + html_table->length);
        fprintf (fout, "Threshold for intolerance is 0.05.<BR>");
	fprintf (fout, "Amino acid color code:  ");
	fprintf (fout, "nonpolar, <font color=green>uncharged polar</font>, ");
	fprintf (fout, "<font color=red>basic</font>, <font color=blue>acidic</font>. <BR>");	
	fprintf (fout, "Capital letters indicate amino acids appearing in the alignment, lower case letters result from prediction. <BR>");
	 fprintf (fout, "<b>'Seq Rep'</b> is the fraction of sequences that\n ");
	
	fprintf (fout, "contain one of the basic amino acids.  A low fraction");
	fprintf(fout, " indicates the position is either severely gapped ");
	fprintf (fout, "or unalignable and has little information.  Expect poor prediction at these positions. <BR>");

	fprintf (fout, "<BR><table cellspacing=0 border=0 width=0 cols=41>");
        fprintf (fout, "<th colspan=20 align=center>Predict Not Tolerated</th>");
	fprintf (fout, "<BR><th>Position</th><th>Seq Rep</th><BR>");
	fprintf (fout, "<th colspan=20 align=center>Predict Tolerated</th>\n");
	for (pos = 0; pos < html_table->length && pos + start_pos < max_length; pos++) {
                        fprintf (fout, "<tr>");
		/* print out intolerant subst */
		for (row = 19;row >=0; row--) {
			color_char (fout, html_table->intolerant_table[row][pos]);
		}
		/* print out pos */
		fprintf (fout, "<th>%d%c</th>", start_pos + pos +1, aa_btoa[html_table->ref_seq->sequence[start_pos + pos]]);
		fprintf (fout, "<td>%.2f</td>", fraction_stored[start_pos +pos]);
		/* print out the allowed subst */
		for (row = 0; row < 20 ; row++) {
			color_char (fout, html_table->table[row][pos]);
                }
                fprintf (fout, "</tr>\n");
        }
        fprintf (fout, "</table>\n");
	fprintf (fout, "</HTML>\n");
}

 
void
color_char (FILE* fout, char c)
{
	switch (c) {
/* amino acids with nonpolar side chain */
		case 'A': case 'a':
		case 'V': case 'v':
		case 'F': case 'f':
		case 'P': case 'p':
		case 'M': case 'm':
		case 'I': case 'i':
		case 'L': case 'l':
		case 'W': case 'w':
		case 'G': case 'g':
			fprintf (fout, "<td><font color=black>%c</font></td>", c);
			break;
		case 'D': case 'd':
		case 'E': case 'e':
			fprintf (fout, "<td><font color=blue>%c</font></td>", c);
			break;
		case 'R': case 'r':
		case 'H': case 'h':
		case 'K': case 'k':
			fprintf (fout, "<td><font color=red>%c</font></td>", c);
			break;
		case 'S': case 's':
		case 'T': case 't':
		case 'Y': case 'y':
		case 'C': case 'c':
		case 'N': case 'n':
		case 'Q': case 'q':
			fprintf (fout, "<td><font color=green>%c</font></td>", c);
			break;
		default:
			fprintf (fout, "<td><font color=gray>%c</font></td>", c);
			break;
	} /* end of switch*/
} 

void
make_coordinates ( int pos_shift,
			Matrix* no_pseudo_pssm, Matrix* pseudo_pssm, 
			HtmlTable* html_table)

{
	int pos;
	struct aa_score list[20];
	int list_size, j, k;

	for (pos = pos_shift ; pos < pos_shift + html_table->length && pos < pseudo_pssm->width;
	 pos++) {
		copy_aa_scores (list, no_pseudo_pssm, pseudo_pssm, pos);
		sort_list (list); /* list is sorted from highest to lowest */

		/* for file to have highest scoring printed out on top,
		have it printed out last */

	        list_size = 19;
	        while (list[list_size].score < INTOLERANCE_PROB_THRESHOLD) {
                        list_size--;

       		 }
		/* lower scoring amino acids at the bottom of the table */	
		/* assign allowed pos */
		 html_table->allowed_subst[pos - pos_shift] = list_size;
		 for (j = 0,k = list_size ; k >= 0; k--, j++) {
	/* observed, print as a capital.  if from pseudo, print it as a lower
	letter */
			if (list[k].bool_variable) {
				html_table->table[j][pos - pos_shift] = aa_btoa[list[k].aa];
			} else {
				html_table->table[j][pos - pos_shift] = tolower(aa_btoa[list[k].aa]);
			}
		}	

		/* intolerant pos */
		for (j = 0, k = list_size + 1;k < 20; k++, j++) {
			if (list[k].bool_variable) {
                                html_table->intolerant_table[j][pos - pos_shift] = aa_btoa[
list[k].aa];
                        } else {
                                html_table->intolerant_table[j][pos - pos_shift] = tolower(
aa_btoa[list[k].aa]);
                        }
                }
	
	} /* end of for pos */
	
}

void
getargs (int argc, char* argv[], char outfilename[LARGE_BUFF_LENGTH], 
	FILE** seqfp, double* threshold, int* gap_option, int* exp_option, 
	int* seq_identity)
{
	char seqfile[LARGE_BUFF_LENGTH];

	if (argc < 5) {
		printf ("allowed_subst_ps.c\n");
		printf ("prints out allowed substitutions\n"); 
	}
 	if (argc > 1) strcpy (outfilename, argv[1]);
	else {
		printf ("Enter  outfile\n");
		fgets (outfilename, LARGE_BUFF_LENGTH, stdin);
	}
	
        strcpy (errorfilename, outfilename);
        strcat (errorfilename, ".error");
        if ((errorfp = fopen (errorfilename, "w")) == NULL) {
                fprintf (stderr, "couldn't open file %s\n", errorfilename);
                exit (-1);
        }

	if (argc > 2) strcpy (seqfile, argv[2]);
	else {
		printf ("Enter aligned sequences file \n");
		fgets (seqfile, LARGE_BUFF_LENGTH, stdin);
	}
	if ((*seqfp = fopen (seqfile, "r")) == NULL) {
		fprintf (stderr, "Cannot open file %s \n", seqfile);
		exit (-1);
	}
        if (argc > 3) *threshold=atof(argv[3]);
        else
        {
                printf ("Enter threshold for intolerance\n");
                scanf ("%lf", threshold);
        }

        if (argc > 4) *gap_option = atoi (argv[4]);
        else {
                printf ("Enter gap option option 1 for allow everything\n");
                printf ("                        0 for ignore gaps\n");
                scanf ("%d", gap_option);
        }
        if (argc > 5) *exp_option = atoi (argv[5]);
        else {
                printf ("Enter exp option 0 - # diffaase\n");
                printf ("                1 - similarity scale\n");
                scanf ("%d", exp_option);
        }

	if (argc > 6) *seq_identity = atoi (argv[6]);
	else {
		*seq_identity = 100;
	}

} /*end of getargs */

Sequence*
read_sequence_from_filename (char filename[LARGE_BUFF_LENGTH])
{
        FILE* fp;
        Sequence* seq;

        if ((fp = fopen (filename, "r")) == NULL) {
                printf ("couldn't open %s\n", filename);
                exit (-1);
        }
        seq = read_a_sequence (fp, FASTA, AA_SEQ);
        fclose (fp);
        return seq;

} /* end of read_sequence_from_filename */

