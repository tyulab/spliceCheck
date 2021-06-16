/* (C) Copyright 2000, Fred Hutchinson Cancer Research Center */
/* Use, modification or distribution of these programs is subject to */
/* the terms of the non-commercial licensing agreement in license.h. */   

#ifndef _PSIBLAST_H_
#define _PSIBLAST_H_

#define EXTERN
#include <assert.h>
#include <string.h>

extern FILE* errorfp;

struct alignment_pair {
	Sequence* query;
	Sequence* subject;
	int score;
	double evalue;
	int length;
	struct alignment_pair* next; /* just in case there's a second alignment
					for this particular subject */
};

typedef struct alignment_pair Aligned_Pair;

void write_sequences_to_msa (char alignfile[LARGE_BUFF_LENGTH], Block* block);

void add_sequence_to_msa (char query_msa_filename[LARGE_BUFF_LENGTH], Sequence* seq);


void psiblast_system_call_first (char query_seq_file[], char database[],
                         char result_file[], FILE* outfp);

void psiblast_system_call_not_first (char database[], char result_file[],
                char alignment_in_filename[],char query_filename[], FILE* outfp);

Aligned_Pair*
read_psiblast_entry ( FILE* fp, char Buffer[LARGE_BUFF_LENGTH]);

/*read 4_alignment_lines*/
/* pass in first line as Buffer, reads Query, subject alignment, and the 
   following newline.  reads the next line into Buffer */

void
read_4_alignment_lines (char Buffer[LARGE_BUFF_LENGTH], Aligned_Pair* alignment,
			FILE* fp, int get_start_pos);


/* void process_sequence_line(Sequence* seq, char* buff); */

Aligned_Pair* initialize_Aligned_Pair (char Buffer[LARGE_BUFF_LENGTH]);

void reading_alignment_at_score_line (Aligned_Pair* alignment, 
				char Buffer[LARGE_BUFF_LENGTH], FILE* fp);

/* fp is where the chekcpoint file is written to, outfp is for comments*/
void block_to_checkpoint_file (Block* block, FILE* fp, FILE* outfp);

void psiblast_system_call (char chkpoint_filename[LARGE_BUFF_LENGTH],
			   char database[LARGE_BUFF_LENGTH],
			   char result_file[LARGE_BUFF_LENGTH],
			   char query_seq_file[LARGE_BUFF_LENGTH],
			   FILE* outfp);

void formatdb_system_call (char database[LARGE_BUFF_LENGTH]);

/*int read_aa_sequence (Sequence* seq, int seq_type, int start_pos, char* str);
*/
#endif

