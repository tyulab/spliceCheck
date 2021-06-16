#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LSTR 32768

//--------------------------------------------------------------------
int main(int argc,  char *argv[])
{

	FILE *infile;
	int  print=0;
	char *s, *idpos;


	if(argc!=3) {
		fprintf(stderr, "Usage: %s database seqid\n", argv[0]);
		exit(1);
	} // end if


   if( !(s=(char*)malloc(LSTR)) ) {
      printf("Error: alloc fail\n");
      exit(1);
   } // end if

   if( !(infile=fopen(argv[1], "r"))) {
      printf("Error: can't open %s\n", argv[1]);
      exit(1);
   } // end if

	while(fgets(s,LSTR,infile)) {

		if(s[0]=='>') {
			if(print) {
				break;
			} // end if
			else if(
				(idpos=strstr(s, argv[2]))!=NULL && // match [\|\>]PATTERN[\s\|]
				(*(idpos-1)=='|' || *(idpos-1)=='>') &&
				(*(idpos+strlen(argv[2]))=='|' || isspace(*(idpos+strlen(argv[2]))) )
			) { // not just strstr(): imagine YZ_HUMAN and XYZ_HUMAN
				print=1;
			} // end else if
		} // end if '>'

		if(print) {
			printf("%s", s);
		} // end if print

	} // end while

	fclose(infile);

	free(s);

// return nothing if seq. not found
/*
	if(!print) {
		fprintf(stderr, "Error: %s not found in %s\n", argv[2], argv[1]);
	} // end if
*/
	exit(0);
} // end main()
//--------------------------------------------------------------------
