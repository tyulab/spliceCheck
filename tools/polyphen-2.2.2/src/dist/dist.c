#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>

//--------------------------------------------------------------------

#define ONE_CHAIN     '_'
#define HET_CHAIN     '*'

#define MAX_ATOMS    150     // max atoms in res
#define MAX_RESIDUES   20000 // max residues in all chains

#define LSTR   120
#define OK     1 // return value
#define FAIL  -1

#define DEBUG 0
#define SKIP_H 1
#define FAST_MODE 1 // if set to 1, first look at Ca-Ca dist
#define UNDEF_SIZE -1 // size for het groups and unknown aa
#define SAFETY_MARGIN 3 // if residues are smaller than mean size

//--------------------------------------------------------------------


typedef struct {
   double x;
   double y;
   double z;
} point;

typedef struct {
   int    num;
   char   type[3];
   char   aa[4];
   char   chain;
	int    is_mainchain;
   char   resnum[6];  // !! may be 2C, 2B, 2A, 2, 3,...
   point  coord;
   double occup;
   double bfac;
} atom;

typedef struct {
   char   num[6];
   char   chain;
   char   name[4];
   int    natoms;
   atom   *atoms;
	int    het;
	int    ter;
	double size;
} Res;

typedef struct {
	char   name[4];
	double size;
} AaSize;

typedef struct {
	int    a1; // atom number from 1st res
	int    a2; // atom number from 2nd res
	double d;  // distance betw these two atoms
} Contact;

AaSize aa_size[] = {
 {"ALA", 3.385}, {"ARG", 7.864}, {"ASN", 5.101}, {"ASP", 5.137}, {"CYS", 4.310},
 {"GLN", 6.055}, {"GLU", 6.148}, {"GLY", 3.239}, {"HIS", 6.118}, {"ILE", 5.284},
 {"LEU", 5.186}, {"LYS", 7.102}, {"MET", 5.784}, {"PHE", 6.756}, {"PRO", 4.147},
 {"SER", 3.887}, {"THR", 4.405}, {"TRP", 7.716}, {"TYR", 8.117}, {"VAL", 4.477}
};

//-- globals: --------------------------------------------------------


int tot_het_res=0;

int    internal_dist; // calculate contacts within chain
//char   file[120];
char *file;
char   pdbid[20];
char   resnum[20];
char   chain;
int    all_res; // calc dist for all residues, set to 0 if resnum given
double thresh;
int    mainchain_cont;
int    all_cnt_below_thresh;

char format[80];

//--------------------------------------------------------------------

int   parseAtomString(char *s, atom *a);
void  Print_Atom(atom x);
// FILE  *openPDB(char *p);
void  flushRes(Res *p);
Res   *newRes(void);
int   readAtoms(char *file, Res **rr);
int   moveHetRes(Res **ra, int n, Res **rh);
void  assignMainchain(Res *r);
void  file2pdbid(char *file, char *pdbid);  

double DistSquared(point p1, point p2);
int setParams(int argc,  char *argv[]);
int countContacts(Res *r1, Res *r2, Contact *c); 

//--------------------------------------------------------------------

int main(int argc,  char *argv[])
{

  Res **aares, **hetres = NULL, *r1, *r2;
	Contact *pcnt;
	int cnum;
	if ((file = (char *) calloc( (strlen(argv[argc - 1]) + 1), sizeof(char))) == NULL) {
	  printf("No memory left\n");
	  exit(1);
	}
	strcpy(file, argv[argc-1]);
   int naares, nhetres;
	int i,j, k;
	int k1, k2;
	int found=0;

	if( DEBUG ) setbuf(stdout, NULL);

	if(!setParams(argc, argv))
		return 1;

	strcpy( format,

		"%-6s %3s "
		"%5s%c %3s  %-2s %-5d %c"
		" >...< "
		"%c %5s%c %3s  %-2s %-5d"
		"   %6.3f\n"

	);


/*
printf("\nid=%s allres=%d resnum=%s chain=%c intdist=%d thresh=%f",
	pdbid, all_res, resnum, chain, internal_dist, thresh);
return 0;
*/


   if(!(aares=(Res**)calloc(MAX_RESIDUES, sizeof(Res*)))) {
      printf("Error: can't allocate memory for aa residues\n");
      return 1;
   } // end if

   if(!(pcnt=(Contact*)calloc(MAX_ATOMS*MAX_ATOMS, sizeof(Contact)))) {
      printf("Error: can't allocate memory for contacts\n");
      return 1;
   } // end if

   if(!(naares=readAtoms(file, aares))) {   // sets global tot_het_res
      printf("Error: can't find residues\n");
      return 1;
	} // end if

	if(DEBUG) printf("\n%s\naares=%d hetres=%d", pdbid, naares, tot_het_res);

	nhetres=0;

	if( tot_het_res ) {

	   if(!(hetres=(Res**)calloc(tot_het_res, sizeof(Res*)))) {
	      printf("Error: can't allocate memory for HET residues\n");
	      return 1;
	   } // end if

		nhetres = moveHetRes(aares, naares, hetres);

	} // end if tot_het_res

	if(DEBUG) printf("\nhetres moved=%d\n", nhetres);

	for(i=0; i<naares; i++) {
		if( aares[i]==NULL )
			continue;
		assignMainchain(aares[i]);

//for(j=0; j<aares[i]->natoms; j++) Print_Atom(aares[i]->atoms[j]); 
	} // end for

	////////////////////////////////////////////////////
	//////// main cycle
	///////////////////////////////////////////////////

	for(i=0; i<naares; i++) {

		if( aares[i]==NULL )
			continue;

		if( !all_res && (strcmp(resnum, aares[i]->num) || chain !=aares[i]->chain ))
			continue;

		found=1;

		for(j=0; j<nhetres; j++) {

				r1 =  aares[i];
				r2 = hetres[j];

				if((cnum=countContacts(r1, r2, pcnt))==0)
					continue;

				for(k=0; k<cnum; k++) { // cycle over contacts

					k1 = pcnt[k].a1;
					k2 = pcnt[k].a2;

					printf(
						format,
	
						pdbid, "HET",
						r1->num, r1->chain, r1->name,
						r1->atoms[k1].type,
						r1->atoms[k1].num,
						r1->atoms[k1].is_mainchain?'m':'s',
	
						r2->atoms[k2].is_mainchain?'m':'s',
						r2->num, r2->chain, r2->name,
						r2->atoms[k2].type,
						r2->atoms[k2].num,
	
						pcnt[k].d	

				   ); // end printf 

				} // end for 

		} // end for hetres

		for(j=all_res?i+1:0; j<naares; j++) {

				if( aares[j]==NULL || i==j )
					continue;

				r1 =  aares[i];
				r2 =  aares[j];

				if( !internal_dist && r1->chain==r2->chain) continue;
				if(  internal_dist && r1->chain==r2->chain
					 && j==i+1 ) continue;

				if((cnum=countContacts(r1, r2, pcnt))==0)
					continue;

				for(k=0; k<cnum; k++) { // cycle over contacts

					k1 = pcnt[k].a1;
					k2 = pcnt[k].a2;

					printf(
	
						format,
	
						pdbid, (r1->chain==r2->chain)?"INT":"EXT",
						r1->num, r1->chain, r1->name,
						r1->atoms[k1].type,
						r1->atoms[k1].num,
						r1->atoms[k1].is_mainchain?'m':'s',

						r2->atoms[k2].is_mainchain?'m':'s',
						r2->num, r2->chain, r2->name,
						r2->atoms[k2].type,
						r2->atoms[k2].num,

						pcnt[k].d

				   ); // end if
				} // end for

		} // end for aares


	} // end for aares

	if( !found ) {
      printf("Error: can't find residue %s%c in %s\n", resnum, chain, file);
	} // end if
	free(file);
	return 0;

} // end main
//--------------------------------------------------------------------
int readAtoms(char *file, Res **rr)
{

   char s[LSTR];
	FILE *fp;
   atom a;
   //Res  **rr, *r;
	Res  *r;
   int  tot_res, i, j; // residue and atom counters
	int  k;
	double size;

//	if(!(fp=openPDB(id))) return NULL;
	if(!(fp=fopen(file, "r"))) {
		printf("Error: can't open file %s\n", file);     
		return 0;
	} // end if

   i = j = 0;
   r = newRes();

   while(fgets(s, LSTR, fp)) {

      if(
				strstr(s,"ATOM")==s ||
				strstr(s,"HETATM")==s && !strstr(s, "HOH") && !strstr(s, "DIS")
		) {

         if(parseAtomString(s, &a)==FAIL) {
            printf("Error: can't parse %s\n", s);
         } // end if
         else { // atom parse ok

				if( SKIP_H && a.type[0]=='H' )
					continue;

            if( j>0 && (strcmp(a.resnum, r->num) || a.chain!=r->chain )) { // new residue
//             if( /* check for completeness */ ) {}

               r->natoms=j;
               rr[i++] = r; // push residue
					if(r->het) tot_het_res++;

					if(i==MAX_RESIDUES-1) { // -1 to store current res.
		            printf("Warning: max number of residues (%d) exceeded\n", MAX_RESIDUES);
						break;
					} // end if

               r = newRes();
               j=0;

            } // end if new residue

            r->atoms[j++]=a;

				if(j==MAX_ATOMS) {
		         printf("Warning: max number of atoms (%d) in residue %s %s exceeded\n", MAX_ATOMS, a.resnum, a.aa);
					break;
				} // end if

            strcpy(r->num, a.resnum);
            strcpy(r->name,  a.aa);
            r->chain = a.chain;
				r->het = (strstr(s,"HETATM")==s)?1:0;

         } // end else atom parsed ok

      } // end if ATOM
		else if( strstr(s, "TER")==s ) {
			r->ter=1;
		} // end else if
		else if( strstr(s, "ENDMDL")==s) {
			break;
		} // end else if

   } // end while

	fclose(fp);

	if( j>0 ) { // to prevent pushing if there are no residues actually
	   r->natoms=j;
	   rr[i++] = r;
		if(r->het) tot_het_res++;
	} // end if j>0

	tot_res =i;

	// assign res. sizes
	for(i=0; i<tot_res; i++) {

		size = UNDEF_SIZE; // dummy value for het and unknown residues

		for(k=0; k<20; k++) {
			if( strcmp(rr[i]->name, aa_size[k].name)==0) {
				size = aa_size[k].size;
				break;
			} // end if
		}  // end for

		rr[i]->size = size;

	} // end for

   return tot_res;

} // end readAtoms()
//--------------------------------------------------------------------
int parseAtomString( char *s, atom *a )
{
/*
from http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
COLUMNS        DATA TYPE       FIELD         DEFINITION
---------------------------------------------------------------------------------
 1 -  6        Record name     "ATOM  "
 7 - 11        Integer         serial        Atom serial number.
13 - 16        Atom            name          Atom name.
17             Character       altLoc        Alternate location indicator.
18 - 20        Residue name    resName       Residue name.
22             Character       chainID       Chain identifier.
23 - 26        Integer         resSeq        Residue sequence number.
27             AChar           iCode         Code for insertion of residues.
31 - 38        Real(8.3)       x             Orthogonal coordinates for X
39 - 46        Real(8.3)       y             Orthogonal coordinates for Y
47 - 54        Real(8.3)       z             Orthogonal coordinates for Z
55 - 60        Real(6.2)       occupancy     Occupancy.
61 - 66        Real(6.2)       tempFactor    Temperature factor.
73 - 76        LString(4)      segID         Segment identifier, left-justified.
77 - 78        LString(2)      element       Element symbol, right-justified.
79 - 80        LString(2)      charge        Charge on the atom.

ATOM     90  NE2BGLN A  14      21.817   8.764   6.180  0.38 21.41           N
ATOM    286  CB  SER    40A     35.314  32.534  20.539  1.00100.00      103L 384
ATOM     35 1HG  GLU     2      27.255  23.117  15.549  1.00 21.11      2SN3 134
ATOM     91  H   GLN A  14      19.939   6.870   2.319  1.00 18.10           H
ATOM    332  CA  LYS A  27       3.311  -6.909   2.519  1.00  9.29           C
ATOM    511  HE2 HIS A  35      -7.660  -6.633   3.137  1.00 14.00           H
ATOM    517  CB BLEU A  36      -8.323  -2.519   6.531  0.47 15.21           C
ATOM  11459  CB  SER H  36A    -19.325  15.256  17.735  1.00 54.57      1PPB 680

01234567890123456789012345678901234567890123456789012345678901234567890123456789
0         1         2         3         4         5         6         7

*/

   if(

       sscanf(s+6,  "%d%*2c%2s", &a->num, a->type) +
       sscanf(s+17, "%3s", a->aa ) +
       sscanf(s+22, "%5s", (char *) &a->resnum ) +
       sscanf(s+30, "%lf%lf%lf%4lf%lf",
          &(a->coord.x), &(a->coord.y), &(a->coord.z), &a->occup, &a->bfac
       ) != 9
     ) return FAIL;

   if((a->chain=s[21])==' ')
      a->chain=ONE_CHAIN;

   return OK;

} // end parseAtomString()
//--------------------------------------------------------------------
void Print_Atom(atom x)
{
// printf("%5d-%2s-%3s-%2s-%5s-%.2lf-%6.2lf\n",
   printf("%d-%s-%s-%c-%s-%d\n",
      x.num, x.type, x.aa, x.chain, x.resnum, x.is_mainchain);

} // end Print_Atom()
//--------------------------------------------------------------------
void assignMainchain(Res *r)
{

	int i,j;
	char *types[] = { "N", "CA", "C", "O" };

	// late initializing but necessary
	for(i=0; i<r->natoms; i++) {
	  r->atoms[i].is_mainchain=0;
	} // end for

	for(j=0; j<4; j++) {
		for(i=0; i<4; i++) {
			if(strcmp(r->atoms[i].type, types[j])==0) {
				r->atoms[i].is_mainchain = 1;
				break;
			} // end if
		} // end for i
	} // end for j

	

} // end assignMainchain()
//--------------------------------------------------------------------
/*
FILE *openPDB(char *id)
{
	FILE *fp;

//#define PDB "E:\\data\\pdb"
//#define PDB "C:\\data\\pdb"
//#define PDB "/data/pdb"


   char s[120];

	if( PQS ) {
   	sprintf(s, "/bork/coot3/ramensky/PQS/%s.pqs", id);
	}// end if
	else {
   	sprintf(s, "/data/pdb/%s.brk", id);
	} // end else

//   sprintf(s, "E:\\data\\pdb\\pdb%s.ent", id);
//   sprintf(s, "%s\\%s.pdb", PDB, id);
//   sprintf(s, "%s/%s.brk", PDB, id);
//   sprintf(s, "%s.pdb",  p);

	if(!(fp=fopen(s, "r")))
      printf("Error: can't open file %s\n", s);

   return fp;

} // end Open_PDB()
*/
//--------------------------------------------------------------------
void flushRes(Res *p)
{

   free(p->atoms);
   free(p);

} // end allocRes()
//--------------------------------------------------------------------
Res *newRes(void)
{

   Res *p;

   if((p=(Res*)calloc(1, sizeof(Res)))==NULL) {
      printf("Error: can't allocate memory for residue\n");
      return NULL;
   } // end if

   if((p->atoms=(atom*)calloc(MAX_ATOMS, sizeof(atom)))==NULL ) {
      printf("Error: can't allocate memory for atoms\n");
      return NULL;
   } // end if

   return p;

} // end newRes()
//--------------------------------------------------------------------
int countContacts(Res *r1, Res *r2, Contact *c)
{

	int n1, n2; // number of atoms
	int i1, i2, k=0;
	double dd, mindd=1000;

	if(
			FAST_MODE &&
			r1->size!=UNDEF_SIZE && r2->size!=UNDEF_SIZE
	) {

		if( r1->natoms>1 && r2->natoms>1 )
			dd = DistSquared(r1->atoms[1].coord, r2->atoms[1].coord);
		else 
			dd = DistSquared(r1->atoms[0].coord, r2->atoms[0].coord);

		if(sqrt(dd) > r1->size+r2->size+thresh+SAFETY_MARGIN) {
			return 0;
		} // end if

	} // end if

	n1 = r1->natoms;
	n2 = r2->natoms;

	for(i1=0; i1<n1; i1++) {
		for(i2=0; i2<n2; i2++) {

			
		  if (
		      ((dd=DistSquared(r1->atoms[i1].coord, r2->atoms[i2].coord)) >= thresh*thresh)
		      ||
		      (!mainchain_cont && r1->atoms[i1].is_mainchain && r2->atoms[i2].is_mainchain)
		      )
		    continue;

			if(all_cnt_below_thresh) {
				(c+k)->a1 = i1;
				(c+k)->a2 = i2;
				(c+k)->d  = sqrt(dd);
				k++;
			} // end if
			else if(dd<mindd) {
				c->a1 = i1;
				c->a2 = i2;
				c->d  = sqrt(dd);
				mindd = dd;
				k=1;
			} // end else 

		} // end for i2
	} // end for i1


	return k;


} // end sub
//--------------------------------------------------------------------
/*
double ResDistInf(Res *r1, Res *r2, int *x1, int *x2)
{

	double d, mindist=10000;
	int n1x, n2x;


	for(i1=0; i1<n1; i1++) {
		for(i2=0; i2<n2; i2++) {

			if(!mainchain_cont && r1->atoms[i1].is_mainchain && r2->atoms[i2].is_mainchain)
				continue; 
			if((d=DistSquared(r1->atoms[i1].coord, r2->atoms[i2].coord))<mindist){
				mindist=d;
				n1x=i1; n2x=i2;
			}
		} // end for
	} // end for

	*x1=n1x; *x2=n2x;
	return 
		sqrt(mindist);

} // end RedDist()
*/
//--------------------------------------------------------------------
double DistSquared(point p1, point p2)
{

		return
			(p1.x-p2.x)*(p1.x-p2.x)+
			(p1.y-p2.y)*(p1.y-p2.y)+
			(p1.z-p2.z)*(p1.z-p2.z);

} // end DistSquared()
//--------------------------------------------------------------------
int moveHetRes(Res **ra, int n, Res **rh )
{

	int i, k=0; // k: het res counter
	char chain = '*'; // dummy starting value
	int move_it = 0;

	for(i=n-1; i>=0; i-- ) {

		if( chain!=ra[i]->chain ) move_it=1;

		// !switching on TER may be an open question...
		if( ra[i]->ter==1 || ra[i]->het==0 ) move_it=0;

		chain = ra[i]->chain;

		if( ra[i]->het && move_it ) {
			rh[k++] = ra[i];
			if(DEBUG) printf("\n residue %s%c moved", ra[i]->num, ra[i]->chain);
			ra[i] = NULL;
		} // end if

	} // end for

	return k;

} // end moveHetRes()
//--------------------------------------------------------------------------
int setParams(int argc,  char *argv[])
{

	int i;
	int input_ok=1;


	all_cnt_below_thresh=0;
   internal_dist=0;
	mainchain_cont=0;
   all_res=1;
	resnum[0] = '\0';
	chain = '\0';

	if(argc<3) {
		input_ok=0;
	} // end if
	else {

	  //	  strcpy(file, argv[argc-1]);
	  file2pdbid(file, pdbid);

		thresh = atof(argv[argc-2]);

		i=1;

		while(i<argc-2) {

			if( argv[i][0]!='-' ) {
				input_ok=0;
				break;
			} // end if
			else if( argv[i][1]=='i' ) {
				internal_dist=1;
			}
			else if( argv[i][1]=='m' ) {
				mainchain_cont=1;
			}
			else if( argv[i][1]=='a' ) {
				all_cnt_below_thresh=1;
			}
			else if( argv[i][1]=='r' ) {
				i++; // -r num chain, so jump 
				strcpy(resnum, argv[i]);
				i++;
				chain = argv[i][0];

				all_res=0;
			}
			else { // unknown option
				input_ok=0;
            break;     
			} // end else

			i++;

		} // end while

	} // end else 


	if( !input_ok ) {
		printf(
			"\nUsage: %s  [options] thresh file"
			"\n where options are"
			"\n  -r resnum chain, %c if one chain"
			"\n  -i calculate intrachain contacts"
			"\n  -a show all contacts below threshold"
			"\n  -m consider mainchain-mainchain distances\n\n",

			argv[0], ONE_CHAIN
		);
		return 0;
	} // end if

	return OK;

} // end setParams()
//--------------------------------------------------------------------------
// ...dir/dir/file.ext -> file
void file2pdbid(char *file, char *pdbid)
{

	char *dot, *slash;

	char *tmp;
	if ((tmp = (char *) calloc(strlen(file) + 1, sizeof(char))) == NULL) {
	  printf("No memory left");
	  exit(1);
	}


	strcpy(tmp, file);


	if((dot=strrchr(tmp, '.'))) {
		*dot = '\0';
	} // end if

	if((slash=strrchr(tmp, '/'))) {   
		strcpy(pdbid, slash+1);
	} // end if
	else {
		strcpy(pdbid, tmp);
	} // end else 

} // end file2pdbid()
//--------------------------------------------------------------------------
