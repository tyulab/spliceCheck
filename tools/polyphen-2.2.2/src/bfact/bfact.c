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


//-- globals: --------------------------------------------------------


int tot_het_res=0;

char   pdbid[200];
char   resnum[20];
char   chain;

//--------------------------------------------------------------------

int   parseAtomString(char *s, atom *a);
void  flushRes(Res *p);
Res   *newRes(void);
int   readAtoms(char *pdbid, Res **rr);
int   moveHetRes(Res **ra, int n, Res **rh);

void Mean_And_Disp(double *a, int n, double *p_mean, double *p_rmsd);
//--------------------------------------------------------------------

int main(int argc,  char *argv[])
{

  Res **aares, **hetres;

  int naares, nhetres;
  int i,j;
  int found=0;

  double mean, rmsd;
  double mean2, rmsd2;
  double *Bfac, *Bfac2;
  int tot_at, tot_at2;

  if(argc!=4) {
    printf("Usage: %s file resnum chain\n", argv[0]);
    return 0;
  } // end if
  else {
    strcpy(pdbid,  argv[1]);
    strcpy(resnum, argv[2]);
    chain = argv[3][0];
  } // end else


  //----- allocate -----------------------------------------------
  if(!(aares=(Res**)calloc(MAX_RESIDUES, sizeof(Res*)))) {
    printf("Error: can't allocate memory for aa residues\n");
    return 0;
  } // end if

  if(!(naares=readAtoms(pdbid, aares))) {   // sets global tot_het_res
    printf("Error: can't find residues for chain %c in %s\n", chain, pdbid);
    return 0;
  } // end if

  if(!(Bfac=(double*)calloc(naares*MAX_ATOMS, sizeof(double)))) {
    printf("Error: can't allocate memory for B-factors\n");
    return 0;
  } // end if

  if(!(Bfac2=(double*)calloc(MAX_ATOMS, sizeof(double)))) {
    printf("Error: can't allocate memory for B-factors\n");
    return 0;
  } // end if

  //----- move hetres --------------------------------------------
  nhetres=0;

  if( tot_het_res ) {

    if(!(hetres=(Res**)calloc(tot_het_res, sizeof(Res*)))) {
      printf("Error: can't allocate memory for HET residues\n");
      return 0;
    } // end if

    nhetres = moveHetRes(aares, naares, hetres);

  } // end if tot_het_res


  //----- calc mean/rmsd -----------------------------------------
  tot_at=0;

  for(i=0; i<naares; i++) {
    if( aares[i]==NULL || chain!=aares[i]->chain )
      continue;

    //  assignMainchain(aares[i]);

    for(j=0; j<aares[i]->natoms; j++) {
      Bfac[tot_at++] = aares[i]->atoms[j].bfac;
    } // end for
  } // end for

  if(!tot_at) {
    printf("Error: zero length of chain %c\n", chain);
    return 0;
  } // end if

  Mean_And_Disp(Bfac, tot_at, &mean, &rmsd);

  printf("Chain: mean=%-6.2f rmsd=%-6.2f", mean, rmsd);


  ////////////////////////////////////////////////////
  //////// main cycle
  ///////////////////////////////////////////////////

    tot_at2 = 0;

    for(i=0; i<naares; i++) {

      if( aares[i]==NULL )
	continue;

      if( strcmp(resnum, aares[i]->num) )
	continue;

      found=1;

      for(j=0; j<aares[i]->natoms; j++) {
	Bfac2[tot_at2++] = aares[i]->atoms[j].bfac;
      } // end for

    } // end for aares

    if( !found ) {
      printf("Error: can't find residue %s%c in %s\n", resnum, chain, pdbid);
    } // end if

    Mean_And_Disp(Bfac2, tot_at2, &mean2, &rmsd2);


    printf("\tResidue: mean=%-6.2f rmsd=%-6.2f", mean2, rmsd2);

    printf("\tNormedB: %-6.2f\n", (rmsd>0)?(mean2-mean)/rmsd:0 );

    return 0;
} // end main
//--------------------------------------------------------------------
int readAtoms(char *id, Res **rr)
{

  char s[LSTR];
  FILE *fp;
  atom a;
  //Res  **rr, *r;
  Res  *r;
  int  tot_res, i, j; // residue and atom counters

  if(!(fp=fopen(pdbid, "r"))) {
    printf("Error: can't open %s\n", pdbid);
    return 0;
  }


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

	if( a.chain!=chain)
	  continue;

	if( j>0 && (strcmp(a.resnum, r->num) || a.chain!=r->chain )) { // new residue
	  //             if( /* check for completeness */ ) {}

	  r->natoms=j;
	  rr[i++] = r; // push residue
	  if(r->het) tot_het_res++;

	  if(i==MAX_RESIDUES-1) { // -1 to store current res.
	    printf("Warning: max number (%d) of residues exceeded\n", MAX_RESIDUES);
	    break;
	  } // end if

	  r = newRes();
	  j=0;

	} // end if new residue

	r->atoms[j++]=a;

	if(j==MAX_ATOMS) {
	  printf("Warning: max number (%d) of atoms in residue %s %s  exceeded\n", MAX_ATOMS, a.resnum, a.aa);
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
     sscanf(s+22, "%5s", &a->resnum ) +
     sscanf(s+30, "%lf%lf%lf%4lf%lf",
	    &(a->coord.x), &(a->coord.y), &(a->coord.z), &a->occup, &a->bfac
	    ) != 9
     ) return FAIL;

  if((a->chain=s[21])==' ')
    a->chain=ONE_CHAIN;

  return OK;

} // end parseAtomString()
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

  if((p=malloc(sizeof(Res)))==NULL) {
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
      if(DEBUG) printf("\n resdiue %s%c moved", ra[i]->num, ra[i]->chain);
      ra[i] = NULL;
    } // end if

  } // end for

  return k;

} // end moveHetRes()
//--------------------------------------------------------------------------
// rmsd = sqrt(dispersion)
void Mean_And_Disp(double *a, int n, double *p_mean, double *p_rmsd)
{

  int i;
  register double tmp, tmp_mean=0, tmp_mean_of_square=0;

  for(i=0; i<n; i++) {
    tmp_mean += a[i];
    tmp_mean_of_square += a[i]*a[i];
  }

  tmp_mean /= n;
  tmp_mean_of_square /= n;

  *p_mean = tmp_mean;
  tmp = tmp_mean_of_square - tmp_mean*tmp_mean;

  if( tmp < 0 )
    tmp = 0;

  *p_rmsd = sqrt(tmp);

} // end Mean_And_Disp();
//--------------------------------------------------------------------

