/*
 *   Prest-Plus Version 4.09
 *  
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>

#include <algorithm>
#include <map>
#include <string>

using namespace std;

#include "options.cpp"
#include "utils.cpp"
#include "em-newton.cpp"



class pedmemb {
 public:
  int data[7];

  void operator = (int a[7]) { for(int i=0;i<7;i++) data[i] = a[i]; }
  int& operator[] (int i) { return data[i]; }
};

typedef vector<pedmemb> pedigree;

vector<pedigree> ped;



typedef int geno_t;


#define MANCE 2000   /* max no. of ancestors for an individual */
#define MAXSTR 5000  /* max characters for a string */

#define NRMLRT  11         /* no. of relationships considered in set A for
			      MLRT excluding the null relationship */

#define NR 10               /* no. of relationships considered as the null */


int debug = 0; /* Debugging flag , also sets srand(0); */

int MC = 200;              /* max no. of chiasma   */

double MYEPS = 0.5e-6; /* EM termination threshold */

double PTHRESH = 0.2;      /* threshold for EIBD, AIBS and IBS */
int MITER = 5000;          /* max no. of iteration in the EM algorithm */

int NREP = 100000;         /* maximum number of replicates for MLRT */

double ALPHA = 0.1; /* Precision of MLRT in relation with p-value */
/* (if 0 then the old step-wise method is used) */
/* The MLRT results will be 95% of the time between
 * (1-2*ALPHA)*p and (1+2*ALPHA)*p
 */

double errorrate = 0.01;   /* approximate genotyping error rate */

int OPTION = 1; /* If 2 then perform MLRT */

int APED = 0; /* If 1 then check across pedigrees */
int WPED = 0; /* If 1 then check within pedigrees */

int MAP = 0; /* If 0 then no map (spacings) is provided */
int ALLELEFREQ = 0; /* If 0 then no allele frequencies in map file */

int all_pairs = 0; /* For pairs that are related but the relationship is not modelled,
 * setting this to 1 will treat them as unrelated, otherwise there will be no output */

/* The following variables set which individuals are included in internal allele frequency computation */
int include_founders = 1;
int include_nonfounders = 0;

int include_cases = 1;
int include_controls = 1 ;
int include_missing = 1;


/* If we want to write the internally computed frequency files */
int write_frequency_file = 0;


/* Options to speedup EM algorithm */
int ibd_useNewton = 1;
int ibd_check10 = 0;


int stat_ibs = 1, stat_ibd = 1;

/* For string alignment in output file */

int max_fid = 4, max_iid = 4; // maximum length of fid and iid
int max_mark = 4;

int mmark;

/* Prefix for output files */
char output_prefix[512] = "prest";

// File names for output files
char file_summary[1000], file_results[2000], file_warnings[1000], file_errors[1000];


/* p012: kinship coef.,   0, 1, 2 IBD sharing */
static double p012[20][3+1]=     {   {0      , 0   , 0   , 0   },
				     {.25    , .25 , .5  , .25 }, // 1
				     {.125   , .5  , .5  , 0   }, // 2
				     {.125   , .5  , .5  , 0   }, // 3
				     {.125   , .5  , .5  , 0   }, // 4
				     {.0625  , .75 , .25 , 0   }, // 5
				     {0      , 1   , 0   , 0   }, // 6
				     {.0625  , .75 , .25 , 0   }, // 7
				     {.03125 , .875, .125, 0   }, // 8
				     {.1875  , .375, .5  , .125}, // 9
				     {.25    , 0   , 1   , 0   }, // 10
				     {.5     , 0   , 0   , 1   }};// 11


/* pedigree file name */
char ped_name[256] = "(nofile)";

/* marker and genotype data files */
char chromnameidx[100][100] = { "(nofile)" }, chromnameped[100][100] = { "(nofile)" };

/* nochrom: no. of chromosomes input files */
int nochrom = 0; 


/* for pedigree info. */

/* noped: no. of pedigrees in the pedgree file
   nomem: pedigree id, size of the ped
   ped: id, father id, mother id,
   recoded fatherid, recoded motherid, gender
   norel: no of each relationship considered
   (1=full-sib, 2=half-sib, 3=grandparent-grandchild,
   4=avuncular, 5=f-cousin, 6=unrelated,
   7=half-avuncular, 8=half-first-cousin, 9=half-sib+first-cousin,
   10 = parent-offspring)
   rel: id1 id2 reltype */
int    pedid_now, noped, numpedid;
int    **nomem,  **norel, ***rel;
//int    ***ped;



/* for marker info */

/*
   nomark: no of markers
   noalle: no of allels
   q: allele frequencies
   qcum: cumulated allele frequencies
   theta: spacing between markers in recombination fraction
   cenm: spacing in centi Morgan
   maptype: 1 = no-interference model, 2 = Kosambi */
int **noalle, maptype;
double ***q, ***qcum, **theta, **cenm, chromlength[100];

int nomark[100];
int *ngeno, ***genotype;

char *** names; // names of markers

// The following will hold the conditional probabilities for every pair, used in the IBD calculations

double *a, *b, *c;

int countmlrt = 0;
geno_t    *simugeno1, *simugeno2;
int **nopairsvalid;

/* for check the allele freq. */
double checkmarker[NR+1][2+1], checkallefreq[3+1][NR+1][5+1];

int    *valid,  **norelwdata;
int *ibs1, *ibs2;

/* for EIBD, AIBS, IBS */

/* for variance calulation */

double ***eeibdibd, ***eeibdeibdibd;
double ***eaibsibd, ***eaibsaibsibd;
double **eibsibd, **eibsibsibd;

/* statistic, mu, std, z_obs, p-value(2-sided) */
double statresult[3+1][5+1];


/* String buffer used for file reading */
char *buf;
int buf_len;

/* NA string */
char NA[5] = "NA";


/* Output files */
FILE *fp_summary, *fp_results, *fp_warnings;

// Support for running prest on a cluster
int cluster_n = 0, cluster_id = 0, cluster_join = 0;


int alternativeRel = 0; // alternative relationship specified

// Assume uniform marker spacings
double uniformSPC = -1;


int setUniformFreq = 0;

int _sim_pedno = 1;
int _sim_rel = 0;
int _sim_n = 0;

double pMLRT;


#include "prest.h"

#include "readped.cpp"
#include "readmap.cpp"


double pIBS0 = 0;
double pIBS2 = 0;
double pIBS0_th = 0;
double pIBS2_th = 0;

double pMarkov(int n,int k, double p)
{
  vector<double> P(k+1), P2(k+1);
  int i,j,t;

  P[0] = 1;
  for(t = 1; t<=k; t++) P[t] = 0;

  for(i = 1; i<=n; i++) {
    double s = 0;
    for(t = 0; t<k; t++) s = s + (1-p)*P[t];
    P2[0] = s;

    for(t = 1; t<k; t++) P2[t] = p*P[t-1];

    P2[k] = P2[k] + p*P2[k-1];

    for(t = 0; t<=k; t++) P[t] = P2[t];
  }

  return P2[k];
}

void estimate_meanIBS()
{
    int n=0;
    int i, j;
    double varp0=0, varp2=0;

  for (i = 1; i <= nochrom; i++)
    for (j = 1; j <= nomark[i]; j++) {

      double p0=0, p2=0;

      for (int a1 = 1; a1 <= noalle[i][j]; a1++)
      for (int a2 = 1; a2 <= noalle[i][j]; a2++) {

          double PA = q[i][j][a1]*q[i][j][a2];

          for (int b1 = 1; b1 <= noalle[i][j]; b1++) if(b1!=a1 && b1!=a2)
          for (int b2 = 1; b2 <= noalle[i][j]; b2++) if(b2!=a1 && b2!=a2)
              p0 += q[i][j][a1]*q[i][j][a2]*q[i][j][b1]*q[i][j][b2];

          for (int b1 = 1; b1 <= noalle[i][j]; b1++)
          for (int b2 = 1; b2 <= noalle[i][j]; b2++)
              if( (b1==a1 && b2==a2) || (b1==a2 && b2==a1))
                  p2 += q[i][j][a1]*q[i][j][a2]*q[i][j][b1]*q[i][j][b2];

          n++;
      }

      pIBS0 += p0;
      pIBS2 += p2;
      varp0 += p0*(1-p0);
      varp2 += p2*(1-p2);
    }

    int N = nomark[0];

    pIBS0 /= nomark[0];
    pIBS2 /= nomark[0];

    varp0 = varp0/(N*N);
    varp2 = varp2/(N*N);


    printf("Expected IBS for unrelated pairs: p.IBS0 = %6.4f  p.IBS1 = %6.4f  p.IBS2 = %6.4f\n",
        pIBS0,1-pIBS0-pIBS2,pIBS2);

    //cout << pMarkov(5000, 100, 0.95) << "\n";
    
    double p = 1-pIBS0;
    
    int r = 0;
    for(r = 1; r < N; r++) {
     double perr = 1 - pow( 1- pow(p,r), N-r );
     //double pm = pMarkov(N,r,p);
     //printf("%d %.2g %.2g %f\n",r,perr, pm , log(perr/pm) );
     if(perr<1e-5) break;
    }
    pIBS0_th = r;


    p = pIBS2;
    for(r = 1; r < N; r++) {
     double perr = 1 - pow( 1- pow(p,r), N-r );
     //double pm = pMarkov(N,r,p);
     //printf("%d %.2g %.2g %f\n",r,perr, pm , log(perr/pm) );
     if(perr<1e-5) break;
    }
    
    pIBS2_th = r;
    if(debug) printf("Thresholds for defining regions as IBD1 = %.0f, IBD2 = %.0f\n",pIBS0_th,pIBS2_th);
}


int main(int argc, char *argv[])
{
  int i, j, k, l, kk, kkk;
  int bugexit, numread, maxpedsize;
  int lineid;

  int id1 = 0, id2 = 0, idped1, idped2; // for testing a specific pair
 
  /* for a pair */
  int reltype, findit, index1 = 0, index2 = 0;
  /* statistic, mu, std, z_obs, p-value(2-sided) */

  /* genodata of all the members for a given pedigree */
  geno_t **geno;
  geno_t **genoB;

  FILE *fp_in1;
  char str1[100], str2[100], str3[100], str4[100];

  buf_len = 64000;
  buf = (char *) malloc( buf_len);

  stringmap_init();
  stringmap_lookup("0"); /* add '0', so that '0' maps to integer 0 (for missing data) */

  printf(
  "-------------------------------------------------------------\n"
  "-    PREST-Plus (Pedigree RElationship Statistical Test)    -\n"
  "-            version 4.09 - October 1st, 2010               -\n"
  "-------------------------------------------------------------\n\n"
  );


  /***** Reading of command line options ****/

  options opts(argc, argv);

  if (opts.get("--config", 1)) read_config_file(opts.arg(0));
  
  if (opts.get("--file", 1)) {
      nochrom = 1;
      strcpy( ped_name, opts.arg(0));
      strcpy( chromnameped[1], opts.arg(0));
  }

   if (opts.get("--geno", 1)) {
      nochrom = 1;
      strcpy( chromnameidx[1] , "-");
      strcpy( chromnameped[1], opts.arg(0));
  }

  if (opts.get("--ped", 1)) {
      nochrom = 1;
      strcpy( ped_name, opts.arg(0));
  }

  if (opts.get("--map", 1)) {
      strcpy( chromnameidx[1], opts.arg(0));
      ALLELEFREQ = 1;
  }

  if (opts.get("--freq", 1)) {
      strcpy( chromnameidx[1], opts.arg(0));
      ALLELEFREQ = 1;
  }

  if (opts.get("--out", 1)) strcpy(output_prefix, opts.arg(0) );
  
  if (opts.get("--cm", 0)) MAP = 1;
  if (opts.get("--uniform",0)) uniformSPC = 1;
  if (opts.get("--uniformfreq",0)) setUniformFreq = 1;

  if (opts.get("--mlrt", 0)) OPTION = 2;
  if (opts.get("--mlrtrep", 1)) { NREP = atoi(opts.arg(0)); ALPHA = 0; }

  if (opts.get("--wped", 0)) WPED = 1;
  if (opts.get("--aped", 0)) APED = 1;
  if (opts.get("--newton", 1)) ibd_useNewton = atoi(opts.arg(0));
  if (opts.get("--all", 0)) all_pairs = 1;

  if (opts.get("--debug", 1)) debug = atoi(opts.arg(0));
  if (opts.get("--eps", 1)) MYEPS = atof(opts.arg(0));

  if (opts.get("--miter", 1)) MITER = atoi(opts.arg(0));

  if (opts.get("--error.rate", 1)) errorrate = atof(opts.arg(0));

  if (opts.get("--sim", 2)) { _sim_rel = atoi(opts.arg(0)); _sim_n = atoi(opts.arg(1)); }
  
  if (opts.get("--writefreq", 0)) write_frequency_file = 1;


  if (opts.get("--pair", 4)) {
      idped1 = stringmap_lookup(opts.arg(0));
      id1 = stringmap_lookup(opts.arg(1));
      idped2 = stringmap_lookup(opts.arg(2));
      id2 = stringmap_lookup(opts.arg(3));
  }
    
 if (opts.get("--reltype", 1)) alternativeRel = atoi(opts.arg(0));

  if ( nochrom == 0 ) {
    printf("\nTypical usage: prest --geno data.ped --map data.map --ped pedigree --wped\n\n");
    printf("See documentation for an explanation of options.\n");
    exit(1);
  }
  
    if (opts.get("--join", 1)) {
          printf("Joining the results of a cluster run..\n");
          cluster_join = atoi(opts.arg(0));
          cluster_id = -999; // signifies to join files
    }

    if (opts.get("--cluster", 1)) {
        int n = sscanf(opts.arg(0), "%d/%d", &cluster_id, &cluster_n);

        if(n != 2) { printf("Malformed cluster string\n"); exit(0); }

        printf("\nRunning prest on a cluster. %d out of %d nodes. \n",cluster_id,cluster_n);
        cluster_id --;
    }
 
  if (strcmp(chromnameidx[1],"") == 0) {
    printf("Error: No map file specified. Use option --map to specify the map file.\n");
    exit(1);
  }

  if(strcmp(ped_name,"(nofile)") == 0) {
    strcpy(ped_name, chromnameped[1]);
  }

  if ( (fp_in1 = fopen(ped_name,"r")) == 0 ) {
    printf("\nUnable to open input pedigree file: %s\n\n", ped_name);
    exit(0);
  }



  printf("... Input files ... \n");
  printf("Pedigree file : %s (--ped) \n", ped_name );
  printf("Genotype file : %s (--file) \n",chromnameped[1] );
  printf("     Map file : %s (--map)\n",chromnameidx[1] );
  printf("\n\n");
  printf("...  Output files ...\n");
  printf("    Main results file : %s_results (--out)\n",output_prefix );
  printf(" Mendelian error file : %s_mendel\n",output_prefix );
  printf("\n\n");

  printf(
"Relationship types considered:\n\n"
"NUM          NAME           CODE    KINSHIP     E(IBD0)  E(IBD1)  E(IBD2)\n"
"-------------------------------------------------------------------------\n"
" 1       Full-Sibling       FS      0.250       0.25     0.50     0.25\n"
" 2       Half-Sibling       HS      0.125       0.50     0.50     0\n"
" 3 Grandparent-grandchild   GPC     0.125       0.50     0.50     0\n"
" 4        Avuncular         AV      0.125       0.50     0.50     0\n"
" 5       First-Cousin       FC      0.0625      0.75     0.25     0\n"
" 6        Unrelated         U       0           1        0        0\n"
" 7      Half-Avuncular      HA      0.0625      0.75     0.25     0\n"
" 9  Half-Sib+First-Cousin   HSFC    0.1875      0.375    0.5      0.125\n"
"10    Parrent-Offspring     PO      0.25        0        1        0\n"
"11         MZ-Twins         MZ      0.5         0        0        1\n"
"\n\n"  );
  
  /* -------------------------------------------------------------- */


  /* -------------------------------------------------------------- */
  open_output_files();

  if(cluster_id == -999) { join_cluster_files(); exit(0); }


  bugexit = 0;
  maptype = 2;
  int menderrors = 0;

  srand((long)time(NULL));
  if(debug) srand(0);

  /* -------------------------------------------------------------- */

  /* Reading of pedigree file */

  /* noped: no. of the pedigrees in the file */
  int itemp[10];

  noped = 0;
  pedid_now = -1;

  lineid = 1;

  while ( fgets_adj(&buf, &buf_len, fp_in1) ) {

    numread = sscanf(buf, "%s %s %s %s %d %d", str1, str2,
		     str3, str4, &itemp[5], &itemp[6]);

    if ( numread != 6 ) {
      printf("\nInvalid line %d in file %s\n", lineid, ped_name);
      printf("Correct format is: pedid id f_id m_id sex affect_status\n\n");
      exit(0);
    }

    itemp[1] = stringmap_lookup(str1);
    itemp[2] = stringmap_lookup(str2);
    itemp[3] = stringmap_lookup(str3);
    itemp[4] = stringmap_lookup(str4);

    /* Allow -9 for missing sex and affection status */
    if(itemp[5] == -9) itemp[5] = 0;
    if(itemp[6] == -9) itemp[6] = 0;

    if(max_fid < strlen(str1)) max_fid = strlen(str1);
    if(max_fid < strlen(str3)) max_fid = strlen(str3);
    if(max_iid < strlen(str2)) max_iid = strlen(str2);
    if(max_iid < strlen(str4)) max_iid = strlen(str4);

    if ( itemp[1] != pedid_now ) {
      ++noped;
      pedid_now = itemp[1];
    }
    ++lineid;
  }

  nomem = new_int2d(noped+1, 2);
  nomem[noped+1][2] = 0;

  /* nomem: pedid, size of the ped*/
  rewind(fp_in1);
  pedid_now = -1;
  numpedid = 0;

  while ( fgets_adj(&buf, &buf_len, fp_in1) ) {

    numread = sscanf(buf, "%s %s %s %s %d %d", str1, str2,
		     str3, str4, &itemp[5], &itemp[6]);

    itemp[1] = stringmap_lookup(str1);
    itemp[2] = stringmap_lookup(str2);
    itemp[3] = stringmap_lookup(str3);
    itemp[4] = stringmap_lookup(str4);

    /* Allow -9 for missing sex and affection status */
    if(itemp[5] == -9) itemp[5] = 0;
    if(itemp[6] == -9) itemp[6] = 0;


    if (itemp[1] != pedid_now) {
      ++numpedid;
      nomem[numpedid][1] = pedid_now = itemp[1];
      nomem[numpedid][2] = 1;
    }
    else {
      ++nomem[numpedid][2];
    }
  }

  ped.resize(noped+2);

  for (i = 1; i <= noped; i++) 
    nomem[noped+1][2] += nomem[i][2];

  /* ped: id, fatherid, motherid, recoded fid, recoded mid, gender */
  rewind(fp_in1);

  for (i = 1; i <= noped; i++) {

    ped[i].resize(nomem[i][2]+1);
    
    for(j = 0; j < 6; j++) ped[i][0][j] = 0;
    
    for (j = 1; j <= nomem[i][2]; j++) {
        int sex, aff;

	fgets_adj(&buf, &buf_len, fp_in1);
        numread = sscanf(buf, "%s %s %s %s %d %d", str1, str2,
                         str3, str4, &sex, &aff);

        if(sex == -9) sex = 0;
        if(aff == -9) aff = 0;

        itemp[1] = stringmap_lookup(str1);
        ped[i][j][1] = stringmap_lookup(str2);
        ped[i][j][2] = stringmap_lookup(str3);
        ped[i][j][3] = stringmap_lookup(str4);
        ped[i][j][6] = sex;
        itemp[6] = aff;
    }
  }

  fclose(fp_in1);


  /* -------------------------------------------------------------- */

  /* recoded parents's ids and errors check */
  for (i = 1; i <= noped; i++) {
    for (j = i+1; j <= noped; j++) {
      if ( nomem[i][1] == nomem[j][1] ) {
	printf("Error: Pedigree id %s is used twice (make sure your pedigree file is sorted by pedigree).\n", s(nomem[i][1]));
	++bugexit;
      }
    }

    bugexit += recoding(ped[i], nomem[i]);
  }

  if ( bugexit != 0 ) {
    printf("\nErrors in the pedigree file found. Prest will now exit.\n");
    exit(1);
  }


  /* -------------------------------------------------------------- */

  /* find all the relationship types considered
     1=full-sib, 2=half-sib, 3=grandp-c, 4=avuncular, 5=first-cousin,
     6=unrelated, 7=half-avuncular, 8=half-first-cousin,
     9=half-sib+first-cousin,
     10 = parent-offspring */

  rel = (int ***)malloc((size_t) ((noped+1)*sizeof(int**)));
  for (i = 1; i <= noped; i++) {
    kk = nomem[i][2]*(nomem[i][2]-1)/2;
    rel[i] = new_int2d( kk+1, 3);
  }

  norel = new_int2d(noped+2, NR+2 );

  for (i = 1; i <= (noped+1); i++)
    for (j = 1; j <= (NR+1); j++)
      norel[i][j] = 0;


  for (i = 1; i <= noped; i++) {
    /* input: i_th ped, size
       return: no. of the pairs found, and
       rel: id1, id2, reltype */
    norel[i][NR+1] = findrelation(ped[i], nomem[i][2], rel[i]);
    norel[noped+1][NR+1] += norel[i][NR+1];
    for (j = 1; j <= norel[i][NR+1]; j++) {
      for (k = 1; k <= NR; k++) {
	if ( rel[i][j][3] == k )
	  ++norel[i][k];
      }
    }
  }


  /* -------------------------------------------------------------- */
  /* Reading of map files */
  int have_freq = read_map_files();

  if(MAP == 0 && OPTION == 2) printf("Warning: No genetic map. Cannot perform MLRT test.\n");

  if(have_freq == 0) compute_allele_frequencies();

  /*if(MAP == 1) {
    if(ALLELEFREQ == 1) {
      read_map_files(1);
    }
    else {
      compute_allele_frequencies();
      read_spacings();
    }
  } else {
    if(ALLELEFREQ == 1)
      read_map_files(0);
     else
      compute_allele_frequencies();
  }*/


  if(debug)
  {
    for (int i = 0; i <= nochrom; i++)
      printf("Chromosome %3d length = %f\n",i, chromlength[i]);
  }

  estimate_meanIBS();
  if (opts.get("--pp", 1)) pIBS0_th = atof(opts.arg(0));

  // Index genotype files for faster reading
  indexgenotypes();

  /* -------------------------------------------------------------- */

  /* Prepare Variance of EIBD, AIBS, and IBS */

  /* ngeno: no. of possible genotypes for each marker = n(n+1)/2
     genotype: detailed possible genotype (g1, g2) for each
     possible genotype for each marker */
  ngeno = new_int1d(nomark[0]);

  genotype = (int ***)malloc((size_t) ((nomark[0]+1)*sizeof(int**)));

  kkk = 0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      ngeno[kkk+j] = noalle[i][j]*(noalle[i][j]+1)/2;
      genotype[kkk+j] = (int **)malloc((size_t)
				       ((ngeno[kkk+j]+1)*sizeof(int*)));
      for (k = 1; k <= ngeno[kkk+j]; k++)
	genotype[kkk+j][k] = (int *)malloc((size_t) ((2+1)*sizeof(int)));
    }
    kkk += nomark[i];
  }

  kkk = 0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      kk = 1;
      for (k = 1; k <= noalle[i][j]; k++) {
	for (l = k; l <= noalle[i][j]; l++) {
	  genotype[kkk+j][kk][1]=k;
	  genotype[kkk+j][kk][2]=l;
	  ++kk;
	}
      }
    }
    kkk += nomark[i];
  }

  /* for EM IBD calculations */

  a = (double *)malloc((size_t) ((nomark[0]+1)*sizeof(double)));
  b = (double *)malloc((size_t) ((nomark[0]+1)*sizeof(double)));
  c = (double *)malloc((size_t) ((nomark[0]+1)*sizeof(double)));

  /* for EIBD variance */

  eeibdibd = (double ***)malloc((size_t) ((NR+1)*sizeof(double**)));
  eeibdeibdibd = (double ***)malloc((size_t) ((NR+1)*sizeof(double**)));
  eaibsibd = (double ***)malloc((size_t) ((NR+1)*sizeof(double**)));
  eaibsaibsibd = (double ***)malloc((size_t) ((NR+1)*sizeof(double**)));

  eibsibd = new_double2d(nomark[0], 3);
  eibsibsibd = new_double2d(nomark[0], 3);

  for (i = 1; i <= NR; i++) {
	  eeibdibd[i] = new_double2d(nomark[0],3);
	  eaibsibd[i] = new_double2d(nomark[0],3);
	  eeibdeibdibd[i] = new_double2d(nomark[0],3);
	  eaibsaibsibd[i] = new_double2d(nomark[0],3);
  }

  forEIBDvar(eeibdibd, eeibdeibdibd, nochrom, nomark,
	     ngeno, genotype, q, p012);

  /* for AIBS variance */
  forAIBSvar(eaibsibd, eaibsaibsibd, nochrom, nomark,
	     ngeno, genotype, q, p012);

  /* for IBS variance
     notice that it does not depend on reltype
     eibsibd[i][1,2,3]: E(IBS at locus i|ibd=0) E(IBS_i|ibd=1)
     E(IBS_i|ibd=2) */

  forIBSvar(eibsibd, eibsibsibd, nochrom, nomark, q, noalle);


  /* -------------------------------------------------------------- */

  /* for each pedigree, calculat EIBD, AIBS, IBS statistics,
     mu, and std. Use Normal approximation to obtain 2-sided p-value.
     If users want to perform MLRT test, further simulate genotype
     data for each pair, and calculate empirical p-value for EIBD,
     AIBS, IBS and MLRT */

  valid = new_int1d( nomark[0]);
  ibs1 = new_int1d( nomark[0]);
  ibs2 = new_int1d( nomark[0]);
  
  norelwdata = new_int2d( noped+1, NR+1 );
  nopairsvalid = new_int2d( noped+1, NR+1 );

  for (i = 0; i <= (noped+1); i++)
    for (j = 0; j <= (NR+1); j++) {
      norelwdata[i][j] = 0;
      nopairsvalid[i][j] = 0;
    }


  maxpedsize = nomem[1][2];
  for (i = 2; i <= noped; i++) if( nomem[i][2] > maxpedsize) maxpedsize = nomem[i][2];


  geno = new_genot2d( maxpedsize, 2 * nomark[0]);

  if(APED == 1 || id1 > 0)
    genoB = new_genot2d( maxpedsize, 2 * nomark[0]);
  else
    genoB = geno;


  simugeno1 =  (geno_t *)malloc((size_t)  ((2*nomark[0]+1)*sizeof(geno_t)));
  simugeno2 =  (geno_t *)malloc((size_t)  ((2*nomark[0]+1)*sizeof(geno_t)));


  max_fid += 1;
  max_iid += 1;

  mmark = (int)log10((double)nomark[0]) + 1; /* max field size of mmark */
  
  if(mmark < 7) mmark = 7;

  char   ctemp1[MAXSTR+1];

  sprintf(ctemp1, "%%%ds %%%ds %%%ds %%%ds %%7s %%%ds ",
		  max_fid, max_iid, max_fid, max_iid, mmark);

  ///fprintf(fp_warnings, ctemp1, "FID1", "IID1", "FID2", "IID2", "reltype", "commark");
  ///fprintf(fp_warnings, "    p012_1   p012_2   p012_3    p.IBD0   p.IBD1   p.IBD2\n");

  fprintf(fp_results,  ctemp1, "FID1", "IID1", "FID2", "IID2", "reltype", "commark");


  fprintf(fp_results, "   EIBD   AIBS    IBS");
  fprintf(fp_results, " p.IBS0 p.IBS1 p.IBS2");
  
  if(stat_ibd == 1)
    fprintf(fp_results, " p.IBD0 p.IBD1 p.IBD2");


  if(MAP == 1)
    fprintf(fp_results, " p.normal.EIBD p.normal.AIBS p.normal.IBS");

  if (OPTION == 2)
    fprintf(fp_results, " p.empirical.EIBD  p.empirical.AIBS  p.empirical.IBS  p.empirical.MLRT");

  if(stat_ibs == 1) {
    fprintf(fp_results, "  %8s", "RMAX1");
    fprintf(fp_results, " %8s", "RAVG1");
    // fprintf(fp_results, " %8s", "RIBD1");
    fprintf(fp_results, "  %8s", "RMAX2");
    fprintf(fp_results, " %8s", "RAVG2");
    // fprintf(fp_results, " %8s", "RIBD2");
  }

  fprintf(fp_results, "\n");

  /* END: Printing of file headers */

  // For running in a cluster, only test every n-th pair
  int count = 0;


  // Test a specific pair of individuals 
  if(id1 > 0) {

      int ped1 = 0, ped2 = 0, j1 = 0, j2 = 0;

      for (i = 1; i <= noped; i++)
        if(nomem[i][1] == idped1)
          for (j = 1; j <= nomem[i][2]; j++) if(ped[i][j][1] == id1) { ped1 = i; j1 = j; }
      
      for (i = 1; i <= noped; i++) 
        if(nomem[i][1]==idped2)
          for (j = 1; j <= nomem[i][2]; j++) if(ped[i][j][1] == id2) { ped2 = i; j2 = j; }
      
   
      if(j1 == 0 || j2 == 0) { printf("Pair not found\n"); exit(1); }

      bugexit += findpedgeno(ped[ped1], nomem[ped1], nomark, noalle, geno);

      if(ped1 == ped2)
          genoB = geno;
      else
         bugexit += findpedgeno(ped[ped2], nomem[ped2], nomark, noalle, genoB);

      int reltype = 6;
      
      if(ped1 == ped2) 
          for (j = 1; j <= norel[ped1][NR+1]; j++) {
              int t1 = rel[ped1][j][1], t2 = rel[ped2][j][2];
              if(t1==id1 && t2==id2) reltype = rel[ped1][j][3];
              if(t1==id2 && t2==id1) reltype = rel[ped1][j][3];
          }
        printf("\n* Testing pair: %s %s - %s %s\n", s(idped1), s(id1), s(idped2), s(id2));

	for(reltype=1;reltype<=10;reltype++)
        for (int index1 = 1; index1 <= nomem[ped1][2]; index1 ++) if(ped[ped1][index1][1] == id1)
          for (int index2 = 1; index2 <= nomem[ped2][2]; index2 ++) if(ped[ped2][index2][1] == id2) {
              //printf("ped=%d %d index=%d %d\n",ped1,ped2,index1,index2);
                    test_pair(geno, ped1, index1, genoB, ped2, index2, reltype);
                    printf("%s %s - %s %s reltype=%2d p.MLRT=%f\n", s(idped1), s(id1), s(idped2), s(id2), reltype, pMLRT);

          }


      exit(0);
  }

  /* File for output of mendelian errors */
  string fp_mend_file = ((string)output_prefix+"_mendel");
  FILE *fp_mend = fopen(  fp_mend_file.c_str(),"w");

  printf("Checking for mendelian errors (output file %s).\n",fp_mend_file.c_str());

  fprintf(fp_mend,"FID\tIID\tMARKER\tGENO\tFATHER\tMOTHER\n");
  
  /* Check within pedigrees */

  if(WPED == 0 && APED == 0) {
    for (i = 1; i <= noped; i++) {

      /* find geno data for a given pedigree */
      bugexit += findpedgeno(ped[i], nomem[i], nomark, noalle, geno);

      /* check Mendelian errors between parent-child */
      menderrors += checkMend(fp_mend, ped[i], geno, nomem[i], nomark);
    }
  }

  printf("Number of mendelian errors: %d\n\n",menderrors);


  if(_sim_rel > 0) {

    printf("\n Simulating %d pairs of reltype %d\n",_sim_n, _sim_rel);
    
    int reltype = _sim_rel;
    char fname[500];
    sprintf(fname,"%s.ped",output_prefix);
    FILE *out = fopen(fname, "w");

    for (k = 1; k <= nomark[0]; k++) valid[k] = 0;

    for(int jj=0; jj<_sim_n; jj++) {

    simugeno1geno2(q, qcum, nochrom, nomark, noalle, reltype,
      chromlength, cenm, simugeno1, simugeno2);

    simugeno1[0] = geno[index1][0];
    simugeno2[0] = genoB[index2][0];

    {
       printf("simugeno: %d\n", reltype);

       fprintf(out, "SIM%d_%d P1 0 0 0 %d ", reltype, _sim_pedno, reltype);
       {for(int i=1;i<=2*nomark[0];i++) fprintf(out, "%d ", simugeno1[i]); fprintf(out, "\n");}

       fprintf(out, "SIM%d_%d P2 0 0 0 %d ", reltype,_sim_pedno, /*s(nomem[j][1]), s(genoB[index2][0]),*/reltype );

       _sim_pedno ++;
       {for(int i=1;i<=2*nomark[0];i++) fprintf(out, "%d ", simugeno2[i]); fprintf(out, "\n");}

       double phat[20];

        EMforp012(simugeno1,simugeno2, valid,
                  q, nochrom, nomark, reltype, phat);

        printf("EM: %f %f %f \n",phat[1],phat[2],phat[3]);

    }

    }

    exit(0);
  }



  if (WPED == 1) {
    printf("Within pedigree testing\n");
    printf("-----------------------\n");
  }

  if (WPED == 1)
    for (i = 1; i <= noped; i++) {
      if(norel[i][NR+1] == 0) continue;
      printf("Testing pedigree %s. Pairs to be tested: %d\n", s(nomem[i][1]), norel[i][NR+1]); fflush(stdout);

      /* find geno data for a given pedigree */
      bugexit += findpedgeno(ped[i], nomem[i], nomark, noalle, geno);

      /* check Mendelian errors between parent-child */
     // menderrors += checkMend(fp_mend,ped[i], geno, nomem[i], nomark);

      //printf("Pairs to be tested: %d\n",norel[i][NR+1]);

      /* for each pair */
      for (j = 1; j <= norel[i][NR+1]; j++) {
	/* Here i and j are the pair to be tested */

	reltype = rel[i][j][3];

	findit = 0;
	k = 1;
     
	while (findit != 2) {
	  if (geno[k][0] == rel[i][j][1]) {
	    index1 = k;
	    ++findit;
	  }

	  if (geno[k][0] == rel[i][j][2]) {
	    index2 = k;
	    ++findit;
	  }
	  ++k;
	}

        int runit = 1;
        count ++;
        if(cluster_n > 0 && (count % cluster_n) != cluster_id) runit = 0;

	if(runit) test_pair(geno, i, index1, geno, i, index2, reltype);
      }


    }

    if (APED == 1) {
     printf("\nAcross pedigree testing\n");
     printf(  "-----------------------\n");
    }

    if(APED == 1)
    for (i = 1; i <= noped; i++) {
      /* find geno data for a given pedigree */
      bugexit += findpedgeno(ped[i], nomem[i], nomark, noalle, geno);


      /* check Mendelian errors between parent-child */
      if(WPED==0) menderrors += checkMend(fp_mend,ped[i], geno, nomem[i], nomark);

      for(j = i+1; j<= noped;j++) {
        printf("Testing pedigree %s x %s.\n", s(nomem[i][1]), s(nomem[j][1])); fflush(stdout);

        /* find geno data for a given pedigree */
        bugexit += findpedgeno(ped[j], nomem[j], nomark, noalle, genoB);

        /* for each pair ACROSS PEDIGREES*/
        for (index1 = 1; index1 <= nomem[i][2]; index1 ++)
          for (index2 = 1; index2 <= nomem[j][2]; index2 ++) {

              int runit = 1;
              count ++;
              if(cluster_n > 0 && (count % cluster_n) != cluster_id) runit = 0;

	      if(runit) test_pair(geno, i, index1, genoB, j, index2, 6);
          }
        }

      }

  printf("\n");

  fclose(fp_mend);

  if (menderrors != 0) {
    printf("Warning: %d mendelian errors detected.\n", menderrors);

    printf("The errors are reported in file %s_mendel.\n",output_prefix);
    printf("In the current implementation of prest, the presence of Mendelian errors ");
    printf("is not taken into account in the hypothesis testing.\n\n");
  }


  
  for (k = 1; k <= 3; k++) {
    for (j = 1; j <= NR; j++) {
      checkallefreq[k][j][1] =
	checkallefreq[k][j][1]/norelwdata[noped+1][j];
      checkallefreq[k][j][2] =
	checkallefreq[k][j][2]/norelwdata[noped+1][j];
      checkallefreq[k][j][4] =
	checkallefreq[k][j][4]/norelwdata[noped+1][j];
      checkallefreq[k][j][3] =
	sqrt((checkallefreq[k][j][3] -
	      norelwdata[noped+1][j]*checkallefreq[k][j][2]*checkallefreq[k][j][2])/(norelwdata[noped+1][j]-1));
      checkallefreq[k][j][5] =
	sqrt((checkallefreq[k][j][5] -
	      norelwdata[noped+1][j]*checkallefreq[k][j][4]*checkallefreq[k][j][4])/(norelwdata[noped+1][j]-1));
    }
  }
  for (j = 1; j <= NR; j++) {
    checkmarker[j][1] = checkmarker[j][1]/norelwdata[noped+1][j];
    checkmarker[j][2] = sqrt((checkmarker[j][2] -
			      norelwdata[noped+1][j]*checkmarker[j][1]*checkmarker[j][1])/(norelwdata[noped+1][j]-1));
  }

  /* delete later */
  /*  if (flag == 2) { */
/*      for (i = 1; i <= NR; i++) { */
/*        checkmu[i] = checkmu[i]/(nopairsvalid[noped+1][i]*NREP); */
/*        printf("%f\n",checkmu[i]); */
/*      } */
/*    } */


  /* -------------------------------------------------------------- */

  /* output general information */

  fprintf(fp_summary, "# Number of pedigrees: \n");
  fprintf(fp_summary, "    %d\n\n", noped);

  fprintf(fp_summary, "# Total number of individuals: \n");
  fprintf(fp_summary, "    %d\n\n", nomem[noped+1][2]);

  fprintf(fp_summary, "# Number of chromosomes: \n");
  fprintf(fp_summary, "    %d\n\n", nochrom);

  fprintf(fp_summary, "# Total number of markers: \n");
  fprintf(fp_summary, "    %d\n\n", nomark[0]);

  fprintf(fp_summary, "# Number of markers on each chromosome\n\n");
  for (i = 1; i <= nochrom; i++)
    fprintf(fp_summary, "    %s: %5d \n", chromnameidx[i], nomark[i]);
  fprintf(fp_summary, "\n");

  fprintf(fp_summary, "# Number of genotyped pairs found in each of the\n");
  fprintf(fp_summary, "# relationship categories analyzed by PREST: \n");
  fprintf(fp_summary, "# (Note that unrelated pairs are checked within each pedigree.)\n\n");

  fprintf(fp_summary, "    Total:  %d\n\n", norelwdata[noped+1][NR+1]);
  fprintf(fp_summary, "    %7d full-sib pairs\n",
	  norelwdata[noped+1][1]);
  fprintf(fp_summary, "    %7d half-sib pairs \n",
	  norelwdata[noped+1][2]);
  fprintf(fp_summary, "    %7d grandparent-grandchild pairs\n",
	  norelwdata[noped+1][3]);
  fprintf(fp_summary, "    %7d avuncular pairs\n",
	  norelwdata[noped+1][4]);
  fprintf(fp_summary, "    %7d first-cousin pairs\n",
	  norelwdata[noped+1][5]);
  fprintf(fp_summary, "    %7d unrelated pairs\n",
	  norelwdata[noped+1][6]);
  fprintf(fp_summary, "    %7d half-avuncular pairs\n",
	  norelwdata[noped+1][7]);
  fprintf(fp_summary, "    %7d half-first-cousin pairs\n",
	  norelwdata[noped+1][8]);
  fprintf(fp_summary, "    %7d half-sib-plus-first-cousin pairs\n",
	  norelwdata[noped+1][9]);
  fprintf(fp_summary, "    %7d parent-offspring pairs\n\n",
	  norelwdata[noped+1][10]);

  fprintf(fp_summary, "# Summary of genotyped pairs found in each pedigree\n\n");
  fprintf(fp_summary, "    pedid size f.sib h.sib g.p.c avun. f.cous. unrel. h.avun h.f.cous. h.sib+f.cous. p.o.   Total\n");
  for (i = 1; i <= noped; i++) {
    fprintf(fp_summary, " %8s %4d ", s(nomem[i][1]), nomem[i][2]);
    for (j = 1; j <= NR; j++)
      fprintf(fp_summary, " %4d", norelwdata[i][j]);
    fprintf(fp_summary, " %8d\n", norelwdata[i][NR+1]);
  }
  fprintf(fp_summary, "\n");


  if( OPTION == 1) {
    fprintf(fp_summary,"# Number of pairs that have p-value of either EIBD or AIBS less than 0.2\n");
    fprintf(fp_summary,"# (for unrelated, p-value of IBS less than 0.2).\n");
    fprintf(fp_summary,"# This is the number of pairs for which PREST would perform the MLRT.\n");
    fprintf(fp_summary,"# (The MLRT will actually be performed only if option 2 is chosen.)\n");
    fprintf(fp_summary,"# For each pair, performance of the MLRT, including 100,000\n");
    fprintf(fp_summary,"# simulated realizations takes approximately 5 minutes on \n");
    fprintf(fp_summary,"# a Sun Ultra II with 360-MHz processor.\n\n");
    fprintf(fp_summary,"    %d\n\n", countmlrt);
  }


  if ( OPTION == 2 ) {
    fprintf(fp_summary, "# Number of pairs checked by MLRT\n\n");
    fprintf(fp_summary, "    Total:  %d\n\n", nopairsvalid[noped+1][NR+1]);
    fprintf(fp_summary, "    %7d full-sib pairs\n",
	    nopairsvalid[noped+1][1]);
    fprintf(fp_summary, "    %7d half-sib pairs\n",
	    nopairsvalid[noped+1][2]);
    fprintf(fp_summary, "    %7d grandparent-grandchild pairs\n",
	    nopairsvalid[noped+1][3]);
    fprintf(fp_summary, "    %7d avuncular pairs\n",
	    nopairsvalid[noped+1][4]);
    fprintf(fp_summary, "    %7d first-cousin pairs\n",
	    nopairsvalid[noped+1][5]);
    fprintf(fp_summary, "    %7d unrelated pairs\n",
	    nopairsvalid[noped+1][6]);
    fprintf(fp_summary, "    %7d half-avuncular pairs\n",
	    nopairsvalid[noped+1][7]);
    fprintf(fp_summary, "    %7d half-first-cousin pairs\n",
	    nopairsvalid[noped+1][8]);
    fprintf(fp_summary, "    %7d half-sib-plus-first-cousin pairs\n",
	    nopairsvalid[noped+1][9]);
    fprintf(fp_summary, "    %7d parent-offspring pairs\n\n",
	    nopairsvalid[noped+1][10]);

    fprintf(fp_summary, "    In each pedigree\n\n");
    fprintf(fp_summary, "    pedid size f.sib h.sib g.p.c avun. f.cous. unrel. h.avun h.f.cous. h.sib+f.cous. p.o.   Total\n");
    for (i = 1; i <= noped; i++) {
      fprintf(fp_summary, " %8s %4d ", s(nomem[i][1]), nomem[i][2]);
      for (j = 1; j <= NR; j++)
	fprintf(fp_summary, " %4d", nopairsvalid[i][j]);
      fprintf(fp_summary, " %8d\n", nopairsvalid[i][NR+1]);
    }
    fprintf(fp_summary, "\n");
  }

  fprintf(fp_summary, "# Assuming that there are not many misclassified relative pairs in the\n");
  fprintf(fp_summary, "# data set and assuming that all the pairs are typed on roughly the\n");
  fprintf(fp_summary, "# same markers (because which markers are typed affects the null means\n");
  fprintf(fp_summary, "# of the AIBS and IBS statistics and affects the null variances of the\n");
  fprintf(fp_summary, "# EIBD, AIBS and IBS statistics), then the sample averaged observed\n");
  fprintf(fp_summary, "# statistics for each relationship type should not be significantly\n");
  fprintf(fp_summary, "# different from what is expected.  If they are significantly\n");
  fprintf(fp_summary, "# different, if may indicate that the allele frequencies are not\n");
  fprintf(fp_summary, "# accurately estimated, or that the assumptions of population\n");
  fprintf(fp_summary, "# homogeneity and Hardy-Weinberg equilibrium are not appropriate.\n");
  fprintf(fp_summary, "# It could also indicate a problem with the way the data are coded \n");
  fprintf(fp_summary, "# (see item #7 in 'Tips' file).  The following output contains \n");
  fprintf(fp_summary, "# the information of sample summary of the statistics\n\n");

  fprintf(fp_summary, "# obs      = sample averaged observed statistic\n");
  fprintf(fp_summary, "# mu(std)  = sample averaged null mean (sample std. dev. of null mean)\n");
  fprintf(fp_summary, "#            (Note that the null mean of EIBD depends only on\n");
  fprintf(fp_summary, "             the null relationship, not on the data, so for EIBD \n");
  fprintf(fp_summary, "#            we give the null mean in place of the sample averaged null mean\n");

  fprintf(fp_summary, "# std(std) = sample averaged standard deviation\n");
  fprintf(fp_summary, "#            (sample std. dev. of null standard deviation)\n");
  fprintf(fp_summary, "# no.mark  = sample averaged number of typed markers for a pair\n");
  fprintf(fp_summary, "#            (sample std. dev.)\n\n");

  fprintf(fp_summary, "              no.pairs no.mark(std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_summary, "    full-sib   %7d      NA(    NA)\n",
	    norelwdata[noped+1][1]);
  else
    fprintf(fp_summary, "    full-sib   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][1], checkmarker[1][1],
	    checkmarker[1][2]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_summary, "    half-sib   %7d      NA(    NA)\n",
	    norelwdata[noped+1][2]);
  else
    fprintf(fp_summary, "    half-sib   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][2], checkmarker[2][1],
	    checkmarker[2][2]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_summary, "    grandp-c   %7d      NA(    NA)\n",
	    norelwdata[noped+1][3]);
  else
    fprintf(fp_summary, "    grandp-c   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][3], checkmarker[3][1],
	    checkmarker[3][2]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_summary, "   avuncular   %7d      NA(    NA)\n",
	    norelwdata[noped+1][4]);
  else
    fprintf(fp_summary, "   avuncular   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][4], checkmarker[4][1],
	    checkmarker[4][2]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_summary, "    f-cousin   %7d      NA(    NA)\n",
	    norelwdata[noped+1][5]);
  else
    fprintf(fp_summary, "    f-cousin   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][5], checkmarker[5][1],
	    checkmarker[5][2]);

  if (norelwdata[noped+1][6] == 0)
    fprintf(fp_summary, "   unrelated   %7d      NA(    NA)\n",
	    norelwdata[noped+1][6]);
  else
    fprintf(fp_summary, "   unrelated   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][6], checkmarker[6][1],
	    checkmarker[6][2]);

  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_summary, " h-avuncular   %7d      NA(    NA)\n",
	    norelwdata[noped+1][7]);
  else
    fprintf(fp_summary, " h-avuncular   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][7], checkmarker[8][1],
	    checkmarker[7][2]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_summary, "  h-f-cousin   %7d      NA(    NA)\n",
	    norelwdata[noped+1][8]);
  else
    fprintf(fp_summary, "  h-f-cousin   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][8], checkmarker[8][1],
	    checkmarker[8][2]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_summary, "h-sib+f-cous   %7d      NA(    NA)\n",
	    norelwdata[noped+1][9]);
  else
    fprintf(fp_summary, "h-sib+f-cous   %7d  %6.1f(%6.1f)\n",
	    norelwdata[noped+1][9], checkmarker[9][1],
	    checkmarker[9][2]);

  if (norelwdata[noped+1][10] == 0)
    fprintf(fp_summary, "pa-offspring   %7d      NA(    NA)\n",
	    norelwdata[noped+1][10]);
  else
    fprintf(fp_summary, "pa-offspring   %7d  %6.1f(%6.1f)\n\n",
	    norelwdata[noped+1][10], checkmarker[10][1],
	    checkmarker[10][2]);

  fprintf(fp_summary, "# Sample summary of EIBD \n\n");
    fprintf(fp_summary, "                   obs     mu    std(     std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_summary, "       full-sib     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][1] == 1)
    fprintf(fp_summary, "       full-sib %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][1][1], checkallefreq[1][1][2],
	    checkallefreq[1][1][4]);
  else
    fprintf(fp_summary, "       full-sib %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][1][1],
	    checkallefreq[1][1][2],
	    checkallefreq[1][1][4], checkallefreq[1][1][5]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_summary, "       half-sib     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][2] == 1)
    fprintf(fp_summary, "       half-sib %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][2][1], checkallefreq[1][2][2],
	    checkallefreq[1][2][4]);
  else
    fprintf(fp_summary, "       half-sib %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][2][1],
	    checkallefreq[1][2][2],
	    checkallefreq[1][2][4], checkallefreq[1][2][5]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_summary, "       grandp-c     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][3] == 1)
    fprintf(fp_summary, "       grandp-c %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][3][1],
	    checkallefreq[1][3][2],checkallefreq[1][3][4]);
  else
    fprintf(fp_summary, "       grandp-c %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][3][1],
	    checkallefreq[1][3][2],
	    checkallefreq[1][3][4], checkallefreq[1][3][5]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_summary, "      avuncular     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][4] == 1)
    fprintf(fp_summary, "      avuncular %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][4][1],
	    checkallefreq[1][4][2], checkallefreq[1][4][4]);
  else
    fprintf(fp_summary, "      avuncular %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][4][1],
	    checkallefreq[1][4][2],
	    checkallefreq[1][4][4], checkallefreq[1][4][5]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_summary, "       f-cousin     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][5] == 1)
    fprintf(fp_summary, "       f-cousin %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][5][1],
	    checkallefreq[1][5][2],checkallefreq[1][5][4]);
  else
    fprintf(fp_summary, "       f-cousin %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][5][1],
	    checkallefreq[1][5][2],
	    checkallefreq[1][5][4], checkallefreq[1][5][5]);

  fprintf(fp_summary, "      unrelated     NA     NA     NA(      NA)\n");

  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_summary, "    h-avuncular     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][7] == 1)
    fprintf(fp_summary, "   h-avuncular %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][7][1],
	    checkallefreq[1][7][2],checkallefreq[1][7][4]);
  else
    fprintf(fp_summary, "   h-avuncular  %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][7][1],
	    checkallefreq[1][7][2],
	    checkallefreq[1][7][4], checkallefreq[1][7][5]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_summary, "     h-f-cousin     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][8] == 1)
    fprintf(fp_summary, "     h-f-cousin %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][8][1],
	    checkallefreq[1][8][2],checkallefreq[1][8][4]);
  else
    fprintf(fp_summary, "     h-f-cousin %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][8][1],
	    checkallefreq[1][8][2],
	    checkallefreq[1][8][4], checkallefreq[1][8][5]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_summary, " h-sib+f-cousin     NA     NA     NA(      NA)\n");
  else if (norelwdata[noped+1][9] == 1)
    fprintf(fp_summary, " h-sib+f-cousin %6.4f %6.4f %6.4f(      NA)\n",
	    checkallefreq[1][9][1],
	    checkallefreq[1][9][2],checkallefreq[1][9][4]);
  else
    fprintf(fp_summary, " h-sib+f-cousin %6.4f %6.4f %6.4f(%8.6f)\n",
	    checkallefreq[1][9][1],
	    checkallefreq[1][9][2],
	    checkallefreq[1][9][4], checkallefreq[1][9][5]);

  fprintf(fp_summary, "   pa-offspring     NA     NA     NA(      NA)\n\n");


  fprintf(fp_summary, "# Sample summary of AIBS \n\n");
  fprintf(fp_summary, "                   obs     mu(     std)    std(     std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_summary, "       full-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][1] == 1)
    fprintf(fp_summary, "       full-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][1][1], checkallefreq[2][1][2],
	    checkallefreq[2][1][4]);
  else
    fprintf(fp_summary, "       full-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][1][1],
	    checkallefreq[2][1][2], checkallefreq[2][1][3],
	    checkallefreq[2][1][4], checkallefreq[2][1][5]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_summary, "       half-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][2] == 1)
    fprintf(fp_summary, "       half-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][2][1],
	    checkallefreq[2][2][2],checkallefreq[2][2][4]);
  else
    fprintf(fp_summary, "       half-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][2][1],
	    checkallefreq[2][2][2], checkallefreq[2][2][3],
	    checkallefreq[2][2][4], checkallefreq[2][2][5]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_summary, "       grandp-c     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][3] == 1)
    fprintf(fp_summary, "       grandp-c %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][3][1],
	    checkallefreq[2][3][2],checkallefreq[2][3][4]);
  else
    fprintf(fp_summary, "       grandp-c %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][3][1],
	    checkallefreq[2][3][2], checkallefreq[2][3][3],
	    checkallefreq[2][3][4], checkallefreq[2][3][5]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_summary, "      avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][4] == 1)
    fprintf(fp_summary, "      avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][4][1],
	    checkallefreq[2][4][2],checkallefreq[2][4][4]);
  else
    fprintf(fp_summary, "      avuncular %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][4][1],
	    checkallefreq[2][4][2], checkallefreq[2][4][3],
	    checkallefreq[2][4][4], checkallefreq[2][4][5]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_summary, "       f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][5] == 1)
    fprintf(fp_summary, "       f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][5][1],
	    checkallefreq[2][5][2],checkallefreq[2][5][4]);
  else
    fprintf(fp_summary, "       f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][5][1],
	    checkallefreq[2][5][2], checkallefreq[2][5][3],
	    checkallefreq[2][5][4], checkallefreq[2][5][5]);

  fprintf(fp_summary, "      unrelated     NA     NA(      NA)     NA(      NA)\n");

  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_summary, "    h-avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][7] == 1)
    fprintf(fp_summary, "    h-avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][7][1],
	    checkallefreq[2][7][2],checkallefreq[2][7][4]);
  else
    fprintf(fp_summary, "    h-avuncular %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][7][1],
	    checkallefreq[2][7][2], checkallefreq[2][7][3],
	    checkallefreq[2][7][4], checkallefreq[2][7][5]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_summary, "     h-f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][8] == 1)
    fprintf(fp_summary, "     h-f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][8][1],
	    checkallefreq[2][8][2],checkallefreq[2][8][4]);
  else
    fprintf(fp_summary, "     h-f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][8][1],
	    checkallefreq[2][8][2], checkallefreq[2][8][3],
	    checkallefreq[2][8][4], checkallefreq[2][8][5]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_summary, " h-sib+f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][9] == 1)
    fprintf(fp_summary, " h-sib+f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][9][1],
	    checkallefreq[2][9][2],checkallefreq[2][9][4]);
  else
    fprintf(fp_summary, " h-sib+f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[2][9][1],
	    checkallefreq[2][9][2], checkallefreq[2][9][3],
	    checkallefreq[2][9][4], checkallefreq[2][9][5]);

  if (norelwdata[noped+1][10] == 0)
    fprintf(fp_summary, "   pa-offsprin      NA     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][10] == 1)
    fprintf(fp_summary, "   pa-offspring %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[2][10][1],
	    checkallefreq[2][10][2],checkallefreq[2][10][4]);
  else
    fprintf(fp_summary, "   pa-offspring %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n\n",
	    checkallefreq[2][10][1],
	    checkallefreq[2][10][2], checkallefreq[2][10][3],
	    checkallefreq[2][10][4], checkallefreq[2][10][5]);


  fprintf(fp_summary, "# Sample summary of IBS \n\n");
  fprintf(fp_summary, "                   obs     mu(     std)    std(     std)\n");

  if (norelwdata[noped+1][1] == 0)
    fprintf(fp_summary, "       full-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][1] == 1)
    fprintf(fp_summary, "       full-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][1][1],
	    checkallefreq[3][1][2],checkallefreq[3][1][4]);
  else
    fprintf(fp_summary, "       full-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][1][1],
	    checkallefreq[3][1][2], checkallefreq[3][1][3],
	    checkallefreq[3][1][4], checkallefreq[3][1][5]);

  if (norelwdata[noped+1][2] == 0)
    fprintf(fp_summary, "       half-sib     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][2] == 1)
    fprintf(fp_summary, "       half-sib %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][2][1],
	    checkallefreq[3][2][2],checkallefreq[3][2][4]);
  else
    fprintf(fp_summary, "       half-sib %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][2][1],
	    checkallefreq[3][2][2], checkallefreq[3][2][3],
	    checkallefreq[3][2][4], checkallefreq[3][2][5]);

  if (norelwdata[noped+1][3] == 0)
    fprintf(fp_summary, "       grandp-c     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][3] == 1)
    fprintf(fp_summary, "       grandp-c %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][3][1],
	    checkallefreq[3][3][2],checkallefreq[3][3][4]);
  else
    fprintf(fp_summary, "       grandp-c %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][3][1],
	    checkallefreq[3][3][2], checkallefreq[3][3][3],
	    checkallefreq[3][3][4], checkallefreq[3][3][5]);

  if (norelwdata[noped+1][4] == 0)
    fprintf(fp_summary, "      avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][4] == 1)
    fprintf(fp_summary, "      avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][4][1],
	    checkallefreq[3][4][2], checkallefreq[3][4][4]);
  else
    fprintf(fp_summary, "      avuncular %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][4][1],
	    checkallefreq[3][4][2], checkallefreq[3][4][3],
	    checkallefreq[3][4][4], checkallefreq[3][4][5]);

  if (norelwdata[noped+1][5] == 0)
    fprintf(fp_summary, "       f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][5] == 1)
    fprintf(fp_summary, "      f -cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][5][1],
	    checkallefreq[3][5][2],checkallefreq[3][5][4]);
  else
    fprintf(fp_summary, "       f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][5][1],
	    checkallefreq[3][5][2], checkallefreq[3][5][3],
	    checkallefreq[3][5][4], checkallefreq[3][5][5]);

  if (norelwdata[noped+1][6] == 0)
    fprintf(fp_summary, "      unrelated     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][6] == 1)
    fprintf(fp_summary, "      unrelated %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][6][1],
	    checkallefreq[3][6][2],checkallefreq[3][6][4]);
  else
    fprintf(fp_summary, "      unrelated %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][6][1],
	    checkallefreq[3][6][2], checkallefreq[3][6][3],
	    checkallefreq[3][6][4], checkallefreq[3][6][5]);

  if (norelwdata[noped+1][7] == 0)
    fprintf(fp_summary, "    h-avuncular     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][7] == 1)
    fprintf(fp_summary, "   h-avuncular %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][7][1],
	    checkallefreq[3][7][2],checkallefreq[3][7][4]);
  else
    fprintf(fp_summary, "   h-avuncular  %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][7][1],
	    checkallefreq[3][7][2], checkallefreq[3][7][3],
	    checkallefreq[3][7][4], checkallefreq[3][7][5]);

  if (norelwdata[noped+1][8] == 0)
    fprintf(fp_summary, "     h-f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][8] == 1)
    fprintf(fp_summary, "     h-f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][8][1],
	    checkallefreq[3][8][2],checkallefreq[3][8][4]);
  else
    fprintf(fp_summary, "     h-f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][8][1],
	    checkallefreq[3][8][2], checkallefreq[3][8][3],
	    checkallefreq[3][8][4], checkallefreq[3][8][5]);

  if (norelwdata[noped+1][9] == 0)
    fprintf(fp_summary, " h-sib+f-cousin     NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][9] == 1)
    fprintf(fp_summary, " h-sib+f-cousin %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][9][1],
	    checkallefreq[3][9][2],checkallefreq[3][9][4]);
  else
    fprintf(fp_summary, " h-sib+f-cousin %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][9][1],
	    checkallefreq[3][9][2], checkallefreq[3][9][3],
	    checkallefreq[3][9][4], checkallefreq[3][9][5]);

  if (norelwdata[noped+1][10] == 0)
    fprintf(fp_summary, "   pa-offsprin      NA     NA(      NA)     NA(      NA)\n");
  else if (norelwdata[noped+1][10] == 1)
    fprintf(fp_summary, "   pa-offspring  %6.4f %6.4f(      NA) %6.4f(      NA)\n",
	    checkallefreq[3][10][1],
	    checkallefreq[3][10][2],checkallefreq[3][10][4]);
  else
    fprintf(fp_summary, "   pa-offspring %6.4f %6.4f(%8.6f) %6.4f(%8.6f)\n",
	    checkallefreq[3][10][1],
	    checkallefreq[3][10][2], checkallefreq[3][10][3],
	    checkallefreq[3][10][4], checkallefreq[3][10][5]);

  /* -------------------------------------------------------------- */
  stringmap_free();


  fclose(fp_summary);
  fclose(fp_results);
  ///fclose(fp_warnings);
  // fclose(fp_error);

  //free(ped);
  free(nomem);
  free(rel);
  free(norel);
  free(norelwdata);

  free(q);
  free(qcum);
  free(theta);
  free(cenm);
  free(noalle);

  free(valid);
  free(ngeno);
  free(genotype);
  free(eeibdibd);
  free(eeibdeibdibd);
  free(eaibsibd);
  free(eaibsaibsibd);
  free(eibsibd);
  free(eibsibsibd);

  free(geno);
  free(simugeno1);
  free(simugeno2);
  free(nopairsvalid);


  free(buf);

  return 0;

} /* end of main */



/* Perform all necessary tests on a pair of individuals */
double pibs[3];


void test_pair(geno_t **geno, int i, int index1, geno_t **genoB, int j, int index2,  int reltype)
{
  int k, n;
  int  commark;

  double dtemp[10], phat[3+1];
  double stat_obs[3+1], stat_simu[3+1];
  double loglh_obs[NRMLRT+1], loglh_simu[NRMLRT+1];
  double mlrt_obs, mlrt_simu;
  int print_warning = 0; /* Specifies if we need to print a line to file out3 */
  int domlrt, countbig[4+1];
  char ctemp1[MAXSTR], ctemp2[MAXSTR];

  if(alternativeRel != 0) reltype = alternativeRel;

  pMLRT = -1;
  /* Test geno[index1], genoB[index2] */


  /* commark: no. of markers typed in both individuals */
  /* valid: which markers are typed in both */
  commark = 0;
  for (k = 1; k <= nomark[0]; k++) {
    if ( geno[index1][2*k-1] != 0 && geno[index1][2*k] != 0 &&
	genoB[index2][2*k-1] != 0 && genoB[index2][2*k] != 0 ) {
      ++commark;
      valid[k] = 1;
    }
    else
      valid[k] = 0;
  }

  // printf("%s %s %s %s %d commX\n", s(nomem[i][1]), s(geno[index1][0]), s(nomem[j][1]), s(genoB[index2][0]), commark);
  
  
  /* testing only for the pair that commark != 0 */
  if (commark == 0) return; 

  domlrt = 0;

  ++norelwdata[i][NR+1];
  ++norelwdata[noped+1][NR+1];

  for (k = 1; k <= NR; k++)
    if ( reltype == k ) {
      ++norelwdata[i][k];
      ++norelwdata[noped+1][k];
    }

  if(debug) printf("testing pair = %s %s  reltype = %d common typed markers = %d\n",
    s(geno[index1][0]),s(genoB[index2][0]), reltype, commark);

    /* statresult: statistic. mu, std, z_obs, p-value of
       EIBD, AIBS, IBS */

    if(MAP == 1)
      get3statandmustd(geno[index1], genoB[index2], valid, commark,
		     reltype, nochrom, nomark, p012[reltype],
		     q, cenm, maptype, eeibdibd, eeibdeibdibd,
		     eaibsibd, eaibsaibsibd, eibsibd, eibsibsibd,
		     statresult, pibs);
    else
      get3stat(geno[index1], genoB[index2], valid, commark,
		     reltype, nochrom, nomark,
		     p012[reltype], q, stat_simu, pibs);

    sprintf(ctemp1, "%%%ds %%%ds %%%ds %%%ds %%7d %%%dd ",
	    max_fid, max_iid, max_fid, max_iid, mmark);

    sprintf(ctemp2, ctemp1,
	    s(nomem[i][1]), s(geno[index1][0]), s(nomem[j][1]), s(genoB[index2][0]),
	    reltype, commark);

    fprintf(fp_results, "%s", ctemp2);

    /* Use EM algorithm to estimate p0, p1, and p2
       (prob. of 0-1-2 IBD sharing) */
    EMforp012(geno[index1], genoB[index2], valid,
	      q, nochrom, nomark, reltype, phat);

    if(phat[1] == -1) {
      printf("Warning: EM for IBD did not converge for pair: %s %s and  %s %s\n",
              s(nomem[i][1]), s(geno[index1][0]), s(nomem[j][1]), s(genoB[index2][0]));
    }

    /* Printing to file prest_errors and increasig countmlrt*/

    /* for the others, use EIBD and AIBS tests, and also
       use p1 > 0.85 to avoid the true is parent-offspring, but
       null is not, especially the null is full-sib */

    /* EIBD, AIBS are not applicable for unrelated = 6
       use IBS for screening */
    /* EIBD is not applicable for parent-offspring = 10 and MZ
       twins = 11, and AIBS and IBS are not powerful for
       distinguish parent-offspring from full-sib, so use estimate
       of p0, p1, p2 to judge*/

    #if 0
    print_warning = 0; /* Specify if we need to print a line to warnings file out3 */

    if (reltype == 6  && phat[2] > 0.75 ) print_warning = 1;
    if (reltype == 10 && phat[2] < 0.9 ) { ++countmlrt; print_warning = 1; }
    if (reltype == 11 && phat[3] < 0.95) { ++countmlrt; print_warning = 1; }

    if( reltype != 6 && reltype != 10 && reltype != 11)
      if( phat[2] > 0.75 ) print_warning = 1;
    #endif

    if (reltype == 6)
      if ((statresult[3][5] < PTHRESH && statresult[3][4] > 0) || phat[2] > 0.75 )
	++countmlrt;

    if( reltype != 6 && reltype != 10 && reltype != 11) {
      if(statresult[1][5] < PTHRESH || statresult[2][5] < PTHRESH || phat[2] > 0.75 )  ++countmlrt;
    }

    if(debug > 1) printf("%s %6d %6d %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n\n ",ctemp2,
	    0, 0, stat_simu[3],
	    pibs[0], pibs[1], pibs[2],
	    phat[1], phat[2], phat[3]);


    #if 0
    if ( print_warning == 1 ) {
      fprintf(fp_warnings, "%s", ctemp2);
      fprintf(fp_warnings,"  %8.4f %8.4f %8.4f  %8.4f %8.4f %8.4f\n",
	  p012[reltype][1],p012[reltype][2],p012[reltype][3],
	  phat[1], phat[2], phat[3]);

      fflush(fp_warnings);
    }

    #endif
    /* END: Printing to file prest_errors */


    /* Print Results to file prest_results */
    if(MAP == 1) {
      if (reltype == 6)
	fprintf(fp_results," %6d %6d %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f  %6s %6s %6.4f ",
		0, 0, statresult[3][1],
		pibs[0], pibs[1], pibs[2],
		phat[1], phat[2], phat[3],
		NA, NA, statresult[3][5]);

      else if (reltype == 10 || reltype == 11)
	fprintf(fp_results," %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f  %6s %6.4f %6.4f ",
		statresult[1][1],statresult[2][1],statresult[3][1],
		pibs[0], pibs[1], pibs[2],
		phat[1], phat[2], phat[3],
		NA, statresult[2][5], statresult[3][5]);

      else
	fprintf(fp_results," %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f  %6.4f %6.4f %6.4f ",
		statresult[1][1],statresult[2][1],statresult[3][1],
		pibs[0], pibs[1], pibs[2],
		phat[1], phat[2], phat[3],
		statresult[1][5], statresult[2][5], statresult[3][5]);

    }

    if (MAP == 0) {
      if (reltype == 6)

	fprintf(fp_results," %6d %6d %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
	    0, 0, stat_simu[3],
	    pibs[0], pibs[1], pibs[2],
	    phat[1], phat[2], phat[3]);

      else

	fprintf(fp_results," %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f ",
	    stat_simu[1], stat_simu[2], stat_simu[3],
	    pibs[0], pibs[1], pibs[2],
	    phat[1], phat[2], phat[3]);

    }


    /* for checking allele freq. */
    if( MAP ==  1)
      for (k = 1; k <= 3; k++) {
	checkallefreq[k][reltype][1] += statresult[k][1];
	checkallefreq[k][reltype][2] += statresult[k][2];
	checkallefreq[k][reltype][3] += statresult[k][2]*statresult[k][2];
	checkallefreq[k][reltype][4] += statresult[k][3];
	checkallefreq[k][reltype][5] += statresult[k][3]*statresult[k][3];
      }

    checkmarker[reltype][1] += commark;
    checkmarker[reltype][2] += commark*commark;


    /* if user want to do MLRT test */
    if ( OPTION == 2 && MAP == 1 ) {

      /* We want to do MLRT for all pairs (as per Lei's suggestion) */
      domlrt = 1;

      /* if the pair is qualified */
      if ( domlrt == 1 ) {

	for (k = 1; k <= 4; k++) countbig[k] = 0;

	++nopairsvalid[i][NR+1];
	++nopairsvalid[noped+1][NR+1];

	for (k = 1; k <= NR; k++) if (reltype == k) {
	    ++nopairsvalid[i][k];
	    ++nopairsvalid[noped+1][k];
	}


	/* observed EIBD, AIBS, IBS statistics */
	for (k = 1; k <= 3; k++)
	  stat_obs[k] = statresult[k][1];

	/* observed log( likelihood ) under each reltype */
	getapprloglikelihood(geno[index1], genoB[index2],
			     valid, q, cenm, nochrom, nomark, maptype,
			     loglh_obs);

	/* MLRT = max(log(likelihood_alt)) -
	   log(likelihood_null) */
	if ( reltype != 1 )
	  dtemp[1] = loglh_obs[1];
	else
	  dtemp[1] = loglh_obs[2];


	for (k = 1; k <= NRMLRT; k++)
	  if ( k != reltype ) {
	    if ( loglh_obs[k] > dtemp[1])
	      dtemp[1] = loglh_obs[k];
	  }



	mlrt_obs = dtemp[1] - loglh_obs[reltype];

	/* APO TESTING */
	/* double *temp = new_double1d(NREP); */

	/* simulate genotype data for the pair and calculate the
	   statistics */

	if(debug) srand(0);

	for (n = 1; n <= NREP; n++) {
	  double phat;

	  simugeno1geno2(q, qcum, nochrom, nomark, noalle, reltype,
			 chromlength, cenm, simugeno1, simugeno2);

	  get3stat(simugeno1, simugeno2, valid,
		   commark, reltype, nochrom, nomark,
		   p012[reltype], q, stat_simu, pibs);

	  /* delete later */
	  /*  checkmu[reltype] += stat_simu[1]; */

	  getapprloglikelihood(simugeno1, simugeno2,
			       valid, q, cenm, nochrom, nomark, maptype,
			       loglh_simu);

	  if (reltype != 1)
	    dtemp[2] = loglh_simu[1];
	  else
	    dtemp[2] = loglh_simu[2];

	  for (k = 1; k <= NRMLRT; k++) {
	    if (k != reltype) {
	      if ( loglh_simu[k] > dtemp[2])
		dtemp[2] = loglh_simu[k];
	    }
	  }

	  mlrt_simu = dtemp[2] - loglh_simu[reltype];

	  /* keep track of no. obs < simu or obs > simu  for the
	     calculation of empirical p-value */
	  for (k = 1; k <= 3; k++) {
	    if (stat_simu[k] > stat_obs[k])
	      ++countbig[k];
	  }

	  if ( mlrt_simu > mlrt_obs )
	    ++countbig[4];

	  /* APO TESTING */
	  /* printf("%8.9f ", 2*(countbig[k]/(n+0.0))); */

	  /* Estimate of p based on the first n simulations */
	  /* if(debug) continue; */

	  phat = (double)countbig[4]/n;
	  if(phat > 0.5) phat = 1-phat;
	  phat = 2*phat;

	  /* Check if phat has enough precision */

	  if(ALPHA > 0 && (n%10)==0) {
	     /* printf("ALPHA=%f n=%d nrep=%f\n", ALPHA,n, (1 / (ALPHA*ALPHA)) * (1-phat)/phat ); */
	  }

	  // Check if we have enough samples for the required precision
	  if(ALPHA > 0 && n > 100) if( (1 / (ALPHA*ALPHA)) * (1-phat)/phat <= n) break;

	  // Old procedure , should be removed
	  if(ALPHA == 0) {
            //Old procedure

	    /*if(n ==     100 && phat >= 0.5) break;
	    if(n ==     400 && phat >= 0.2) break;
	    if(n ==     900 && phat >= 0.1) break;
	    if(n ==    1900 && phat >= 0.05) break;
	    if(n ==    9900 && phat >= 0.01) break;
	    if(n ==   19900 && phat >= 0.005) break;
	    if(n ==   99900 && phat >= 0.001) break;
	    if(n ==  999900 && phat >= 0.0001) break;
	    if(n == 9999900 && phat >= 0.00001) break;*/
	  }

	}/* loop for NREP replicates simulation */


	//if(debug) printf("MLRT %d %d %d %d %d %d\n", n, countbig[1],countbig[2],countbig[3],countbig[4],n);

	/*
	double p = (countbig[4]/(NREP+0.0));
	   estimate of iterations needed
	double estN = (2.6*2.6)*p*(1-p)/(0.1*0.1);

	if(debug) printf("\nfinal=%f est=%f \n",p,estN );
	if(debug) if(estN < NREP && estN >0) {
		printf("diff = %f\n", temp[(int)estN] - p);
	}
	*/

	for (k = 1; k <= 4; k++) {
	    double v;

	    /* EIBD, AIBS are not applicable for unrelated = 6 */

	    if(reltype == 6 && (k==1 || k==2))
	      { fprintf(fp_results, "       NA "); continue; }

	    if((reltype == 10 || reltype == 11) && (k==1))
	      { fprintf(fp_results, "       NA "); continue; }

	    v = countbig[k]/((double)n+0.0);
	    if ( v > 0.5) { v = 1-v; }
	    v = 2*v;

	    fprintf(fp_results, " %8.6f ", v);
            pMLRT = v;

	    if(debug && k==4) printf("MLRT stopped after %d iterations: p=%f\n",n,v );
	}

      }

      /* if the pair is not qualified */
      else {
	for (k = 1; k <= 4; k++)
	  fprintf(fp_results, "       NA ");
      }

    }

    ibs_stats(geno,i,index1,genoB,j,index2,reltype);

    fprintf(fp_results, "\n");
    fflush(fp_results);
}


 void ibs_compute(double result[6], int *m, int n, double th, int *pibs2)
{
  int k;

  double parity1 = 0;

  for (k = 2; k <= n; k++) if(m[k] != m[k-1]) parity1 ++;
  parity1 /= n;

  int max1 = 0;

  for (k = 2; k <= n; k++) {
    if(m[k]) m[k] = m[k-1]+1;
    if(m[k] > max1) max1 = m[k];
  }

  for (k = n; k >= 2; k--) if(m[k-1] && m[k]) m[k-1] = m[k];

  double meanSize = 0;
  for (k = 1; k <= n; k++) meanSize += m[k];
  meanSize = meanSize/n;

  double reg = 0;
  for (k = 1; k <= n; k++) if(m[k] > th && (pibs2 == NULL || pibs2[k] <= th)) reg++;
  //for (k = 1; k <= n; k++) if(m[k] > th) reg++;

  result[1] = max1;
  result[2] = meanSize;
  result[3] = reg/n;
}

void ibs_stats(geno_t **geno, int i, int index1, geno_t **genoB, int j, int index2,  int reltype)
{
  int k;

  int n = 0;

  for (k = 1; k <= nomark[0]; k++) if(valid[k]) {
    int a1 = geno[index1][2*k-1];
    int a2 = geno[index1][2*k];
    int b1 = genoB[index2][2*k-1];
    int b2 = genoB[index2][2*k];

    n ++;

    ibs2[n] = a1==b1 && a2==b2; // IBS=2
    ibs1[n] = a1==b1 || a1==b2 || a2==b1 || a2==b2; // IBS>=1

  }


  //for (k = 2; k <= n; k++) printf("%d",ibs1[k]);printf("\n");
  //for (k = 2; k <= n; k++) printf("%d",ibs2[k]);printf("\n");

  double ibs1par[6], ibs2par[6];
  ibs_compute(ibs2par, ibs2, n, pIBS2_th, NULL);
  ibs_compute(ibs1par, ibs1, n, pIBS0_th, ibs2);

  fprintf(fp_results, " %8d %8.2f", (int)ibs1par[1], ibs1par[2]); // ibs1par[3]);
  fprintf(fp_results, "  %8d %8.2f", (int)ibs2par[1], ibs2par[2]); // ibs2par[3]);
}

/* for find all the relationships considered */

int findrelation(pedigree &ped, int pedsize, int** rel)
{
  int i, j, findit;
  int nofr, inbred;
  int un, ne, np, nup;
  int gp, gc, gcp, gcup;
  int com, ucom;
  int c1p, c2p, c1up, c2up;
  int reltype;

  nofr = 0;
  for (i = 1; i <= (pedsize-1) ; i++) {
    for (j = (i+1); j <= pedsize; j++) {
      /* find relationship between i and j */
      findit = 0;

      if(0) printf("find relation %d,%d,%d,%d,%d <-> %d,%d,%d,%d,%d\n",
              ped[i][1],ped[i][2],ped[i][3],ped[i][4],ped[i][5],
              ped[j][1],ped[j][2],ped[j][3],ped[j][4],ped[j][5]
              );

      /* parent-offspring */
      /* j is the parent of i, and j is not inbred */
      if ( ped[i][4] == j || ped[i][5] == j ) {
	if( inbred_check(ped, pedsize, j) == 0 ) {
	  reltype = 10;
	  findit = 1;
	}
      }
      if (!findit ) {
	/* i is the parent of j, and i is not inbred */
	if ( ped[j][4] == i || ped[j][5] == i ) {
	  if( inbred_check(ped, pedsize, i) == 0 ) {
	    findit = 1;
	    reltype = 10;
	  }
	}
      }
      /* not parent-offspring */
      if ( (!findit) && (! (ped[i][4] == j || ped[i][5] == j ||
			    ped[j][4] == i || ped[j][5] == i) )) {
	/* grandparent-grandchild */
	if ( findit != 1 ) {
	  /* i is the grandparent, through j'th father */
	  gp = i;
	  gc = j;
	  gcp = 4;
	  gcup = 5;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 )
	    reltype = 3;
	}
	if ( findit != 1 ) {
	  /* i is the grandparent, through j'th mother */
	  gp = i;
	  gc = j;
	  gcp = 5;
	  gcup = 4;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 )
	    reltype = 3;
	}
	if ( findit != 1 ) {
	  /* j is the grandparent, through i'th father */
	  gp = j;
	  gc = i;
	  gcp = 4;
	  gcup = 5;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 )
	    reltype = 3;
	}
	if ( findit != 1 ) {
	  /* j is the grandparent, through j'th mother */
	  gp = j;
	  gc = i;
	  gcp = 5;
	  gcup = 4;
	  findit = check_grandpc(gp, gc, gcp, gcup, ped, pedsize);
	  if ( findit == 1 )
	    reltype = 3;
	}
      	/* unrelated */
	if ( findit != 1 ) {
	  inbred = check(ped, pedsize, i, j);
	  if ( inbred == 0 ) {
	    findit = 1;
	    reltype = 6;
	  }
	}
      	/* for full-sib, half-sib, avuncular, first-cousin,
	   half-plus-cousin relative pairs, no founder should be
	   involved */
	if ( findit != 1 &&
	     ped[i][4] != 0 && ped[i][5] != 0 &&
	     ped[j][4] != 0 && ped[j][4] != 0 ) {
	  /* full-sib
	     parents are the same and parents are unrelated,
	     and parents are not inbred themselves*/
	  if ( (ped[i][4] == ped[j][4]) && (ped[i][5] == ped[j][5]) ) {
	    inbred = check(ped, pedsize, ped[i][4], ped[i][5]);
	    inbred = inbred + inbred_check(ped, pedsize, ped[i][4]);
	    inbred = inbred + inbred_check(ped, pedsize, ped[i][5]);
	    if ( inbred == 0 ) {
	      findit = 1;
	      reltype = 1;
	    }
	  }
	  /* avuncular e.g. uncle-nephew */
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th father */
	    un = i;
	    ne = j;
	    np = 4;
	    nup = 5;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 4;
	  }
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th mother */
	    un = i;
	    ne = j;
	    np = 5;
	    nup = 4;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 4;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th father */
	    un = j;
	    ne = i;
	    np = 4;
	    nup = 5;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 4;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th mother */
	    un = j;
	    ne = i;
	    np = 5;
	    nup = 4;
	    findit = check_avuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 4;
	  }
	  /* half-sib */
	  if ( findit != 1 ) {
	    /* same mother */
	    com = 4;
	    ucom = 5;
	    findit = check_halfsib(i, j, com, ucom, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 2;
	  }
	  if (findit != 1) {
	    /* same father */
	    com = 5;
	    ucom = 4;
	    findit = check_halfsib(i, j, com, ucom, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 2;
	  }

	  /* first-cousin */
	  if ( findit != 1 ) {
	    /* i's father is sibling with j's father */
	    c1p = 4;
	    c2p = 4;
	    c1up = 5;
	    c2up = 5;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up,
				  ped, pedsize);
	    if ( findit == 1 )
	      reltype = 5;
	  }
	  if ( findit != 1 ) {
	    /* i's father is sibling with j's mother */
	    c1p = 4;
	    c2p = 5;
	    c1up = 5;
	    c2up = 4;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up,
				  ped, pedsize);
	    if ( findit == 1 )
	      reltype = 5;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is sibling with j's father */
	    c1p = 5;
	    c2p = 4;
	    c1up = 4;
	    c2up = 5;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up,
				  ped, pedsize);
	    if ( findit == 1 )
	      reltype = 5;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is sibling with j's mother */
	    c1p = 5;
	    c2p = 5;
	    c1up = 4;
	    c2up = 4;
	    findit = check_cousin(i, j, c1p, c2p, c1up, c2up,
				  ped, pedsize);
	    if ( findit == 1 )
	      reltype = 5;
	  }

	  /* half-avuncular e.g. half-uncle-nephew */
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th father */
	    un = i;
	    ne = j;
	    np = 4;
	    nup = 5;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 7;
	  }
	  if ( findit != 1 ) {
	    /* if i is the uncle and sibling with j'th mother */
	    un = i;
	    ne = j;
	    np = 5;
	    nup = 4;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 7;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th father */
	    un = j;
	    ne = i;
	    np = 4;
	    nup = 5;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 7;
	  }
	  if ( findit != 1 ) {
	    /* if j is the uncle and sibling with j'th mother */
	    un = j;
	    ne = i;
	    np = 5;
	    nup = 4;
	    findit = check_halfavuncular(un, ne, np, nup, ped, pedsize);
	    if ( findit == 1 )
	      reltype = 7;
	  }

	  /* half-first-cousin */
	  if ( findit != 1 ) {
	    /* i's father is halfsibling with j's father */
	    c1p = 4;
	    c2p = 4;
	    c1up = 5;
	    c2up = 5;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up,
				      ped, pedsize);
	    if ( findit == 1 )
	      reltype = 8;
	  }
	  if ( findit != 1 ) {
	    /* i's father is halfsibling with j's mother */
	    c1p = 4;
	    c2p = 5;
	    c1up = 5;
	    c2up = 4;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up,
				      ped, pedsize);
	    if ( findit == 1 )
	      reltype = 8;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is halfsibling with j's father */
	    c1p = 5;
	    c2p = 4;
	    c1up = 4;
	    c2up = 5;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up,
				      ped, pedsize);
	    if ( findit == 1 )
	      reltype = 8;
	  }
	  if ( findit != 1 ) {
	    /* i's mother is halfsibling with j's mother */
	    c1p = 5;
	    c2p = 5;
	    c1up = 4;
	    c2up = 4;
	    findit = check_halfcousin(i, j, c1p, c2p, c1up, c2up,
				      ped, pedsize);
	    if ( findit == 1 )
	      reltype = 8;
	  }

	  /* half-sib plus first cousin */
	  if (findit != 1) {
	    com = 4;
	    ucom = 5;
	    findit =
	      check_halfpluscousin(i, j, com, ucom, ped, pedsize);
	    if (findit == 1)
	      reltype = 9;
	  }
	  if (findit != 1) {
	    com = 5;
	    ucom = 4;
	    findit =
	      check_halfpluscousin(i, j, com, ucom, ped, pedsize);
	    if (findit == 1)
	      reltype = 9;
	  }
	}
      }

      if ( findit == 0 && all_pairs == 1) {
        printf("Warning: Relationship between %s and %s is not modelled. Setting as unrelated.\n", s(ped[i][1]), s(ped[j][1]) );
        findit = 1;
        reltype = 6;
      }

      if ( findit == 1 ) {
	++nofr;
	rel[nofr][1] = ped[i][1];
	rel[nofr][2] = ped[j][1];
	rel[nofr][3] = reltype;
      }

    }
  }

  return nofr;
}

int check_halfsib(int i, int j, int com, int ucom, pedigree &ped,
		  int pedsize)
{
  int findit, inbred;

  findit = 0;
  if ( ped[i][com] == ped[j][com] ) {
    /* the other parts unrelated, and the common parent is not inbred,
       parents are unrelated */
    inbred = check(ped, pedsize, ped[i][ucom],
				   ped[j][ucom]);
    inbred += inbred_check(ped, pedsize, ped[i][com]);
    inbred += check(ped, pedsize, ped[i][4], ped[i][5]);
    inbred += check(ped, pedsize, ped[j][4], ped[j][5]);
    if ( inbred == 0 )
      findit = 1;
  }

  return findit;
}

int check_grandpc(int gp, int gc, int gcp, int gcup, pedigree &ped,
		  int pedsize)
{
  int findit, inbred;

  findit = 0;
  if ( ped[gc][gcp] != 0 ) {
    if ( gp == ped[ped[gc][gcp]][4] || gp == ped[ped[gc][gcp]][5] ) {
      /* grandparents should not be related, and unrelated with
	 the other parent of the grandchild */
      inbred = check(ped, pedsize,
		       ped[ped[gc][gcp]][4],ped[ped[gc][gcp]][5]);
      inbred += check(ped, pedsize,
		      ped[ped[gc][gcp]][4],ped[gc][gcup]);
      inbred += check(ped, pedsize,
		      ped[ped[gc][gcp]][5],ped[gc][gcup]);
      if (inbred == 0)
	findit = 1;
    }
  }

  return findit;
}

int check_avuncular(int un, int ne, int np, int nup, pedigree &ped,
		    int pedsize)
{
  int findit, inbred;

  findit = 0;
  /* uncles's parents are the same as the parents of nephew's related
     parent */
  if ( (ped[ne][np] != 0 &&
	ped[un][4] == ped[ped[ne][np]][4] &&
	ped[un][5] == ped[ped[ne][np]][5]) ) {
    /* nephew's other parent should not be related with the grandparents,
       and gradparents should not be related, and should not be inbred */
    inbred = check(ped, pedsize, ped[un][4], ped[ne][nup]);
    inbred += check(ped, pedsize, ped[un][5], ped[ne][nup]);
    inbred += check(ped, pedsize, ped[un][4], ped[un][5]);
    inbred += inbred_check(ped, pedsize, ped[un][4]);
    inbred += inbred_check(ped, pedsize, ped[un][5]);
    if (inbred == 0)
      findit = 1;
  }

  return findit;
}

int check_halfavuncular(int un, int ne, int np, int nup, pedigree &ped,
			int pedsize)
{
  int findit, inbred;
  int unrp, unurp;

  findit = 0;
  /* uncles's mother are the same as the mother of nephew's related
     parent */
  unrp = 5;
  unurp = 4;
  if ( ped[ne][np] != 0 &&
       ped[un][unrp] == ped[ped[ne][np]][unrp] ){
    /* uncles's father are unrelatd with the father of of nephew's related
       parent, and they are not inbred and nephew's unrelated
       parent are not related to them */
    inbred = check(ped, pedsize, ped[un][unurp], ped[ped[ne][np]][unurp]);
    inbred += check(ped, pedsize, ped[un][unurp], ped[un][unrp]);
    inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[un][unrp]);
    inbred += check(ped, pedsize, ped[un][unurp], ped[ne][nup] );
    inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[ne][nup] );
    inbred += check(ped, pedsize, ped[un][unrp], ped[ne][nup] );
    inbred += inbred_check(ped, pedsize, ped[un][unrp]);
    inbred += inbred_check(ped, pedsize, ped[un][unurp]);
    inbred += inbred_check(ped, pedsize, ped[ped[ne][np]][unurp]);
    if (inbred == 0)
      findit = 1;
  }
  if(findit == 0) {
    /* uncles's father are the same as the father of nephew's related
       parent */
    unrp = 4;
    unurp = 5;
    if ( ped[ne][np] != 0 &&
	 ped[un][unrp] == ped[ped[ne][np]][unrp] ){
      /* uncles's mother are unrelatd with the mother of of nephew's related
	 parent, and they are not inbred and nephew's unrelated
	 parent are not related to them */
      inbred = check(ped, pedsize, ped[un][unurp], ped[ped[ne][np]][unurp]);
      inbred += check(ped, pedsize, ped[un][unurp], ped[un][unrp]);
      inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[un][unrp]);
      inbred += check(ped, pedsize, ped[un][unurp], ped[ne][nup] );
      inbred += check(ped, pedsize, ped[ped[ne][np]][unurp], ped[ne][nup] );
      inbred += check(ped, pedsize, ped[un][unrp], ped[ne][nup] );
      inbred += inbred_check(ped, pedsize, ped[un][unrp]);
      inbred += inbred_check(ped, pedsize, ped[un][unurp]);
      inbred += inbred_check(ped, pedsize, ped[ped[ne][np]][unurp]);
      if (inbred == 0)
	findit = 1;
    }
  }
  return findit;
}



int check_cousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		 pedigree &ped, int pedsize)
{
  int findit, inbred;

  findit = 0;
  if ( ped[ped[i][c1p]][4] == ped[ped[j][c2p]][4] &&
       ped[ped[i][c1p]][4] != 0 &&
       ped[ped[i][c1p]][5] == ped[ped[j][c2p]][5] &&
       ped[ped[i][c1p]][5] != 0 ) {
    /* grandparents are unrelated, non-inbred, the other parents are
       not related with the grandparents */
    inbred = check(ped, pedsize,
		   ped[ped[i][c1p]][4], ped[ped[i][c1p]][5]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][4]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][5]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[ped[i][c1p]][4]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[ped[i][c1p]][5]);
    inbred += check(ped, pedsize,
		    ped[j][c2up], ped[ped[i][c1p]][4]);
    inbred += check(ped, pedsize,
		    ped[j][c2up], ped[ped[i][c1p]][5]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[j][c2up]);
    if (inbred == 0)
      findit = 1;
  }

  return findit;
}

int check_halfcousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		     pedigree &ped, int pedsize)
{
  int findit, inbred;
  int rp, urp;

  findit = 0;
  /* related grandparent is the grandmother*/
  rp = 5;
  urp = 4;
  if ( ped[ped[i][c1p]][rp] == ped[ped[j][c2p]][rp] &&
       ped[ped[i][c1p]][rp] != 0){
    /* the other grandparents are unrelated, non-inbred, and the other
       parents are not related with each other and not related with
       the three grandparents */
    inbred = check(ped, pedsize,
		   ped[ped[i][c1p]][urp], ped[ped[j][c2p]][urp]);
    inbred = check(ped, pedsize,
		   ped[ped[i][c1p]][urp],ped[ped[i][c1p]][rp]);
    inbred = check(ped, pedsize,
		   ped[ped[i][c1p]][rp], ped[ped[j][c2p]][urp]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][rp]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][urp]);
    inbred += inbred_check(ped, pedsize, ped[ped[j][c2p]][urp]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[j][c2up]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[ped[i][c1p]][rp]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[ped[i][c1p]][urp]);
    inbred += check(ped, pedsize,
		    ped[i][c1up], ped[ped[j][c2p]][urp]);
    inbred += check(ped, pedsize,
		    ped[j][c2up], ped[ped[i][c1p]][rp]);
    inbred += check(ped, pedsize,
		    ped[j][c2up], ped[ped[i][c1p]][urp]);
    inbred += check(ped, pedsize,
		    ped[j][c2up], ped[ped[j][c2p]][urp]);
    if (inbred == 0)
      findit = 1;
  }
  /* related grandparent is the grandfather*/
  if(findit == 0) {
    rp = 4;
    urp = 5;
    if ( ped[ped[i][c1p]][rp] == ped[ped[j][c2p]][rp] &&
	 ped[ped[i][c1p]][rp] != 0){
      /* the other grandparents are unrelated, non-inbred, and the other
	 parents are not related with each other and not related with
	 the three grandparents */
      inbred = check(ped, pedsize,
		     ped[ped[i][c1p]][urp], ped[ped[j][c2p]][urp]);
      inbred = check(ped, pedsize,
		     ped[ped[i][c1p]][urp],ped[ped[i][c1p]][rp]);
      inbred = check(ped, pedsize,
		     ped[ped[i][c1p]][rp], ped[ped[j][c2p]][urp]);
      inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][rp]);
      inbred += inbred_check(ped, pedsize, ped[ped[i][c1p]][urp]);
      inbred += inbred_check(ped, pedsize, ped[ped[j][c2p]][urp]);
      inbred += check(ped, pedsize,
		      ped[i][c1up], ped[j][c2up]);
      inbred += check(ped, pedsize,
		      ped[i][c1up], ped[ped[i][c1p]][rp]);
      inbred += check(ped, pedsize,
		      ped[i][c1up], ped[ped[i][c1p]][urp]);
      inbred += check(ped, pedsize,
		      ped[i][c1up], ped[ped[j][c2p]][urp]);
      inbred += check(ped, pedsize,
		      ped[j][c2up], ped[ped[i][c1p]][rp]);
      inbred += check(ped, pedsize,
		      ped[j][c2up], ped[ped[i][c1p]][urp]);
      inbred += check(ped, pedsize,
		      ped[j][c2up], ped[ped[j][c2p]][urp]);
      if (inbred == 0)
	findit = 1;
    }
  }

  return findit;
}

int check_halfpluscousin(int i, int j, int com, int ucom, pedigree &ped,
			 int pedsize)
{
  int findit, inbred;

  findit = 0;
  if ( ped[i][com] == ped[j][com] &&
       ped[ped[i][ucom]][4] == ped[ped[j][ucom]][4] &&
       ped[ped[i][ucom]][4] != 0 &&
       ped[ped[i][ucom]][5] == ped[ped[j][ucom]][5] &&
       ped[ped[i][ucom]][5] != 0 ) {
    /* grandparents are unrelated, non-inbred, and unrelated with the
       common parent, and the common parent should not be inbred */
    inbred = check(ped, pedsize,
		   ped[ped[i][ucom]][4], ped[ped[i][ucom]][5]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][ucom]][4]);
    inbred += inbred_check(ped, pedsize, ped[ped[i][ucom]][5]);
    inbred += check(ped, pedsize,
		  ped[i][com], ped[ped[i][ucom]][4]);
    inbred += check(ped, pedsize,
		  ped[i][com], ped[ped[i][ucom]][5]);
    inbred += inbred_check(ped, pedsize, ped[i][com]);
    if (inbred == 0)
      findit = 1;
  }

  return findit;
}


void build_anclist(pedigree &ped, int pedsize, int id1, int
		   anclist[MANCE][3], int* totalnum, int flag)
{
  static int numinlist;

  if ( flag == 0 ) {
    numinlist = 1;
    anclist[numinlist][1] = id1;
    anclist[numinlist][2] = 0;
  }
  if ( ped[id1][4] == 0 && ped[id1][5] == 0 ) {
    *totalnum = numinlist;
    return;
  }

  ++numinlist;
  anclist[numinlist][1] = ped[id1][4];
  anclist[numinlist][2] = ped[id1][5];

  build_anclist(ped, pedsize, ped[id1][4], anclist, totalnum, 1);
  build_anclist(ped, pedsize, ped[id1][5], anclist, totalnum, 1);

  *totalnum = numinlist;

  return;

}

int check(pedigree &ped, int pedsize, int id1, int id2)
{
  static int anclist[MANCE][3];
  int totalnum;

  if ( id1 == id2 && id1 != 0 )
    return 1;
  if ( check1_anclist(ped, pedsize, id1) ||
       check1_anclist(ped, pedsize, id2) )
    return 1;

  build_anclist(ped, pedsize, id1, anclist, &totalnum, 0);

  return check_anclist(ped, pedsize, id2, anclist, totalnum);

}

int check_anclist(pedigree &ped, int pedsize, int id2,
		  int anclist[MANCE][3], int totalnum)
{
  int i, c1, c2;

  if ( ped[id2][4] == 0 || ped[id2][5] == 0 ) {
    for (i = 1; i <= totalnum; ++i) {
      if ( anclist[i][1] == id2 || anclist[i][2] == id2 )
	return 1;
    }
    return 0;
  }

  for (i = 1; i <= totalnum; ++i) {
    if ( anclist[i][1] == ped[id2][4] ||
	 anclist[i][2] == ped[id2][5] )
      return 1;
  }

  c1 = check_anclist(ped, pedsize, ped[id2][4], anclist, totalnum);
  if ( c1 == 1 )
    return 1;
  c2 = check_anclist(ped, pedsize, ped[id2][5], anclist, totalnum);
  if ( c2 == 1 )
    return 1;

  return 0;
}

int check1_anclist(pedigree &ped, int pedsize, int id1)
{

  int c1, c2;

  if ( ped[id1][4] == 0 && ped[id1][5] == 0 )
    return 0;
  if ( ped[id1][4] == 0 && ped[id1][5] != 0 )
    return 1;
  if ( ped[id1][5] == 0 && ped[id1][4] != 0 )
    return 1;

  c1 = check1_anclist(ped, pedsize, ped[id1][4]);
  if ( c1 == 1 )
    return 1;
  c2 = check1_anclist(ped, pedsize, ped[id1][5]);
  if ( c2 == 1 )
    return 1;

  return 0;

}

int inbred_check(pedigree &ped, int pedsize, int id)
{

  if ( ped[id][4]==0 && ped[id][5]==0 )
    return 0;
  if ( ped[id][4]==0 && ped[id][5]!=0 )
    return 1;
  if ( ped[id][4]!=0 && ped[id][5]==0 )
    return 1;

  return check(ped,pedsize,ped[id][4],ped[id][5]);

}


int recoding(pedigree &ped, int* nomem)
{
  int j, k, kkk, done, error;

  //printf("Recoding...\n"); fflush(stdout);
  
  error = 0;
  for (j = 1; j <= nomem[2]; j++) {

    /* bug checking */
    done = 0;
    if ( ped[j][2] == ped[j][3] && ped[j][2] != 0 ) {
      printf("Error in pedigree %s: ", s(nomem[1]));
      printf("%s'th parents have the same id\n",s(ped[j][1]));
      ++error;
    }
    for (k = j+1; k <= nomem[2]; k++) {
      if ( ped[j][1] == ped[k][1] ) {
        printf("Error in pedigree %s: ", s(nomem[1]));
	printf("individual id = %s being used\n",s(ped[j][1]));
	++error;
      }
      if ( ped[j][2] == ped[k][3] && ped[j][2] != 0 ) {
        printf("Error in pedigree %s: ", s(nomem[1]));
	printf("%s has two genders\n",s(ped[j][2]));
      	++error;
      }
      if ( ped[j][3] == ped[k][2] && ped[j][3] != 0 ) {
        printf("Error in pedigree %s: ", s(nomem[1]));
	printf("%s has two genders\n",s(ped[j][3]));
      	++error;
      }
    }


    ped[j][4] = 0;
    ped[j][5] = 0;
    
    /* recoding for founders */
    if ( ped[j][2] == 0 ) ++done;
    if ( ped[j][3] == 0 ) ++done;

    /* recoding for non-founders */
    kkk = 1;
    while ( done != 2 && kkk <= nomem[2] ) {
      if ( ped[j][2] == ped[kkk][1] ) {
	ped[j][4] = kkk;
	++done;
      }
      if ( ped[j][3] == ped[kkk][1] ) {
	ped[j][5] = kkk;
	++done;
      }
      ++kkk;
    }


    /* bug checking */
    if ( done != 2 ) {
      printf( "Warning: pedigree %s is not complete: ", s(nomem[1]));
      printf( "%s's parents: %s and/or %s are missing from pedigree.\n",
              s(ped[j][1]), s(ped[j][2]), s(ped[j][3]));

    //  if(ped[j][4]==0) ped[j][2]=0;
     // if(ped[j][5]==0) ped[j][3]=0;

//      ++error;
    }

   if(done != 2) {
       if(ped[j][2] != 0 && ped[j][4] == 0) {
         int a[6] = {0, ped[j][2], 0,0,0,0 };

         pedmemb t; t = (a);
         ped.push_back(t);
         ped[j][4] = ped.size()-1;
         nomem[2]++;
         
         t = ped[ ped.size()-1 ];
         cout << "Added missing parent: " << ped.size()-1 << ": " << s(t[1]) << " " << s(t[2]) << " " << s(t[3]) << " " << t[4] << " " << t[5] << "\n";
       }

        if(ped[j][3] != 0 && ped[j][5] == 0) {
         int a[6] = {0, ped[j][3], 0,0,0,0 };

         pedmemb t; t = (a);
         ped.push_back(t);
         ped[j][5] = ped.size()-1;
         nomem[2]++;

         t = ped[ ped.size()-1 ];
         cout << "Added missing parent: " << ped.size()-1 << ": " << s(t[1]) << " " << s(t[2]) << " " << s(t[3]) << " " << t[4] << " " << t[5] << "\n";
 
       }
    }

    if(0) printf("%d %d %d %d %d %d\n",ped[j][0],ped[j][1],ped[j][2],ped[j][3],
            ped[j][4], ped[j][5]);


  }

  return error;
}




int checkMend(FILE *out, pedigree &ped, geno_t** geno, int* nomem, int* nomark)
{
  int i, k, m, kkk;
  int error, findit, ok;
  int s1, s2, f1, f2, m1, m2;

  error = 0;
  for (i = 1; i <= nomem[2]; i++) {

    if ( ped[i][4] != 0 && ped[i][5] != 0) {
      kkk = 0;
      
      for (k = 1; k <= nochrom; k++) {
	for (m = (kkk+1); m <= (kkk+nomark[k]); m++) {
	  findit = 0;

	  s1 = geno[i][2*m-1];
	  s2 = geno[i][2*m];
	  f1 = geno[ped[i][4]][2*m-1];
	  f2 = geno[ped[i][4]][2*m];
	  m1 = geno[ped[i][5]][2*m-1];
	  m2 = geno[ped[i][5]][2*m];

          ok = 1;

          if ( s1 == 0 || s2 == 0 ) continue;

          if ( f1 != 0 && f2 != 0 && m1 != 0 && m2 != 0 ) {
            if ( (s1 == f1 || s1 == f2) && (s2 == m1 || s2 == m2) )
              findit = 1;
            if ( !findit ) {
              if ( (s2 == f1 || s2 == f2) && (s1 == m1 || s1 == m2) )
                findit = 1;
              else { ok = 0; }
            }
          }
          else if ( f1 != 0 && f2 != 0 && m1 == 0 && m2 == 0 ) {
            if ( s1 == f1 || s1 == f2 || s2 == f1 || s2 == f2 )
              findit = 1;
            if ( !findit ) { ok = 0; }
          }
          else if ( f1 == 0 && f2 == 0 && m1 != 0 && m2 != 0 ) {
            if ( s1 == m1 || s1 == m2 || s2 == m1 || s2 == m2 )
              findit = 1;
            if ( !findit ) { ok = 0; }
          }


          if(!ok) {
            fprintf(out,"%s\t%s\t%s\t%d,%d\t%d,%d\t%d,%d\n",
                          s(nomem[1]),s(ped[i][1]), names[k][m-kkk],
                          s1, s2, f1, f2, m1, m2);
            ++error;
          }

	}
      	kkk += nomark[k];
      }
    }
  }

  return error;
}


/* for trans. prob., cond. prob. and other prob.*/

/* cp2( ) is the conditional probability of P(X|2 IBD )
   where X is the genotype data for the pair,
   similar for cp1( ) and cp0( ), unordered genotype for a pair */
double cp2(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  if ( f1 == m1 && m1 == f2 && f2 == m2  )
    temp = q[f1]*q[f1];
  else if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) )
    temp = 2*q[f1]*q[m1];
  else
    temp = 0;

  return temp;
}

double cp1(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  if ( f1 == m1 && m1 == f2 && f2 == m2 )
    temp = q[f1]*q[f1]*q[f1];
  else if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) )
    temp = q[f1]*q[m1]*(q[f1] + q[m1]);
  else if ( (f1 == m1 && m1 == f2) || (f1 == m1 && m1 == m2) )
    temp = 2*q[f1]*q[f2]*q[m2];
  else if ( (f2 == m2 && m2 == f1) || (f2 == m2 && m2 == m1) )
    temp = 2*q[f1]*q[m1]*q[f2];
  else if ( f1 == f2 || f1 == m2 )
    temp = 2*q[m1]*q[f2]*q[m2];
  else if ( m1 == f2 || m1 == m2 )
    temp = 2*q[f1]*q[f2]*q[m2];
  else
    temp = 0;

  return temp;
}

double cp0(int f1, int m1, int f2, int m2, double* q)
{
  double temp;

  temp = q[f1]*q[m1]*q[f2]*q[m2];
  if ( f1 == m1 && m1==f2 && f2 == m2 )
    temp = temp;
  else if ( f1 == m1 &&  f2 == m2 )
    temp = 2*temp;
  else if ( f1 == m1 || f2 == m2 || (f1 == f2 && m1 == m2)
	  || (f1 == m2 && m1 == f2) )
    temp = 4*temp;
  else
    temp = 8*temp;

  return temp;
}

void transprob(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1],
	       int reltype)
{
  int i, j, k;
  double ta, phi, psi;

  ta = mapcmtotheta(dij, maptype);
  phi = ta*ta+(1-ta)*(1-ta);
  psi = 2*ta*(1-ta);

  for (i = 1; i <= NRMLRT; i++)
    for (j = 1; j <= 3; j++)
      for (k = 1; k <= 3; k++)
	pibdiibdj[i][j][k] = 0;

  /* in 0-1-2 IBD to 0-1-2 IBD order */

  /* fullsib */
  if (reltype == 1) {
    pibdiibdj[1][1][1] = pow((ta*ta+(1-ta)*(1-ta)),2);
    pibdiibdj[1][1][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][1][3] = pow((2*ta*(1-ta)),2);
    pibdiibdj[1][2][1] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][2][2] = pow((ta*ta+(1-ta)*(1-ta)),2) + pow((2*ta*(1-ta)),2);
    pibdiibdj[1][2][3] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][3][1] = pow((2*ta*(1-ta)),2);
    pibdiibdj[1][3][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
    pibdiibdj[1][3][3] = pow((ta*ta+(1-ta)*(1-ta)),2);
  }
  /* half-sib */
  else if (reltype == 2) {
    pibdiibdj[2][1][1] = ta*ta+(1-ta)*(1-ta);
    pibdiibdj[2][1][2] = 2*ta*(1-ta);
    pibdiibdj[2][2][1] = 2*ta*(1-ta);
    pibdiibdj[2][2][2] = ta*ta+(1-ta)*(1-ta);
  }
  /* grandparent-grandchild */
  else if (reltype == 3) {
    pibdiibdj[3][1][1] = (1-ta);
    pibdiibdj[3][1][2] = ta;
    pibdiibdj[3][2][1] = ta;
    pibdiibdj[3][2][2] = (1-ta);
  }
  /* avuncular */
  else if (reltype ==4) {
    pibdiibdj[4][2][2] = 1-2.5*ta+4*ta*ta-2*ta*ta*ta;
    pibdiibdj[4][2][1] = 1-pibdiibdj[4][2][2];
    pibdiibdj[4][1][2] = 1-pibdiibdj[4][2][2];
    pibdiibdj[4][1][1] = pibdiibdj[4][2][2];
  }
  /* first-cousin */
  else if (reltype == 5) {
    pibdiibdj[5][2][2] = 1-4*ta+7.5*ta*ta- 6*ta*ta*ta+2*ta*ta*ta*ta;
    pibdiibdj[5][2][1] = 1-pibdiibdj[5][2][2];
    pibdiibdj[5][1][2] = (1-pibdiibdj[5][2][2])/3.0;
    pibdiibdj[5][1][1] = 1-pibdiibdj[5][1][2];
  }
  /* unrelated */
  else if (reltype == 6 || reltype == 12)
    pibdiibdj[6][1][1] = 1;

  /* half-avuncular */
  else if (reltype == 7) {
    pibdiibdj[7][2][2]=(1-ta)*(1-ta)*
      (1-ta)+ta*ta*(1-ta);
    pibdiibdj[7][2][1]=1-pibdiibdj[7][2][2];
    pibdiibdj[7][1][2]=(1-pibdiibdj[7][2][2])/3.0;
    pibdiibdj[7][1][1]=1-pibdiibdj[7][1][2];
  }
  /* half-cousin */
  else if (reltype == 8) {
    pibdiibdj[8][2][2]=(1-ta)*(1-ta)*(1-ta)*(1-ta)+ta*ta*(1-ta)*(1-ta);
    pibdiibdj[8][2][1]=1-pibdiibdj[8][2][2];
    pibdiibdj[8][1][2]=(1-pibdiibdj[8][2][2])/7.0;
    pibdiibdj[8][1][1]=1-pibdiibdj[8][1][2];
  }
  /* half-sib plus first-cousin */
  else if (reltype == 9) {
    pibdiibdj[9][3][3] = pibdiibdj[5][2][2]*phi;
    pibdiibdj[9][3][1] = pibdiibdj[5][2][1]*psi;
    pibdiibdj[9][3][2] = 1-pibdiibdj[9][3][3]-pibdiibdj[9][3][1];

    pibdiibdj[9][1][3] = pibdiibdj[5][1][2]*psi;
    pibdiibdj[9][1][1] = pibdiibdj[5][1][1]*phi;
    pibdiibdj[9][1][2] = 1-pibdiibdj[9][1][3]-pibdiibdj[9][1][1];

    pibdiibdj[9][2][3] = 1.0/3*pibdiibdj[5][2][2]*psi+
      2.0/3*pibdiibdj[5][1][2]*phi;
    pibdiibdj[9][2][1] = 1.0/3*pibdiibdj[5][2][1]*phi+
      2.0/3*pibdiibdj[5][1][1]*psi;
    pibdiibdj[9][2][2] = 1-pibdiibdj[9][2][3]-pibdiibdj[9][2][1];
  }
  /* parent-offspring */
  else if (reltype == 10)
    pibdiibdj[10][2][2] = 1;
  /* MZ twins */
  else if (reltype == 11)
    pibdiibdj[11][3][3] = 1;

}

void transproball(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1])
{
  int i, j, k;
  double ta, phi, psi;

  ta = mapcmtotheta(dij, maptype);
  phi = ta*ta+(1-ta)*(1-ta);
  psi = 2*ta*(1-ta);

  for (i = 1; i <= NRMLRT; i++)
    for (j = 1; j <= 3; j++)
      for (k = 1; k <= 3; k++)
	pibdiibdj[i][j][k] = 0;

  /* in 0-1-2 IBD to 0-1-2 IBD order */

  /* fullsib */
  pibdiibdj[1][1][1] = pow((ta*ta+(1-ta)*(1-ta)),2);
  pibdiibdj[1][1][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][1][3] = pow((2*ta*(1-ta)),2);
  pibdiibdj[1][2][1] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][2][2] = pow((ta*ta+(1-ta)*(1-ta)),2) + pow((2*ta*(1-ta)),2);
  pibdiibdj[1][2][3] = (ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][3][1] = pow((2*ta*(1-ta)),2);
  pibdiibdj[1][3][2] = 2*(ta*ta+(1-ta)*(1-ta))*2*ta*(1-ta);
  pibdiibdj[1][3][3] = pow((ta*ta+(1-ta)*(1-ta)),2);

  /* half-sib */
  pibdiibdj[2][1][1] = ta*ta+(1-ta)*(1-ta);
  pibdiibdj[2][1][2] = 2*ta*(1-ta);
  pibdiibdj[2][2][1] = 2*ta*(1-ta);
  pibdiibdj[2][2][2] = ta*ta+(1-ta)*(1-ta);

  /* grandparent-grandchild */
  pibdiibdj[3][1][1] = (1-ta);
  pibdiibdj[3][1][2] = ta;
  pibdiibdj[3][2][1] = ta;
  pibdiibdj[3][2][2] = (1-ta);

  /* avuncular */
  pibdiibdj[4][2][2] = 1-2.5*ta+4*ta*ta-2*ta*ta*ta;
  pibdiibdj[4][2][1] = 1-pibdiibdj[4][2][2];
  pibdiibdj[4][1][2] = 1-pibdiibdj[4][2][2];
  pibdiibdj[4][1][1] = pibdiibdj[4][2][2];

  /* first-cousin */
  pibdiibdj[5][2][2] = 1-4*ta+7.5*ta*ta- 6*ta*ta*ta+2*ta*ta*ta*ta;
  pibdiibdj[5][2][1] = 1-pibdiibdj[5][2][2];
  pibdiibdj[5][1][2] = (1-pibdiibdj[5][2][2])/3.0;
  pibdiibdj[5][1][1] = 1-pibdiibdj[5][1][2];

  /* unrelated */
  pibdiibdj[6][1][1] = 1;

  /* half-avuncular */
  pibdiibdj[7][2][2]=(1-ta)*(1-ta)*
    (1-ta)+ta*ta*(1-ta);
  pibdiibdj[7][2][1]=1-pibdiibdj[7][2][2];
  pibdiibdj[7][1][2]=(1-pibdiibdj[7][2][2])/3.0;
  pibdiibdj[7][1][1]=1-pibdiibdj[7][1][2];

  /* half-cousin */
  pibdiibdj[8][2][2]=(1-ta)*(1-ta)*(1-ta)*(1-ta)+ta*ta*(1-ta)*(1-ta);
  pibdiibdj[8][2][1]=1-pibdiibdj[8][2][2];
  pibdiibdj[8][1][2]=(1-pibdiibdj[8][2][2])/7.0;
  pibdiibdj[8][1][1]=1-pibdiibdj[8][1][2];

  /* half-sib plus first-cousin */
  pibdiibdj[9][3][3] = pibdiibdj[5][2][2]*phi;
  pibdiibdj[9][3][1] = pibdiibdj[5][2][1]*psi;
  pibdiibdj[9][3][2] = 1-pibdiibdj[9][3][3]-pibdiibdj[9][3][1];

  pibdiibdj[9][1][3] = pibdiibdj[5][1][2]*psi;
  pibdiibdj[9][1][1] = pibdiibdj[5][1][1]*phi;
  pibdiibdj[9][1][2] = 1-pibdiibdj[9][1][3]-pibdiibdj[9][1][1];

  pibdiibdj[9][2][3] = 1.0/3*pibdiibdj[5][2][2]*psi+
    2.0/3*pibdiibdj[5][1][2]*phi;
  pibdiibdj[9][2][1] = 1.0/3*pibdiibdj[5][2][1]*phi+
    2.0/3*pibdiibdj[5][1][1]*psi;
  pibdiibdj[9][2][2] = 1-pibdiibdj[9][2][3]-pibdiibdj[9][2][1];

  /* parent-offspring */
  pibdiibdj[10][2][2] = 1;

  /* MZ twins */
  pibdiibdj[11][3][3] = 1;

}

/* return P(Z<=z) for Z a N(0,1) RV */
double pnorm(double z)
{
  double x = 0.0;

  if (z != 0.0) {
    double y = 0.5 * fabs (z);
    if (y >= 3.0)
      x = 1.0;
    else if (y < 1.0) {
      double w = y*y;
      x = ((((((((0.000124818987 * w
        -0.001075204047) * w +0.005198775019) * w
        -0.019198292004) * w +0.059054035642) * w
        -0.151968751364) * w +0.319152932694) * w
        -0.531923007300) * w +0.797884560593) * y * 2.0;
    }
    else {
      y -= 2.0;
      x = (((((((((((((-0.000045255659 * y
        +0.000152529290) * y -0.000019538132) * y
        -0.000676904986) * y +0.001390604284) * y
        -0.000794620820) * y -0.002034254874) * y
        +0.006549791214) * y -0.010557625006) * y
        +0.011630447319) * y -0.009279453341) * y
        +0.005353579108) * y -0.002141268741) * y
        +0.000535310849) * y +0.999936657524;
    }
  }

  return (z > 0.0 ? ((x + 1.0) * 0.5) : ((1.0 - x) * 0.5));

}


/* map functions */

double mapthetatocm(double theta, int maptype)
{
  double cm = 0;

  /* no-interference model */
  if ( maptype == 1 ) {
    cm = -0.5*log(1-2*theta)*100;
  }
  /* kosambi model */
  else if ( maptype == 2 ) {
    cm = 0.25*log((1+2*theta)/(1-2*theta))*100;
  }

  return cm;
}

double mapcmtotheta(double cm, int maptype)
{
  double theta = 0;

  /* no-interference model */
  if ( maptype == 1 ) {
    theta = 0.5*(1-exp(-.02*cm));
  }
  /* kosambi model */
  else if ( maptype == 2 ) {
    theta = 0.5*(exp(.04*cm)-1)/(exp(.04*cm)+1);
  }

  return theta;
}


/* EM algorithm for estimation of p0, p1, p2 (0-1-2 IBD sharing) */

/*
 * Run one iteration of EM and return the next Phat values
 */
void EM_getnext(double phat[3 + 1], double p0, double p1, double a[],
    double b[], double c[], int commark)
{
  int i;
  double p2;
  double p0_next, p1_next, p2_next;

  p2 = 1 - p0 - p1;

  p0_next = p1_next = p2_next = 0;

  for (i = 1; i <= commark; i++) {
    p0_next += a[i] * p0 / (a[i] * p0 + b[i] * p1 + c[i] * p2);
    p1_next += b[i] * p1 / (a[i] * p0 + b[i] * p1 + c[i] * p2);
  }

  p0_next = p0_next / commark;
  p1_next = p1_next / commark;
  p2_next = 1.0 - p0_next - p1_next;

  phat[1] = p0_next;
  phat[2] = p1_next;
  phat[3] = p2_next;
}


double logLikelyhood(double p0, double p1, int commark)
{
    double L = 0;
    for(int i=1;i<=commark;i++) {
        L = L + log( a[i]*p0 + b[i]*p1 + c[i]*(1-p1-p0));
    }

    return L;
}

/* need to input three coefficents, initial p0, p1, and p2 values, and
   the total marker numbers */
void EMforp012(geno_t* geno1, geno_t* geno2, int* valid,
	       double*** q, int nochrom, int* nomark,
	       int reltype, double phat[3+1])
{
  int i, j, index, commark, niter;
  int f1, m1, f2, m2;
  double p0, p1, p2;
  double p0_next, p1_next, p2_next;
  
  int converge10 = 0;

  if(MYEPS > 0.5) { phat[1] = phat[2] = phat[3] = -1; return; }
  
  /* initialized the value for phat */
  p0 = p1 = p2 = 0;
  
  p0_next = .2;
  p1_next = .6;
  p2_next = .2;

  niter = index = commark = 0;

  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {
      if ( valid[index+j] != 0 ) {

	++commark;

	f1 = geno1[2*(index+j)-1];
	m1 = geno1[2*(index+j)];
	f2 = geno2[2*(index+j)-1];
	m2 = geno2[2*(index+j)];

	a[commark] = cp0(f1, m1, f2, m2, q[i][j]);
	b[commark] = cp1(f1, m1, f2, m2, q[i][j]);
	c[commark] = cp2(f1, m1, f2, m2, q[i][j]);

	/* for parent-offspring, use modified a b c if all of them
	   are 0 in the presence of genotyping error, so that
	   the denomiator a[i]*p0+b[i]*p1+c[i]*p2 will not be 0 in the EM */
	if (reltype == 10 && b[commark] < MYEPS)
        // if(1)
        {
	  b[commark] = (1-errorrate)*(1-errorrate)*b[commark]
	    +(1-(1-errorrate)*(1-errorrate))*a[commark];
	  c[commark] = (1-errorrate)*(1-errorrate)*c[commark]
	    +(1-(1-errorrate)*(1-errorrate))*a[commark];
	}

      }
    }

    index += nomark[i];
  }


  
  // The following is some debugging output
  if(debug == 3) {
    printf(" /* DATA for this pair */\n");
    printf("double p0=%g,p1=%g,p2=%g;\nint commark = %d;\n\n",phat[1],phat[2],phat[3],commark);
    printf("// niter=%d\n",niter);

    printf("double a[%d] = {",commark+10) ; for (i = 0; i <= commark; i++) printf("%7.5f , ",a[i]); printf(" }; \n");
    printf("double b[%d] = {",commark+10) ;for (i = 0; i <= commark; i++) printf("%7.5f , ",b[i]); printf(" }; \n");
    printf("double c[%d] = {",commark+10) ;for (i = 0; i <= commark; i++) printf("%7.5f , ",c[i]); printf(" };\n");
    printf("\n\n");

    printf("# R DATA:\n");
    printf("p0=%g; p1=%g; p2=%g; commark = %d;\n\n",phat[1],phat[2],phat[3],commark);
    printf("a<-c(",commark) ; for (i = 1; i <= commark; i++) printf("%7.5f%s ",a[i],(i==commark)?"":","); printf(" ); \n");
    printf("b<-c(",commark) ; for (i = 1; i <= commark; i++) printf("%7.5f%s ",b[i],(i==commark)?"":","); printf(" ); \n");
    printf("c<-c(",commark) ; for (i = 1; i <= commark; i++) printf("%7.5f%s ",c[i],(i==commark)?"":","); printf(" ); \n");
    printf("#\n");

  }


  //Use newton's method
  if(ibd_useNewton == 1) {
      EM_Newton em;

      double p0 = 0.6, p1 = 0.2;
      
      em.init(a, b, c, commark);
      int ok = em.EM(p0, p1, MYEPS, MITER);

      // Check if EM converged
      if(p0 > -1) {
          phat[1] = p0;
          phat[2] = p1;
          phat[3] = 1-p1-p0;

          if(phat[3] < 0) phat[3] = 0;

          //printf("Lik %f %f %f %f %f\n", p0, p1, logLikelyhood(1,0,commark), logLikelyhood(0.25,0.5,commark) , logLikelyhood(0.5,0.5,commark) );

          return;
      }

      printf("Warning: EM algorithm didn't converge.\n");
  }

  /* Check if EM converges at 1,0 */

  if(ibd_check10 == 1)
  {
      double next[4], p0, p1;
      double e = 1e-3;

      EM_getnext(next,1-2*e,e,a,b,c,commark);

      p0 = next[1]; p1 = next[2];

      /* printf("%.10f,%.10ff -> %.10f,%.10f\n",1-2*e,e,p0,p1); */

      if(p0 > 1-2*e && p1 < e)
      {
	  e = 1e-5;
	  EM_getnext(next,1-2*e,e,a,b,c,commark);
	  p0 = next[1]; p1 = next[2];
	  /* printf("%.10f,%.10ff -> %.10f,%.10f\n",1-2*e,e,p0,p1); */

	  if(p0 > 1-2*e && p1 < e) converge10 = 1;
      }
  }

  if(converge10 == 1) {
      p0_next = p0 = 1;
      p1_next = p1 = 0;
      p2_next = p2 = 0;
  }


  while( (fabs(p0_next-p0) > MYEPS || fabs(p1_next-p1) > MYEPS ) &&
	 niter <= MITER ) {
    p0 = p0_next;
    p1 = p1_next;
    p2 = p2_next;
    p0_next = p1_next = p2_next = 0;

    for (i = 1; i <= commark; i++) {

      p0_next += a[i]*p0/(a[i]*p0+b[i]*p1+c[i]*p2);
      p1_next += b[i]*p1/(a[i]*p0+b[i]*p1+c[i]*p2);

    }

    p0_next = p0_next/commark;
    p1_next = p1_next/commark;

    // Experimental series convergence acceleration procedure

    if(niter < -500) {
        double ap0 = p0 + 10*(p0_next - p0);
        double ap1 = p1 + 10*(p1_next - p1);
        if(ap0 + ap1 < 1 && ap0 > 0 && ap1 > 0) { p0_next = ap0; p1_next = ap1; }
    }

    p2_next = 1.0 - p0_next - p1_next;
    ++niter;
  }

  if(debug > 1) printf("EM: Converged after %d iterations , %.6g, %.6g\n", niter, p0_next, p1_next);

  if (niter <= MITER) {
    phat[1] = p0_next;
    phat[2] = p1_next;
    phat[3] = p2_next;
    /* to avoid print -0 */
    if(fabs(phat[3]) < MYEPS) phat[3] = 0.0;
  }
  else {
    phat[1] = -1;
    phat[2] = -1;
    phat[3] = -1;
  }

}


/* EIBD, AIBS, IBS statistic, mu and std calculation */

void get3stat(geno_t* geno1, geno_t* geno2, int* valid, int commark,
	      int reltype, int nochrom, int* nomark,
	      double p012[4], double*** q,
	      double statresult[3+1], double pibs[3])
{
  int i, n, kkk;
  int f1, f2, m1, m2;
  double p0, p1, p2;
  double eibd, aibs, ibs;


  kkk = 0;
  eibd = aibs = ibs = 0;
  pibs[0] = pibs[1] = pibs[2] = 0;

  for (n = 1; n <= nochrom; n++) {
    for (i = 1; i <= nomark[n]; i++) {
      if ( valid[kkk+i] != 0 ) {
	f1 = geno1[2*(kkk+i)-1];
	m1 = geno1[2*(kkk+i)];
	f2 = geno2[2*(kkk+i)-1];
	m2 = geno2[2*(kkk+i)];
 	/* eibd */
	p2 = cp2(f1, m1, f2, m2, q[n][i]);
	p1 = cp1(f1, m1, f2, m2, q[n][i]);
	p0 = cp0(f1, m1, f2, m2, q[n][i]);

        /* if it is not parent-offspring nor unrelated, nor MZ twins */
        if( !(reltype == 10 || reltype == 6 || reltype == 11) )
          eibd += (2*p012[3]*p2+p012[2]*p1)/
            (p012[3]*p2+p012[2]*p1+p012[1]*p0);

	/* aibs and ibs */
        if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]) +
	    p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          pibs[2] += 1;
          ibs += 2;
        }
        else if ( f1 == f2 || f1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]);
          pibs[1] += 1;
          ibs += 1;
        }
        else if ( m1 == f2 || m1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          pibs[1] += 1;
          ibs += 1;
        }
        else
          pibs[0] += 1;
      }
    }
    kkk += nomark[n];
  }

  pibs[0] /= commark;
  pibs[1] /= commark;
  pibs[2] /= commark;

  statresult[1] = eibd/commark;
  if(reltype == 6)
    statresult[1] = 0;
  else if(reltype == 10)
    statresult[1] = 1;
  else if(reltype == 11)
    statresult[1] = 2;
  statresult[2] = aibs/commark;
  if(reltype == 6)
    statresult[2] = 0;
  statresult[3] = ibs/commark;

}

void get3statandmustd(geno_t* geno1, geno_t* geno2, int* valid, int commark,
		      int reltype, int nochrom, int* nomark,
		      double p012[4],
		      double*** q, double** cenm, int maptype,
		      double*** eeibdibd, double*** eeibdeibdibd,
		      double*** eaibsibd, double*** eaibsaibsibd,
		      double** eibsibd, double** eibsibsibd,
		      double statresult[3+1][5+1], double pibs[3])
{
  int i, j, k, l, n, kkk;
  int f1, f2, m1, m2;
  double dij, p2, p1, p0;
  double eibd, veibd, aibs, vaibs, ibs, vibs;
  double dtemp[3+1], dtemp1[3+1], pibdiibdj[NRMLRT+1][3+1][3+1];
  double **exxchr, **exchr;

  exxchr = (double **)malloc((size_t) ((3+1)*sizeof(double*)));
  exchr = (double **)malloc((size_t) ((3+1)*sizeof(double*)));
  for (i = 1; i <= 3; i++) {
    exxchr[i] = (double *)malloc((size_t) ((nochrom+1)*sizeof(double)));
    exchr[i] = (double *)malloc((size_t) ((nochrom+1)*sizeof(double)));
  }

  eibd = aibs = ibs = 0;
  veibd = vaibs = vibs = 0;
  pibs[0] = pibs[1] = pibs[2] = 0;

  for (i = 1; i <= 3; i++) {
    dtemp1[i] = 0;
    for (n = 1; n <= nochrom; n++)
      exxchr[i][n] = exchr[i][n] = 0;
  }

  kkk = 0;
  for (n = 1; n <= nochrom; n++) {
    for (i = 1; i <= nomark[n]; i++) {
      /* statistic calculation */
      if ( valid[kkk+i] != 0 ) {
	f1 = geno1[2*(kkk+i)-1];
	m1 = geno1[2*(kkk+i)];
	f2 = geno2[2*(kkk+i)-1];
	m2 = geno2[2*(kkk+i)];
 	/* eibd */
	p2 = cp2(f1,m1,f2,m2,q[n][i]);
	p1 = cp1(f1,m1,f2,m2,q[n][i]);
	p0 = cp0(f1,m1,f2,m2,q[n][i]);

	/* if it is not parent-offspring nor unrelated, nor MZ twins */
	if( !(reltype == 10 || reltype == 6 || reltype == 11) )
	  eibd += (2*p012[3]*p2+p012[2]*p1)/
	    (p012[3]*p2+p012[2]*p1+p012[1]*p0);
        
	/* aibs and ibs */
        if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]) +
	    p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          ibs += 2;
          pibs[2] += 1;
        }
        else if ( f1 == f2 || f1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][f1]);
          ibs += 1;
          pibs[1] += 1;
        }
        else if ( m1 == f2 || m1 == m2 ) {
          aibs += p012[0]/(p012[0] + (1-p012[0])*q[n][i][m1]);
          ibs += 1;
          pibs[1] += 1;
        }
        else
          pibs[0] += 1;
      }

      /* variance and expection calculation */
      for (j = i; j <= nomark[n]; j++) {

	if (i == j && valid[kkk+i] == 1) {
	  for (k = 1; k <= 3; k++)
	    dtemp[k] = 0;
	  for (k = 1; k <= 3; k++) {
	    dtemp[1] += eeibdeibdibd[reltype][kkk+i][k]*p012[k];
            dtemp[2] += eaibsaibsibd[reltype][kkk+i][k]*p012[k];
	    dtemp[3] += eibsibsibd[kkk+i][k]*p012[k];
 	  }
	  for (k = 1; k <= 3; k++)
	    exxchr[k][n] += dtemp[k];
	}

	else if ( j > i && valid[kkk+i] == 1 && valid[kkk+j] == 1) {
	  dij = 0;
	  for (k = (i+1); k <= j; k++)
	    dij += cenm[n][k];
          
	  transprob(dij, pibdiibdj, reltype);

	  for (k = 1; k <= 3; k++) dtemp[k] = 0;
	  for (k = 1; k <= 3; k++)
	    for (l = 1; l <=3; l++) {
	      dtemp[1] += eeibdibd[reltype][kkk+i][k]*
		eeibdibd[reltype][kkk+j][l]*p012[k]*
		pibdiibdj[reltype][k][l];
	      dtemp[2] += eaibsibd[reltype][kkk+i][k]*
		eaibsibd[reltype][kkk+j][l]*p012[k]*
		pibdiibdj[reltype][k][l];
	      dtemp[3] += eibsibd[kkk+i][k]*
		eibsibd[kkk+j][l]*p012[k]*
		pibdiibdj[reltype][k][l];
	    }
	  for (l = 1; l <=3; l++)
	    exxchr[l][n] += 2*dtemp[l];
	}
      }

      if ( valid[kkk+i] == 1 ) {
	exchr[1][n] += 4*p012[0];
	dtemp[2] = 0;
	for (k = 1; k <= 3; k++)
          dtemp[2] += eaibsibd[reltype][kkk+i][k]*p012[k];
        exchr[2][n] += dtemp[2];
	dtemp[3] = 0;
        for (k = 1; k <= 3; k++)
          dtemp[3] += eibsibd[kkk+i][k]*p012[k];
        exchr[3][n] += dtemp[3];
      }

    }

    for (k =  1; k <= 3; k++)
      dtemp1[k] += exchr[k][n];

    if (!(reltype == 10 || reltype == 6 || reltype == 11) )
      veibd += exxchr[1][n] - exchr[1][n]*exchr[1][n];
    if (!(reltype == 6) )
      vaibs += exxchr[2][n] - exchr[2][n]*exchr[2][n];
    vibs += exxchr[3][n] - exchr[3][n]*exchr[3][n];
    kkk += nomark[n];

  }

  pibs[0] /= commark;
  pibs[1] /= commark;
  pibs[2] /= commark;

  statresult[1][1] = eibd/commark;
  if(reltype == 6)
    statresult[1][1] = 0;
  else if(reltype == 10)
    statresult[1][1] = 1;
  else if(reltype == 11)
    statresult[1][1] = 2;

  /* E(EIBD) = kinship-coef*4 */
  statresult[1][2] =  p012[0]*4;
  /* delete later */
  if (fabs(dtemp1[1]/commark-statresult[1][2])>0.001) {
    printf("program is wrong for EIBD\n");
    exit(0);
  }

  statresult[1][3] = sqrt(veibd/(commark*commark));
  if(reltype == 6 || reltype == 10 || reltype == 11)
    statresult[1][3] = 0;

  statresult[1][4] = (statresult[1][1] - statresult[1][2])/
    statresult[1][3];
  statresult[1][5] = 2*pnorm(-fabs(statresult[1][4]));

  statresult[2][1] = aibs/commark;
  if(reltype == 6)
    statresult[2][1] = 0;

  statresult[2][2] = dtemp1[2]/commark;
  statresult[2][3] = sqrt(vaibs/(commark*commark));
  if(reltype == 6)
    statresult[2][3] = 0;

  statresult[2][4] = (statresult[2][1] - statresult[2][2])/
    statresult[2][3];
  statresult[2][5] = 2*pnorm(-fabs(statresult[2][4]));

  statresult[3][1] = ibs/commark;
  statresult[3][2] = dtemp1[3]/commark;
  statresult[3][3] = sqrt(vibs/(commark*commark));
  statresult[3][4] = (statresult[3][1] - statresult[3][2])/
    statresult[3][3];
  statresult[3][5] = 2*pnorm(-fabs(statresult[3][4]));

  for (i = 1; i <= 3; i++) {
    free(exchr[i]);
    free(exxchr[i]);
  }

  free(exchr);
  free(exxchr);
}



void forEIBDvar(double*** eeibdibd, double*** eeibdeibdibd,
		int nochrom, int* nomark, int* ngeno, int*** genotype,
		double*** q, double p012[NRMLRT+1][3+1])
{
  int i, j, k, l, m, kkk;
  int f1, f2, m1, m2;
  double p2, p1, p0;
  double eibd;

  kkk=0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {

      for (k = 1; k <= NR; k++) {
	for (l = 1; l <= 3; l++)
	  eeibdibd[k][kkk+j][l] = eeibdeibdibd[k][kkk+j][l] = 0;
      }

      for (k = 1; k <= ngeno[kkk+j]; k++) {
	for (l = k; l <= ngeno[kkk+j]; l++) {
	  f1 = genotype[kkk+j][k][1];
	  m1 = genotype[kkk+j][k][2];
	  f2 = genotype[kkk+j][l][1];
	  m2 = genotype[kkk+j][l][2];

	  p2 = cp2(f1,m1,f2,m2,q[i][j]);
	  p1 = cp1(f1,m1,f2,m2,q[i][j]);
	  p0 = cp0(f1,m1,f2,m2,q[i][j]);

	  for (m = 1; m <= NR; m++) {
	    eibd = (2*p012[m][3]*p2+p012[m][2]*p1)/
	      (p012[m][3]*p2+p012[m][2]*p1+p012[m][1]*p0);

	    eeibdibd[m][kkk+j][3] += eibd*p2;
	    eeibdeibdibd[m][kkk+j][3] += eibd*eibd*p2;

	    eeibdibd[m][kkk+j][2] += eibd*p1;
	    eeibdeibdibd[m][kkk+j][2] += eibd*eibd*p1;

	    eeibdibd[m][kkk+j][1] += eibd*p0;
	    eeibdeibdibd[m][kkk+j][1] += eibd*eibd*p0;
	  }
	}
      }
    }
    kkk += nomark[i];
  }

}

void forAIBSvar(double*** eaibsibd, double*** eaibsaibsibd,
		int nochrom, int* nomark, int* ngeno, int*** genotype,
		double*** q, double p012[NRMLRT+1][3+1])
{
  int i, j, k, l, m, kkk;
  int f1, f2, m1, m2;
  double p2, p1, p0;
  double aibs;

  kkk=0;

  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {

      for (k = 1; k <= NR; k++) {
	for (l = 1; l <= 3; l++)
	eaibsibd[k][kkk+j][l] = eaibsaibsibd[k][kkk+j][l] = 0;
      }

      for (k = 1; k <= ngeno[kkk+j]; k++) {
	for (l = k; l <= ngeno[kkk+j]; l++) {

	  f1 = genotype[kkk+j][k][1];
	  m1 = genotype[kkk+j][k][2];
	  f2 = genotype[kkk+j][l][1];
	  m2 = genotype[kkk+j][l][2];

	  p2 = cp2(f1,m1,f2,m2,q[i][j]);
	  p1 = cp1(f1,m1,f2,m2,q[i][j]);
	  p0 = cp0(f1,m1,f2,m2,q[i][j]);

	  for (m = 1; m <= NR; m++) {
	    if ( (f1 == f2 && m1 == m2) || (f1 == m2 && m1 == f2) )
	      aibs =
		p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][f1]) +
		p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][m1]);
	    else if ( f1 == f2 || f1 == m2 )
	      aibs = p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][f1]);
	    else if ( m1 == f2 || m1 == m2 )
	      aibs = p012[m][0]/(p012[m][0]+(1-p012[m][0])*q[i][j][m1]);
	    else
	      aibs = 0;

	    eaibsibd[m][kkk+j][3] += aibs*p2;
	    eaibsaibsibd[m][kkk+j][3] += aibs*aibs*p2;

            eaibsibd[m][kkk+j][2] += aibs*p1;
            eaibsaibsibd[m][kkk+j][2] += aibs*aibs*p1;

            eaibsibd[m][kkk+j][1] += aibs*p0;
            eaibsaibsibd[m][kkk+j][1] += aibs*aibs*p0;
	  }
	}
      }
    }
    kkk += nomark[i];
  }

}

void forIBSvar(double** eibsibd, double** eibsibsibd,
	       int nochrom, int* nomark, double*** q, int** noalle)
{
  int i, j, k, l, m, kkk;
  double dtemp, dtemp1, dtemp2, dtemp3, dtemp4;

  kkk=0;
  for (i = 1; i <= nochrom; i++) {
    for (j = 1; j <= nomark[i]; j++) {

      eibsibd[kkk+j][3] = 2;
      eibsibsibd[kkk+j][3] = 4;

      dtemp = 0;
      for (k = 1; k <= noalle[i][j]; k++)
        dtemp += q[i][j][k]*q[i][j][k];
      eibsibd[kkk+j][2] = dtemp*2 + (1-dtemp);
      eibsibsibd[kkk+j][2] = dtemp*4 + (1-dtemp);

      dtemp1 = dtemp2 = dtemp3 = dtemp4 = 0;
      for (k = 1; k <= noalle[i][j]; k++) {
        dtemp1 += q[i][j][k]*q[i][j][k]*q[i][j][k]*q[i][j][k];
        for (l = 1; l <= noalle[i][j];l++) {
          if (l != k) {
            dtemp2 += 2*q[i][j][k]*q[i][j][k]*q[i][j][l]*q[i][j][l];
            dtemp3 += 4*q[i][j][k]*q[i][j][k]*q[i][j][k]*q[i][j][l];
            for (m = 1; m <= noalle[i][j]; m++) {
              if (m != k && m != l)
                dtemp4 +=
		  4*q[i][j][k]*q[i][j][k]*q[i][j][l]*q[i][j][m];
            }
          }
        }
      }
      eibsibd[kkk+j][1] = 2*(dtemp1+dtemp2)+(dtemp3+dtemp4);
      eibsibsibd[kkk+j][1] = 4*(dtemp1+dtemp2)+(dtemp3+dtemp4);
    }
    kkk += nomark[i];
  }

}

/* for likelihood calculation */
void getapprloglikelihood(geno_t* geno1, geno_t* geno2, int* valid,
			  double*** q, double** cenm,
			  int nochrom, int* nomark, int maptype,
			  double loglhresult[NRMLRT+1])
{
  int i, k, n, kkk;
  int validmin, validmax;
  int f1, m1, f2, m2;
  double s0, s1, s2, h0, h1, g0, g1;
  double au0, au1, ac0, ac1, unr0;
  double s0_next, s1_next, s2_next, h0_next, h1_next, g0_next, g1_next;
  double au0_next, au1_next, ac0_next, ac1_next, unr0_next;
  double po1, po1_next, mz2, mz2_next;

  double ahu0, ahu1, ahc0, ahc1, ahplusc0, ahplusc1, ahplusc2;
  double ahu0_next, ahu1_next, ahc0_next, ahc1_next;
  double ahplusc0_next, ahplusc1_next, ahplusc2_next;

  double p0=0, p1=0, p2=0;
  double dij, pij[NRMLRT+1][3+1][3+1];


  kkk = 0;
  for (i = 1; i <= NRMLRT; i++)
    loglhresult[i] = 0;

  for (n = 1; n <= nochrom; n++) {
    /* initialization for HMM for LRT */

    /* full-sib */
    s0 = s2 = 0.25;
    s1 = 0.5;
    /* half-sib */
    h0 = h1 = 0.5;
    /* grandparent-granchild */
    g0 = g1 = 0.5;
    /* avuncular, use approximation */
    au0 = au1 = 0.5;
    /* first-cousin, use approximation */
    ac0 = 0.75;
    ac1 = 0.25;
    /* unrelated */
    unr0 = 1;
    /* half-avuncular */
    ahu0 = 0.75;
    ahu1 = 0.25;
    /* half-cousin */
    ahc0 = 0.875;
    ahc1 = 0.125;
    /* halfsib plus firstcousin */
    ahplusc0 = 0.375;
    ahplusc1 = 0.5;
    ahplusc2 = 0.125;
    /* parent-offspring */
    po1 = 1;
    /* MZ twins */
    mz2 = 1;

    /* get the location for the first typed marker and the last typed
       marker on each chromosome */
    k = 1;
    while( valid[kkk+k] == 0 )
      ++k;
    validmin = k;
    k = nomark[n];
    while( valid[kkk+k] == 0 )
      k = k - 1;
    validmax = k;

    /* for each chromosome */
    for (i = validmin; i <= validmax; i++) {
      if ( valid[kkk+i] == 1 ) {

	f1 = geno1[2*(kkk+i)-1];
        m1 = geno1[2*(kkk+i)];
        f2 = geno2[2*(kkk+i)-1];
        m2 = geno2[2*(kkk+i)];

	p2 = cp2(f1, m1, f2, m2, q[n][i]);
        p1 = cp1(f1, m1, f2, m2, q[n][i]);
        p0 = cp0(f1, m1, f2, m2, q[n][i]);

	if ( i != validmax) {
	  /* get the distrance between i_th marker and the next typed
	     marker */
	  k = i + 1;
	  dij = cenm[n][k];
	  while (valid[kkk+k] == 0) {
	    ++k;
	    dij += cenm[n][k];
	  }

	  /* get the transition probabilty for all the reltype */
	  transproball(dij, pij);

	  s0_next = s0*p0*pij[1][1][1] + s1*p1*pij[1][2][1] +
	    s2*p2*pij[1][3][1];
	  s1_next = s0*p0*pij[1][1][2] + s1*p1*pij[1][2][2] +
	    s2*p2*pij[1][3][2];
	  s2_next = s0*p0*pij[1][1][3] + s1*p1*pij[1][2][3] +
	    s2*p2*pij[1][3][3];
	  h0_next = h0*p0*pij[2][1][1] + h1*p1*pij[2][2][1];
	  h1_next = h0*p0*pij[2][1][2] + h1*p1*pij[2][2][2];
	  g0_next = g0*p0*pij[3][1][1] + g1*p1*pij[3][2][1];
	  g1_next = g0*p0*pij[3][1][2] + g1*p1*pij[3][2][2];
	  au0_next = au0*p0*pij[4][1][1] + au1*p1*pij[4][2][1];
	  au1_next = au0*p0*pij[4][1][2] + au1*p1*pij[4][2][2];
	  ac0_next = ac0*p0*pij[5][1][1] + ac1*p1*pij[5][2][1];
	  ac1_next = ac0*p0*pij[5][1][2] + ac1*p1*pij[5][2][2];

	  unr0_next = unr0*p0;

	  ahu0_next = ahu0*p0*pij[7][1][1] + ahu1*p1*pij[7][2][1];
	  ahu1_next = ahu0*p0*pij[7][1][2] + ahu1*p1*pij[7][2][2];
	  ahc0_next = ahc0*p0*pij[8][1][1] + ahc1*p1*pij[8][2][1];
	  ahc1_next = ahc0*p0*pij[8][1][2] + ahc1*p1*pij[8][2][2];
	  ahplusc0_next = ahplusc0*p0*pij[9][1][1] +
	    ahplusc1*p1*pij[9][2][1] + ahplusc2*p2*pij[9][3][1];
	  ahplusc1_next = ahplusc0*p0*pij[9][1][2] +
	    ahplusc1*p1*pij[9][2][2] + ahplusc2*p2*pij[9][3][2];
	  ahplusc2_next = ahplusc0*p0*pij[9][1][3] +
	    ahplusc1*p1*pij[9][2][3] + ahplusc2*p2*pij[9][3][3];

	  /* for parent-offspring, and MZ twins, take into account of
	     genotyping error */
	  po1_next =
	    po1*((1-errorrate)*(1-errorrate)*p1+(1-(1-errorrate)*(1-errorrate))*p0);

	  mz2_next =
	    mz2*((1-errorrate)*(1-errorrate)*p2+(1-(1-errorrate)*(1-errorrate))*p0);


	  s0 = s0_next;
	  s1 = s1_next;
	  s2 = s2_next;
	  h0 = h0_next;
	  h1 = h1_next;
	  g0 = g0_next;
	  g1 = g1_next;
	  au0 = au0_next;
	  au1 = au1_next;
	  ac0 = ac0_next;
	  ac1 = ac1_next;

	  unr0 = unr0_next;

	  ahu0 = ahu0_next;
	  ahu1 = ahu1_next;
	  ahc0 = ahc0_next;
	  ahc1 = ahc1_next;
	  ahplusc0 = ahplusc0_next;
	  ahplusc1 = ahplusc1_next;
	  ahplusc2 = ahplusc2_next;

	  po1 = po1_next;
	  mz2 = mz2_next;

	}
      }
    }
    loglhresult[1] += log(s0*p0+s1*p1+s2*p2);
    loglhresult[2] += log(h0*p0+h1*p1);
    loglhresult[3] += log(g0*p0+g1*p1);
    loglhresult[4] += log(au0*p0+au1*p1);
    loglhresult[5] += log(ac0*p0+ac1*p1);
    loglhresult[6] += log(unr0*p0);

    loglhresult[7] += log(ahu0*p0+ahu1*p1);
    loglhresult[8] += log(ahc0*p0+ahc1*p1);
    loglhresult[9] += log(ahplusc0*p0+ahplusc1*p1+ahplusc2*p2);

    loglhresult[10] +=
      log(po1*((1-errorrate)*(1-errorrate)*p1+(1-(1-errorrate)*(1-errorrate))*p0));

    loglhresult[11] +=
      log(mz2*((1-errorrate)*(1-errorrate)*p2+(1-(1-errorrate)*(1-errorrate))*p0));

    kkk += nomark[n];
  }
}


/* for simulate genodata of the relaitonships considered */

void simugeno1geno2(double*** q, double*** qcum, int nochrom,
		    int* nomark, int **noalle, int reltype,
		    double* chromlength, double** cenm,
		    geno_t* simugeno1, geno_t* simugeno2)
{
  int i, j, k, kkk;
  int *ibd1, *ibd2, ibd;
  int nc1, nc2, ncs, ncnew1, ncnew2;
  double *pos1 = new_double1d(MC+1);
  double *pos2 = new_double1d(MC+1);
  double *poss = new_double1d(MC+1);
  double *posnew1 = new_double1d(MC+1);
  double *posnew2 = new_double1d(MC+1);
  char *ind1 = new_char1d(MC+1);
  char *ind2 = new_char1d(MC+1);
  char *inds = new_char1d(MC+1);
  char *indnew1 = new_char1d(MC+1);
  char *indnew2 = new_char1d(MC+1);

  /*double pos1[MC+1], pos2[MC+1], poss[MC+1];
  double posnew1[MC+1], posnew2[MC+1];
  char ind1[MC+1], ind2[MC+1], inds[MC+1];
  char indnew1[MC+1], indnew2[MC+1];*/
  double sign;


  kkk = 0;
  for (i = 1; i <= nochrom; i++) {

    ibd1 = (int *)malloc((size_t) ((nomark[i]+1)*sizeof(int)));
    ibd2 = (int *)malloc((size_t) ((nomark[i]+1)*sizeof(int)));

    for (j = 1; j <= nomark[i]; j++)
      ibd1[j] = ibd2[j] = 0;

    /* simulate full-sib */
    if ( reltype == 1 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '1', '2', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd1);

      nc1 = crossover_new((chromlength[i]), pos1,'3','4', ind1);
      nc2 = crossover_new((chromlength[i]), pos2,'3','4', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd2);
    }

    /*simulate half-sib */
    else if ( reltype == 2 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '1', '2', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }

    /* simulate grandparent-grandchild */
    else if ( reltype == 3 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs,
		       ind1, ind2, inds, posnew1, indnew1);
      getibd0(posnew1, ncnew1, indnew1, '1', '2',
	      chromlength[i],cenm[i], nomark[i], ibd1);
    }

    /* simulate avuncular */
    else if ( reltype == 4 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);

      getibd2(pos1, pos2, posnew1, nc1, nc2, ncnew1,
	      ind1, ind2, indnew1, chromlength[i], cenm[i], nomark[i],
	      ibd1);
    }

    /* simulate first-cousin */
    else if ( reltype == 5 ) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew2 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew2, indnew2);

      getibd1(posnew1, posnew2, ncnew1, ncnew2, indnew1, indnew2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }

    /* for unrelated reltype == 6 */
    /* no crossover simulation needed */

    /* for half-avuncular */
    else if (reltype == 7) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);

      getibd1(posnew1, pos1, ncnew1, nc1, indnew1, ind1,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }

    /* for half-cousin */
    else if (reltype == 8) {
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '5', '6', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew2 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew2, indnew2);

      getibd1(posnew1, posnew2, ncnew1, ncnew2, indnew1, indnew2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
    }

    /* for half-sib plus first-cousin */
    else if (reltype == 9) {

      /* half-sib part */
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '1', '2', ind2);
      getibd1(pos1, pos2, nc1, nc2, ind1, ind2,
	      chromlength[i], cenm[i], nomark[i], ibd1);
      /* first-cousin part */
      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew1 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew1, indnew1);

      nc1 = crossover_new((chromlength[i]), pos1, '1', '2', ind1);
      nc2 = crossover_new((chromlength[i]), pos2, '3', '4', ind2);
      ncs = crossover_new((chromlength[i]), poss, 'p', 'm', inds);
      ncnew2 = combine(pos1, pos2, poss, nc1, nc2, ncs, ind1, ind2,
		       inds, posnew2, indnew2);

      getibd1(posnew1, posnew2, ncnew1, ncnew2, indnew1, indnew2,
	      chromlength[i], cenm[i], nomark[i], ibd2);
    }

    /* now based on the ibd value for each marker,
       simulate allele data*/
    for (j = 1; j <= nomark[i]; j++) {

      /* if not parent-offspring or MZ twins*/
      if( !(reltype == 10 || reltype == 11) ) {
	ibd = ibd1[j] + ibd2[j];

	if ( ibd == 0 ) {

	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  sign=1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)] = k;
	  }
	}

	else if ( ibd == 1 ) {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}

	else {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  simugeno2[2*(kkk+j)-1] = simugeno1[2*(kkk+j)-1];
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}
      }

      else if (reltype == 10) {
	/* simulate data for parent-offspring with genotypying rate for
	   each marker = errorrate */
	sign = 1.0 * rand() / RAND_MAX;
	/* if error occurr, sample as if unrelated */
	if( sign <= (1-(1-errorrate)*(1-errorrate))) {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  sign=1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)] = k;
	  }
	}
	/* if no error, sample as IBD = 1 */
	else {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}
      }

      else if (reltype == 11) {
	/* simulate data for MZ twins with genotypying rate for
	   each marker = errorrate */
	sign = 1.0 * rand() / RAND_MAX;
	/* if error occurr, sample as if unrelated */
	if( sign <= (1-(1-errorrate)*(1-errorrate))) {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)-1] = k;
	  }
	  sign=1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno2[2*(kkk+j)] = k;
	  }
	}
	/* if no error, sample as IBD = 2 */
	else {
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)-1] = k;
	  }
	  sign = 1.0 * rand() / RAND_MAX;
	  for (k = 1; k <= noalle[i][j]; k++) {
	    if ( sign >= qcum[i][j][k-1] && sign <= qcum[i][j][k] )
	      simugeno1[2*(kkk+j)] = k;
	  }
	  simugeno2[2*(kkk+j)-1] = simugeno1[2*(kkk+j)-1];
	  simugeno2[2*(kkk+j)] = simugeno1[2*(kkk+j)];
	}
      }

    }
    kkk += nomark[i];

    free(ibd1);
    free(ibd2);
  }

  free(pos1);
  free(pos2);
  free(poss);
  free(posnew1);
  free(posnew2);
  free(ind1);
  free(ind2);
  free(inds);
  free(indnew1);
  free(indnew2);
}

/* input the length of chromosome and position vector
   output the number of crossover and update the position vector */
int crossover_new(double chromlength, double* pos, char c1, char c2,
                  char* ind)
{
  int  i, K, cn, cn_temp;
  double sign, d0, d_next, d_cum;
  double *pos_temp = new_double1d(MC+2);

  /* first crossover K iid exp(1/10), K is uniform{1,2,3,4,5} */
  sign = 1.0 * rand() / RAND_MAX;
  if ( sign == 1.0 )
    K = 5;
  else
    K = (int)floor(sign/0.2)+1;
  d0 = 0;
  for (i = 1; i <= K; i++) {
    sign = 1.0 * rand() / RAND_MAX;
    d0 += ((sign == 1.0) ? 100000 : (-0.1*log(1-sign)));
  }
  d0 = d0*100;
  /* following crossovers Gamma(5,1/10)= 5 iid exp(1/10) */
  cn_temp = 0;
  d_cum = 0;
  d_cum += d0;
  if ( d_cum <= chromlength ) {
    ++cn_temp;
    pos_temp[cn_temp] = d_cum;
  }
  else
    pos_temp[1] = d_cum;

  while( d_cum < chromlength && cn_temp > 0 ) {
    d_next = 0;
    K = 5;
    for (i = 1; i <= K; i++) {
      sign = 1.0 * rand() / RAND_MAX;
      d_next += ((sign == 1.0) ? 100000 : (-0.1*log(1-sign)));
    }
    d_next = d_next*100;
    d_cum += d_next;
    if ( d_cum <= chromlength ) {
      ++cn_temp;
      pos_temp[cn_temp] = d_cum;
    }
  }
  pos_temp[cn_temp+1] = d_cum;

  /* with 1/2 prob. delete each crossover point independently and
     reallocate the leftover points in the position vector */
  cn = 0;
  pos[0] = 0;
  for (i = 1; i <= cn_temp; i++) {
    sign = 1.0 * rand() / RAND_MAX;
    if ( sign >= 0.5 ) {
      ++cn;
      pos[cn] = pos_temp[i];
    }
  }
  pos[cn+1] = pos_temp[cn_temp+1];

  /* with 1/2 prob. decide the first index, then alternate the index
     for following markers */
  sign = 1.0 * rand() / RAND_MAX;
  if ( sign >= 0.5 ) {
    for (i = 0; i <= cn; i = i+2) {
      ind[i] = c1;
      ind[i+1] = c2;
    }
  }
  else{
    for (i = 0; i <= cn; i = i+2) {
      ind[i] = c2;
      ind[i+1] = c1;
    }
  }

  
  {
//    printf("crossover_new: %f (%c%c) ", chromlength, c1,c2);
 //   for (int i = 0; i <= cn; i ++) printf("%c ", ind[i]); printf("\n");
  }

  free(pos_temp);
  return cn;
}

/* combine two crossovered strings to a new one
   based on the crossover information */
int combine(double* p1, double* p2, double* ps, int nc1, int nc2, int ncs,
            char* i1, char* i2, char* is, double* pnew, char* inew)
{
  int currs, curr1, curr2, newpos;
  char currstr;

  curr1 = curr2 = newpos = 0;
  for (currs = 1; currs <= ncs; ++currs) {
    currstr = is[currs - 1];
    if ( currstr == 'p' ) {
      while ( curr1 != nc1 && p1[curr1 + 1] <= ps[currs - 1] )
        ++curr1;
      pnew[newpos] = ps[currs - 1];
      inew[newpos] = i1[curr1];
      newpos++;
      while ( curr1 != nc1 && p1[curr1 + 1] <= ps[currs] ) {
        ++curr1;
        pnew[newpos] = p1[curr1];
        inew[newpos] = i1[curr1];
        newpos++;
      }
    }
    else {
      while ( curr2 != nc2 && p2[curr2 + 1] <= ps[currs - 1] )
        ++curr2;
      pnew[newpos] = ps[currs - 1];
      inew[newpos] = i2[curr2];
      newpos++;
      while ( curr2 != nc2 && p2[curr2 + 1] <= ps[currs] ) {
        ++curr2;
        pnew[newpos] = p2[curr2];
        inew[newpos] = i2[curr2];
        newpos++;
      }
    }
  }
  if ( is[ncs] == 'p' ) {
    while ( curr1 <= nc1 && p1[curr1] <= ps[ncs] )
      ++curr1;
    curr1 = curr1 - 1;
    pnew[newpos] = ps[ncs];
    inew[newpos] = i1[curr1];
    newpos++;
    while ( curr1 < nc1 ) {
      ++curr1;
      pnew[newpos] = p1[curr1];
      inew[newpos] = i1[curr1];
      newpos++;
    }
    pnew[newpos] = p1[nc1+1];
  }
  else {
    while ( curr2 <= nc2 && p2[curr2] <= ps[ncs] )
      ++curr2;
    curr2 = curr2 - 1;
    pnew[newpos] = ps[ncs];
    inew[newpos] = i2[curr2];
    newpos++;
    while ( curr2 < nc2 ) {
      ++curr2;
      pnew[newpos] = p2[curr2];
      inew[newpos] = i2[curr2];
      newpos++;
    }
    pnew[newpos] = p2[nc2+1];
  }
  return newpos - 1;
}

/* getibd2(): compare one crossovered strings with the other two
   crossoved strings */
void getibd2(double* p1,double* p2, double* pnew, int nc1,int nc2,
	     int ncnew, char* i1,char* i2, char* inew,
	     double chromlength, double* cenm, int nomark,
	     int* ibd)
{
  int i, currind, curr1, curr2, currnew;
  double  ipos;

  currind = 1;
  curr1 = curr2 = currnew = 0;
  ipos = 0.0;

  for (i = 1; i <= nomark; i++) {
    while ( curr1 != nc1 && p1[curr1 + 1] <= ipos )
      ++curr1;
    while ( curr2 != nc2 && p2[curr2 + 1] <= ipos )
      ++curr2;
    while ( currnew != ncnew && pnew[currnew + 1] <= ipos )
      ++currnew;
    if ( i1[curr1] == inew[currnew] || i2[curr2] == inew[currnew] )
      ibd[currind] = 1;
    else
      ibd[currind] = 0;

    ++currind;

    if (i < nomark)
      ipos += cenm[i+1];
  }
}

/* getibd1(): compare two crossovered strings
   important: ibd starts from ibd[0] */
void getibd1(double* pnew1, double* pnew2, int ncnew1, int ncnew2,
	     char* inew1, char* inew2, double chromlength,
	     double* cenm, int nomark, int* ibd)
{
  int i, currind, curr1, curr2;
  double  ipos;

  currind = 1;
  curr1 = curr2 = 0;
  ipos = 0.0;

  for (i = 1; i <= nomark; i++) {

    while ( curr1 != ncnew1 && pnew1[curr1 + 1] <= ipos )
      ++curr1;
    while ( curr2 != ncnew2 && pnew2[curr2 + 1] <= ipos )
      ++curr2;
    if ( inew1[curr1] == inew2[curr2] )
      ibd[currind] = 1;
    else
      ibd[currind] = 0;

    ++currind;

    if (i < nomark)
      ipos += cenm[i+1];
  }
}

/* getibd0(): compare one crossovered string with one set of orignial
   two indexes */
void getibd0(double* pnew, int ncnew, char* inew, char c1, char c2,
	     double chromlength, double* cenm, int nomark,
	     int* ibd)
{
  int i, currind, currnew;
  double  ipos;

  currind = 1;
  currnew = 0;
  ipos = 0.0;

  for (i = 1; i <= nomark; i++) {
    while ( currnew != ncnew && pnew[currnew + 1] <= ipos )
      ++currnew;
    if ( inew[currnew]==c1 || inew[currnew]==c2 )
      ibd[currind] = 1;
    else
      ibd[currind] = 0;
    ++currind;

    if (i < nomark)
      ipos += cenm[i+1];
  }
}




#ifdef _GNU_SOURCE

/* fgets_adj: reads a line from input file *in into buffer s of size s_len
 *  adjusts size of buffer if not enough.
 */

// Use the very fast GNU getline when available
char* fgets_adj(char **s, int *s_len, FILE *in)
{
  size_t n = *s_len;

  ssize_t l = getline(s, &n, in);
  *s_len = n;

  if(l <= 0) return NULL;
  return *s;
}

#else

char* fgets_adj(char **s, int *s_len, FILE *in)
{
	int l;
	char *ret;

	fpos_t pos;
	fgetpos(in, &pos);

again:
	ret = fgets(*s, (*s_len), in);
	if(ret == NULL) return NULL;

	l = strlen(*s);

	if((*s)[l-1] != '\n' && l == (*s_len)-1) {
	      *s_len = 2* (*s_len);
	      *s = (char*)realloc(*s,(size_t)(*s_len));
	      fsetpos(in,&pos);

	      goto again;
	}

	return *s;
}

#endif


/*
 * Read next configuration file line and skip comments
 */

char* get_next_conf_line(FILE *fp, char *line, int line_size)
{
    int pos;
    char *ret = fgets (line, line_size, fp);

    if(ret == NULL) { line[0] = 0; return ret; }

    /* remove comments */
    for(pos = 0; line[pos] != 0; pos ++) {
        if(line[pos] == '\n') { line[pos] = 0; break; }
        if(line[pos] == '#') { line[pos] = 0; break; }
    }

    return ret;
}

/*
 * Check if file exists, report an error and exit if it doesn't.
 */

void file_exists(char *name)
{
    FILE *tmpfp;
    tmpfp = fopen(name, "r");
	
    if (tmpfp == NULL) {
        printf("\nError: Unable to open input file: '%s'\n", name);
        exit(1);
    }

    fclose(tmpfp);
}

/*
 * Read configuration file
 */
int read_config_file(char *name)
{
    char line[5000];
    int error = 0, ok = 0, line_no = 0, pos, i;

    FILE *f = fopen(name,"r");
    if(!f) {
        printf("Error: Cannot open configuration file %s\n", name);
        exit(1);
    }

    while(!feof(f)) {
        line_no++;
        get_next_conf_line(f,line,5000);

        /* skip empty lines */
        ok = 0;
        for(pos = 0; line[pos] != 0; pos ++)
          if(line[pos] != ' ' && line[pos] != '\t' && line[pos] != '\n') { ok = 1; break; }

        if(ok==0) { continue; }

        if (sscanf (line, " DEBUG %d", &debug) == 1) continue;

        if (sscanf (line, " WRITEFREQUENCIES %d", &write_frequency_file) == 1) continue;

        if (sscanf (line, " NEWTON %d", &ibd_useNewton) == 1) continue;
        if (sscanf (line, " SPEEDUP %d", &ibd_check10) == 1) continue;

        if (sscanf (line, " OUTPUT %s", output_prefix) == 1) continue;

        if (sscanf (line, " HAVEMAP %d", &MAP) == 1) continue;
		
  	    if (sscanf (line, " HAVEALLELEFREQ %d", &ALLELEFREQ) == 1) continue;

        if (sscanf (line, " MYEPS %lg", &MYEPS) == 1) continue;

        if (sscanf (line, " MC %d", &MC) == 1) continue;


	if (sscanf (line, " APED %d", &APED) == 1) continue;
	if (sscanf (line, " WPED %d", &WPED) == 1) continue;

	if (sscanf (line, " MITER %d", &MITER) == 1) continue;

	if (sscanf (line, " NREP %d", &NREP) == 1) continue;


	if (sscanf (line, " MLRT_PRECISION p/%lf", &ALPHA) == 1) {
	  ALPHA = 1 / ALPHA;
	  printf("Setting MLRT precision to %g*p\n",ALPHA);
	  continue;
	}

	if (sscanf (line, " MLRT_PRECISION p*%lf", &ALPHA) == 1) {
	  printf("Setting MLRT precision to %g*p\n",ALPHA);
	  continue;
	}



        if (sscanf (line, " PTHRESH %lg", &PTHRESH) == 1) continue;

        if (sscanf (line, " MYEPS %lg", &MYEPS) == 1) continue;

        if (sscanf (line, " ERRORRATE %lg", &errorrate) == 1) continue;

        if (sscanf (line, " OPTION %d", &OPTION) == 1) continue;

        if (sscanf (line, " PEDIGREE %s", ped_name) == 1) continue;

        if (sscanf (line, " CHROMOSOME %s %s", chromnameidx[nochrom+1], chromnameped[nochrom+1]) == 2) {
            nochrom ++;
            continue;
        }

        if (sscanf (line, " INCLUDE_FOUNDERS %d", &include_founders) == 1) continue;
        if (sscanf (line, " INCLUDE_NONFOUNDERS %d", &include_nonfounders) == 1) continue;

        if (sscanf (line, " INCLUDE_CASES %d", &include_cases) == 1) continue;
        if (sscanf (line, " INCLUDE_CONTROLS %d", &include_controls) == 1) continue;
        if (sscanf (line, " INCLUDE_MISSING %d", &include_missing) == 1) continue;

        printf("Error: Invalid line in configuration file %d: '%s'\n",line_no, line);
        error = 1;
    }

    if(error == 1) { printf("\nPlease fix configuration file '%s' and rerun.\n", name); exit(1); }

    if ( OPTION != 1 && OPTION != 2 ) {
      printf("\nError: Incorrect option: OPTION %d\n\n", OPTION);
      printf("  use 1, if you only want to apply EIBD, AIBS, IBS\n");
      printf("  use 2, if you also want to apply MLRT\n\n");
      exit(0);
    }

    if(MAP == 0 && OPTION == 2) {
      printf("Warning: Cannot run MLRT when no MAP is provided.\n");
      OPTION = 1;
    }

    /*
    if(APED == 1 && WPED == 1) {
      printf("Error: Only one of APED and WPED has to be specified. Please check your configuartion file.\n");
      exit(0);
    }
    */

    file_exists(ped_name);
    file_exists(chromnameped[1]);

    return 1;
}

/*
 * String <-> integer map functions
 *
 * The following functions take care of the recoding of strings in inputs files
 * to integers. This makes lookup and compare much simpler and faster.
 *
 */


map<string,int> strmap;
map<int, string> strmapR;

int strmap_len = 0;

void stringmap_init() { stringmap_lookup("0"); strmap["-9"] = 0; strmapR[-9] = "0"; }

int stringmap_lookup(const char *s)
{
  string m(s);

  if(strmap.find(m) == strmap.end())
  {
    strmap[m] = strmap_len;
    strmapR[strmap_len] = m;
    strmap_len++;
  }

  return strmap[m];
}

const char * stringmap_get(int a)
{
  return strmapR[a].c_str();
}

int stringmap_max() {
  int max = 3;

  for(map<string,int>::const_iterator it = strmap.begin(); it != strmap.end(); ++it)
     {
      int l = strlen(it->first.c_str());
      if(l > max) max = l;
     }

  /* printf("MAX=%d\n",max); */

  return max;
}

void stringmap_free() { }


// Pick suitable output file names and open them

void open_output_files()
{
  char head[1000];


  char WPorAP[10] = "WP";

  if (WPED == 1 && APED == 1)
    strcpy(WPorAP, "WP+AP");

  else if (APED == 1)
    strcpy(WPorAP, "AP");

  strcpy(head, output_prefix);
  
  sprintf(file_summary, "%s_summary", head);
  sprintf(file_results, "%s_results", head);
  sprintf(file_warnings, "%s_warnings", head);
  sprintf(file_errors, "%s_errors", head);

  if(cluster_n > 0) {
      char s[100];
      sprintf(s,".%d",cluster_id+1);
      strcat(file_summary,s);
      strcat(file_results,s);
      strcat(file_warnings,s);
      strcat(file_errors,s);
  }

  if ((fp_summary = fopen(file_summary, "w")) == 0) {
    printf("\nError: Unable to open output file: %s \n\n", file_summary);
    exit(0);
  }

  if ((fp_results = fopen(file_results, "w")) == 0) {
    printf("\nError: Unable to open output file: %s \n\n", file_results);
    exit(0);
  }

  /*if ((fp_warnings = fopen(file_warnings, "w")) == 0) {
    printf("\nError: nUnable to open output file: %s \n\n", file_warnings);
    exit(0);
  }*/

}

/*
 *
 * Read marker spacings from a specific map file:
 *  modifiles: theta, cenm, chromlength
 *
 */


int read_spacings_from(FILE *fp_idx, int chr)
{
  int i = chr, j; /* chromosome number */
  int numread;
  double dtemp[10];
  int bugexit = 0;

  theta[i][1] = cenm[i][1] = 0;
  for (j = 1; j<= (nomark[i]-1); j++) {
    numread = fscanf(fp_idx,"%lf", &theta[i][j+1]);

    if ( numread != 1 ) {
	fprintf(stdout, "Error of marker spacing in file %s:\n",
		chromnameidx[i]);
	fprintf(stdout, "%d markers, should be (%d - 1) spacings\n\n",
		nomark[i], nomark[i]);
	++bugexit;
    }
    if ( theta[i][j+1] > 1 ) {
	fprintf(stdout, "Error of marker spacing in file %s:\n",
		chromnameidx[i]);
	fprintf(stdout, "Spacing between adjacent markers should be in recombination fraction (theta), not in centi-Morgan\n\n");
	++bugexit;
    }
    if ( theta[i][j+1] < 0 ) {
	fprintf(stdout, "Error of marker spacing in file %s:\n",
		chromnameidx[i]);
	fprintf(stdout, " Spacing cannot be negative value\n\n");
	++bugexit;
    }
    if ( theta[i][j+1] == 0 )
	theta[i][j+1] = 0.0001;

    /* different map function */
    cenm[i][j+1] = mapthetatocm(theta[i][j+1], maptype);
    chromlength[i] += cenm[i][j+1];
  }

  /* n markers n-1 spacings, not n */
  numread = fscanf(fp_idx,"%lf", &dtemp[1]);
  if ( numread == 1 ) {
    fprintf(stdout, "Error of marker spacing in file %s:\n",
	      chromnameidx[i]);
    fprintf(stdout, "%d markers, should be (%d - 1) spacings\n\n",
	      nomark[i], nomark[i]);
    ++bugexit;
  }

  /* add 10 CM for later simulation */
  chromlength[0] += chromlength[i];
  chromlength[i] += 10;

  return bugexit;
}


/*
 * read the spacings only from the map files
 */
void read_spacings()
{
  int errors = 0;
  int i;

  theta = (double **)malloc((size_t) ((nochrom+1)*sizeof(double*)));
  cenm = (double **)malloc((size_t) ((nochrom+1)*sizeof(double*)));

  for (i = 1; i <= nochrom; i++) {
    FILE *idx = fopen(chromnameidx[i], "r");

    theta[i] = (double *)malloc((size_t) ((nomark[i]+1)*sizeof(double)));
    cenm[i] = (double *)malloc((size_t) ((nomark[i]+1)*sizeof(double)));

    errors += read_spacings_from(idx, i);

    if ( errors != 0 ) {
        printf("Error in marker spacings in file %s. Prest will now exit.\n", chromnameidx[i]);
        exit(1);
    }
  }
}




/*
 * read_map_file():
 *   Read all chromosome map files. Computes:
 *   q[i]
 *   qcum[i]
 *   theta[i]
 *   cenm[i]
 *   noalle[i]
 * nomark[i]
 * chromlength[i]
 *
 */

void read_map_files_old(int have_spacings)
{
  int i,j,k, numread;
  FILE *fp_idx;
  int itemp[10];
  int bugexit = 0;

  q = (double ***)malloc((size_t) ((nochrom+1)*sizeof(double**)));
  qcum = (double ***)malloc((size_t) ((nochrom+1)*sizeof(double**)));
  theta = (double **)malloc((size_t) ((nochrom+1)*sizeof(double*)));
  cenm = (double **)malloc((size_t) ((nochrom+1)*sizeof(double*)));
  noalle = (int **)malloc((size_t) ((nochrom+1)*sizeof(int*)));
  //nomark = (int *)malloc((size_t) ((nochrom+1)*sizeof(int)));
  //chromlength = (double *)malloc((size_t) ((nochrom+1)*sizeof(double)));

  for (i = 0; i <= nochrom; i++) chromlength[i] = 0;
  nomark[0] = 0;

  /* open and read .idx files : markers info. */

  for (i = 1; i <= nochrom; i++) {

    printf("Reading map file %s.\n", chromnameidx[i]);
    /* nomark: no. of markers in each chromosome */
    fp_idx = fopen(chromnameidx[i], "r");

    if(!fp_idx) {
        printf("*** Error: Cannot read map files ***\n");
        MAP = 0; ALLELEFREQ = 0; OPTION = 1;
        exit(1);
        return;
    }

    fgets_adj(&buf, &buf_len, fp_idx);

    numread = sscanf(buf,"%d", &nomark[i]);
    if (numread != 1) {
      printf("\nError reading num_mark from file %s\n %s\n",
	     chromnameidx[i], buf);
      exit(0);
    }
    nomark[0] += nomark[i];

    q[i] = (double **)malloc((size_t) ((nomark[i]+1)*sizeof(double*)));
    qcum[i] = (double **)malloc((size_t) ((nomark[i]+1)*sizeof(double*)));
    noalle[i] = (int *)malloc((size_t) ((nomark[i]+1)*sizeof(int)));


    /* noalle: no. of alleles of each marker */
    for (j = 1; j<= nomark[i]; j++) {

      fgets_adj(&buf, &buf_len, fp_idx);

      numread = sscanf(buf,"%d %d", &itemp[1], &noalle[i][j]);
      if (numread != 2) {
	printf("\nError in file %s:\n %s\n", chromnameidx[i], buf);
	exit(1);
      }

      q[i][j] = (double *)malloc((size_t)
				 ((noalle[i][j]+1)*sizeof(double)));
      qcum[i][j] = (double *)malloc((size_t)
				    ((noalle[i][j]+1)*sizeof(double)));
      q[i][j][0] = qcum[i][j][0] = 0;

      /* q: allele freq. of each allele */
      for (k = 1; k <= noalle[i][j]; k++) {
	numread = fscanf(fp_idx, "%lf", &q[i][j][k]);
	if (numread != 1) {
	  printf("\nFATAL ERROR: Error of allele freq. in %d_th marker in file %s:\n",
		 j, chromnameidx[i]);
	  exit(1);
	}
	if (q[i][j][k] <= 0 || q[i][j][k] > 1) {
	  printf("\nFATAL ERROR: Error of allele freq. in %d_th marker in file %s (<=0 or > 1):\n",
		 j, chromnameidx[i]);
	  exit(1);
	}
	qcum[i][j][k] = qcum[i][j][k-1] + q[i][j][k];
      }
      fgets_adj(&buf, &buf_len, fp_idx);

      if ( fabs(qcum[i][j][noalle[i][j]]-1) > 0.01 ) {
	printf("FATAL ERROR: in file %s: marker %d:  ", chromnameidx[i], j);
	printf("sum of allele freq = %f\n\n", qcum[i][j][noalle[i][j]]);
	++bugexit;
      }
      qcum[i][j][noalle[i][j]] = 1;
    }

    /* theta/cenm: marker spacing in r.f. or centi-Morgan,
       important: n markers, n-1 spacings,
       and assign 0 to the first marker */
    if(have_spacings == 1) {
      theta[i] = (double *)malloc((size_t) ((nomark[i]+1)*sizeof(double)));
      cenm[i] = (double *)malloc((size_t) ((nomark[i]+1)*sizeof(double)));
      bugexit += read_spacings_from(fp_idx, i);
    }

    fclose(fp_idx);
  }

  if ( bugexit != 0 ) {
    printf("\nFATAL ERROR - Errors found in marker data.");
    printf("\n  Check prest_errors file, and re-run the program\n\n");
    exit(0);
  }
}


/*
 * int include_individual(father_id, mother_id, affection status):
 *  returns 1 if the specific individual should be included in internal
 *  allele frequency computations.
 *
 */


int
include_individual(int fid, int mid, int aff)
{
  int founder = 0;

  if(fid == 0 && mid == 0) founder = 1;

  if(founder == 0 && include_nonfounders == 0) return 0;
  if(founder == 1 && include_founders == 0) return 0;

  if(aff == 0 && include_missing == 0) return 0;
  if(aff == 1 && include_controls == 0) return 0;
  if(aff == 2 && include_cases == 0) return 0;

  return 1;
}


/*
 * void compute_allele_frequencies(): Computes allele frequencies internally
 *
 */

void compute_allele_frequencies()
{
  int i, j, k;

  q = (double ***) malloc((size_t) ((nochrom + 1) * sizeof(double**)));
  qcum = (double ***) malloc((size_t) ((nochrom + 1) * sizeof(double**)));

  // theta = (double **) malloc((size_t) ((nochrom + 1) * sizeof(double*)));
  // cenm = (double **) malloc((size_t) ((nochrom + 1) * sizeof(double*)));
  noalle = (int **) malloc((size_t) ((nochrom + 1) * sizeof(int*)));
  //nomark = (int *) malloc((size_t) ((nochrom + 1) * sizeof(int)));
  //chromlength = (double *) malloc((size_t) ((nochrom + 1) * sizeof(double)));

  //for (i = 0; i <= nochrom; i++) { nomark[i] = 0; }
  //nomark[0] = 0;

  printf("Estimating allele frequencies from genotype data.\n");

  /* No MAP provided by user such as GWAS data:: compute afreqs internally  */
  /* Calculate allele frequencies internally base on founders. */

  for (i = 1; i <= nochrom; i++)

  {
    FILE *PED;
    fpos_t file_start;
    char str[8][256], *pch;
    int mark = 0;

    file_exists(chromnameped[1]);
    PED = fopen(chromnameped[1], "r");

    if(debug) printf("Chromosome %d: %d markers\n",i, nomark[i]);

    for (i = 1; i <= nochrom; i++) {
      noalle[i] = new_int1d(nomark[i]+1);
      for(j = 0 ; j <= nomark[i] ; j++)  noalle[i][j] = 0;

      q[i] = (double **) malloc((size_t) ((nomark[i] + 1) * sizeof(double*)));
      qcum[i] = (double **) malloc((size_t) ((nomark[i] + 1) * sizeof(double*)));
    }

    /* Now compute number of alleles at each marker */
    
    while(1) {
      int n = fscanf(PED,"%s %s %s %s %s %s",str[0],str[1],str[2],str[3],str[4],str[5]);
      if(feof(PED)) break;

      if(n != 6) {
	printf("Error: PED file %s is not valid.\n", chromnameped[i]);
	exit(1);
      }

      for (i = 1; i <= nochrom; i++)
      for (j = 1; j <= nomark[i]; j++) {
	int a1, a2, n;

	n = fscanf(PED,"%d %d",&a1, &a2);
	if(a1 == -9) a1 = 0;
	if(a2 == -9) a2 = 0;

	if(n != 2) {
	  printf("Error: Cannot read all markers from file %s. Prest accepts only markers coded numerically.\n", chromnameped[i]);
	  exit(1);
	}

	if(a1 > noalle[i][j]) noalle[i][j] = a1;
	if(a2 > noalle[i][j]) noalle[i][j] = a2;
      }
    }


    for (i = 1; i <= nochrom; i++)
    for (j = 0 ; j <= nomark[i] ; j++) {

      /* If all data are missing then set number of alleles to 2 */
      if(noalle[i][j] == 0) noalle[i][j] = 2;
      if(noalle[i][j] == 1) noalle[i][j] = 2;

      /* Initialize q and qcum to zero */
      q[i][j] = new_double1d( noalle[i][j] + 1 );
      qcum[i][j] = new_double1d( noalle[i][j] + 1 );
      for(k = 0 ; k <= noalle[i][j]; k++) { q[i][j][k] = 0; qcum[i][j][k] = 0; }

    }

    /*
     *  Now compute q and qcum
     */
    rewind(PED);
    
    int totalused = 0;

    while(1) {
      int sex, aff, fid, mid, n;
      n = fscanf(PED,"%s %s %s %s %d %d",str[0],str[1],str[2],str[3], &sex, &aff);
      if(feof(PED)) break;

      fid = stringmap_lookup(str[2]);
      mid = stringmap_lookup(str[3]);

      /* Allow -9 for missing*/
      if(mid == -9) mid = 0;
      if(fid == -9) fid = 0;
      if(aff == -9) aff = 0;

      /* Check if individual should be included in sample */
      if( ! include_individual(fid, mid, aff) ) { fgets_adj(&buf, &buf_len, PED); continue; }

      totalused++;

      if(debug) printf("compute_allele_frequencies(): counting individual %s %s\n",str[0],str[1]);


      for (i = 1; i <= nochrom; i++)
      for(j = 1 ; j <= nomark[i] ; j++) {
	int a1,a2, n;

	n = fscanf(PED,"%d %d",&a1, &a2);

	if(a1 == -9) a1 = 0;
	if(a2 == -9) a2 = 0;
	if(a1 == 0 && a2 == 0) continue;

	q[i][j][a1] ++;
	q[i][j][a2] ++;
      }
    }


    /* Normalize frequencies */
    for (i = 1; i <= nochrom; i++) for(j = 1 ; j <= nomark[i] ; j++) {

      double sum = 0;

      for(k = 1 ; k <= noalle[i][j]; k++) sum += q[i][j][k];

      if(sum == 0) { 
          printf("Warning: No founders found for marker %d:%d\n", i, j);
          for(k = 1 ; k <= noalle[i][j]; k++) q[i][j][k] = 1/noalle[i][j];
          
          sum = 1;
      }

      // Give each allele at least 2% frequency
      //for(k = 1 ; k <= noalle[i][j]; k++) q[i][j][k] += (3*sum)/100;
      //sum = 0;
      //for(k = 1 ; k <= noalle[i][j]; k++) sum += q[i][j][k];

      for(k = 1 ; k <= noalle[i][j]; k++) {
	if(sum == 0) { q[i][j][k] = 0; continue; }
	q[i][j][k] = q[i][j][k] / sum;

        // Set a small frequency in case it is 0
        if(q[i][j][k] == 0) q[i][j][k] = 0.001/noalle[i][j];
      }

      // Renormalize
      sum = 0;
      for(k = 1 ; k <= noalle[i][j]; k++) sum += q[i][j][k];
      for(k = 1 ; k <= noalle[i][j]; k++) q[i][j][k] = q[i][j][k] / sum;

      if(setUniformFreq) {
        q[i][j][1] = 0.75;
        q[i][j][2] = 0.25;
        //for(k = 1 ; k <= noalle[i][j]; k++) q[i][j][k] = 1.0/(double)noalle[i][j];
      }
    }


    /* Compute qcum */
    for (i = 1; i <= nochrom; i++) for(j = 1 ; j <= nomark[i] ; j++) {
      q[i][j][0] = 0;
      qcum[i][j][0] = 0;
      for(k = 1 ; k <= noalle[i][j]; k++) qcum[i][j][k] = qcum[i][j][k-1] + q[i][j][k];
    }

    if(debug)
      for (i = 1; i <= nochrom; i++) for(j = 1 ; j <= nomark[i] ; j++) {
	printf("Marker %3d:%3d freq = (%2d) ",i,j,noalle[i][j]);
	for(k = 0 ; k <= noalle[i][j]; k++) {
	  printf("%.5f(%.5f) ",q[i][j][k],qcum[i][j][k]);
	}
	printf("\n");
      }

    /* Write computed allele frequencies to a map file */
    if (write_frequency_file == 1) {
      char fname[500];

      FILE *f;

      sprintf(fname, "freqs.map", i);
      f = fopen(fname, "w");

      for (i = 1; i <= nochrom; i++) for (j = 1; j <= nomark[i]; j++) {
        fprintf(f, "%3d %12s %10f %8d ",i, names[i][j], cenm[i][j], 0);
	fprintf(f, " %4d ", noalle[i][j]);

	for (k = 1; k <= noalle[i][j]; k++) {
	  fprintf(f, "%8.6f ", q[i][j][k]);
	}
	fprintf(f, "\n");

      }

      fclose(f);
    }

    printf("%d individuals used in allele frequency calculations.\n", (int)totalused);

    fclose(PED);
  }


}

void join_files(char *prefix, int skip_header)
{
    FILE *all = fopen(prefix,"w");

    for(int i = 1; i <= cluster_n; i++) {
        char s[100];
        sprintf(s, "%s.%d",prefix, i);

        FILE *f = fopen(s,"r");

        printf("Joining file %s.\n",s);

        if(!f) {
            printf("Error: Cannot open results file %s\n", s);
            exit(1);
        }

        int c = 0;
        while ( fgets_adj(&buf, &buf_len, f) ) {
            if(skip_header && i > 1 && c == 0) { c++; continue; }
            fputs(buf, all);
            c++;
        }


    }

}

void join_cluster_files()
{
    fclose(fp_summary);
    fclose(fp_results);
    // fclose(fp_warnings);

    cluster_n = cluster_join;
    
    //join_files(file_summary, 0);
    join_files(file_results, 1);
    //join_files(file_warnings, 0);
}



/* Allocate 1d arrays */
int* new_int1d(int x) { return (int *)malloc((size_t)  ((x+1)*sizeof(int))); }
double* new_double1d(int x) { return (double *)malloc((size_t)  ((x+1)*sizeof(double))); }
char* new_char1d(int x) { return (char *)malloc((size_t)  ((x+1)*sizeof(char))); }


/* Allocate 2d array of char */
char** new_char2d(int x,int y)
{
    int i;
    char **a =  (char **)malloc((size_t)  ((x+1)*sizeof(char*)));

    for (i = 0; i <= x; i++)
        a[i] = (char *)malloc((size_t) ((y+1)*sizeof(char)));

    return a;
}


/* Allocate 2d array */
double** new_double2d(int x,int y)
{
    int i;
    double **a =  (double **)malloc((size_t)  ((x+1)*sizeof(double*)));

    for (i = 0; i <= x; i++)
        a[i] = (double *)malloc((size_t) ((y+1)*sizeof(double)));

    return a;
}

/* Allocate 2d array of int */
int** new_int2d(int x,int y)
{
    int i;
    int **a =  (int **)malloc((size_t)  ((x+1)*sizeof(int*)));

    for (i = 0; i <= x; i++)
        a[i] = (int *)malloc((size_t) ((y+1)*sizeof(int)));

    return a;
}

geno_t** new_genot2d(int x,int y)
{
    int i;
    geno_t **a =  (geno_t **)malloc((size_t)  ((x+1)*sizeof(geno_t*)));

    for (i = 0; i <= x; i++)
        a[i] = (geno_t *)malloc((size_t) ((y+1)*sizeof(geno_t)));

    return a;
}
