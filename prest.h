


/*
 *  Function forward declerations
 */

void ibs_stats(geno_t **geno, int i, int index1, geno_t **genoB, int j, int index2,  int reltype);

/* Perform tests on a pair of individuals */
void test_pair(geno_t **geno, int i, int index1, geno_t **genoB, int j, int index2,  int reltype);

/* for find all the relationships considered */
int findrelation(pedigree &ped, int pedsize, int** rel);
int check_halfsib(int i, int j, int com, int ucom, pedigree &ped,
		  int pedsize);
int check_grandpc(int gp, int gc, int gcp, int gcup, pedigree &ped,
		  int pedsize);
int check_avuncular(int un, int ne, int np, int nup, pedigree &ped,
		    int pedsize);
int check_halfavuncular(int un, int ne, int np, int nup, pedigree &ped,
			int pedsize);
int check_cousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		 pedigree &ped, int pedsize);
int check_halfcousin(int i, int j, int c1p, int c2p, int c1up, int c2up,
		     pedigree &ped, int pedsize);
int check_halfpluscousin(int i, int j, int com, int ucom,
			 pedigree &ped, int pedsize);
void build_anclist(pedigree &ped, int pedsize, int id1,
		   int anclist[MANCE][3], int* totalnum, int flag);
int check(pedigree &ped, int pedsize, int id1, int id2);
int check_anclist(pedigree &ped, int pedsize, int id2,
		  int anclist[MANCE][3], int totalnum);
int check1_anclist(pedigree &ped, int pedsize, int id1);
int inbred_check(pedigree &ped, int pedsize, int id);
int recoding(pedigree &ped, int* nomem);

/* for find genodata, and Mendelian error checking */
int findpedgeno(pedigree &ped, int* nomem,
		int* nomark, int** noalle,
		geno_t** geno);
int checkMend(FILE *out, pedigree &ped, geno_t** geno, int* nomem,
	      int* nomark);

/* for trans. prob., cond. prob. and other prob.*/
double cp2(int f1, int m1, int f2, int m2, double* q);
double cp1(int f1, int m1, int f2, int m2, double* q);
double cp0(int f1, int m1, int f2, int m2, double* q);

void transprob(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1]);
void transproball(double dij, double pibdiibdj[NRMLRT+1][3+1][3+1]);

double pnorm(double z);

/* map functions */
double mapthetatocm(double theta, int maptype);
double mapcmtotheta(double cm, int maptype);

/* EM algorithm for estimation of p0, p1, p2 (0-1-2 IBD sharing) */
void   EMforp012(geno_t* geno1, geno_t* geno2, int* valid,
		 double*** q, int nochrom, int* nomark,
		 int reltype, double phat[3+1]);

/* EIBD, AIBS, IBS statistic, mu and std calculation */
void get3statandmustd(geno_t* geno1, geno_t* geno2, int* valid,
		      int commark, int reltype, int nochrom,
		      int* nomark, double p012[3+1],
		      double*** q, double** cenm, int maptype,
		      double*** eeibdibd, double*** eeibdeibdibd,
		      double*** eaibsibd, double*** eaibsaibsibd,
		      double** eibsibd, double** eibsibsibd,
		      double statresult[3+1][5+1], double pibs[3]);

void forEIBDvar(double*** eeibdeibd, double*** eeibdeibdibd,
		int nochrom, int* nomark, int* ngeno, int*** genotype,
		double*** q, double p012[NRMLRT+1][3+1]);

void forAIBSvar(double*** eaibsibd, double*** eaibsaibsibd,
		int nochrom, int* nomark, int* ngeno, int*** genotype,
		double*** q, double p012[NRMLRT+1][3+1]);

void forIBSvar(double** eibsibd, double** eibsibsibd,
	       int nochrom, int* nomark, double*** q, int** noalle);

void get3stat(geno_t* geno1, geno_t* geno2, int* valid, int commark,
	      int reltype, int nochrom, int* nomark,
	      double p012[4], double*** q,
	      double statresult[3+1], double pibs[3]);

/* for likelihood calculation, use Markov approximation if IBD
   process is MC */
void getapprloglikelihood(geno_t* geno1, geno_t* geno2, int* valid,
			  double*** q, double** cenm,
			  int nochrom, int* nomark, int maptype,
			  double loglhresult[NRMLRT+1]);

/* for simulate genodata of the relaitonships considered */
void simugeno1geno2(double*** q, double*** qcum, int nochrom,
		    int* nomark, int **noalle, int reltype,
		    double* chromlength, double** cenm,
		    geno_t* simugeno1, geno_t* simugeno2);
int crossover_new(double chromlength, double* pos, char c1, char c2,
                  char* ind);
int combine(double* p1, double* p2, double* ps, int nc1, int nc2,
	    int ncs, char* i1, char* i2, char* is, double* pnew,
	    char* inew);
void getibd2(double* p1,double* p2, double* pnew, int nc1,int nc2,
	     int ncnew, char* i1,char* i2, char* inew,
	     double chromlength, double* cenm, int nomark,
	     int* ibd);
void getibd1(double* pnew1, double* pnew2, int ncnew1, int ncnew2,
	     char* inew1, char* inew2, double chromlength,
	     double* cenm, int nomark, int* ibd);
void getibd0(double* pnew, int ncnew, char* inew, char c1, char c2,
	     double chromlength, double* cenm, int nomark,
	     int* ibd);


void open_output_files();
void read_map_files(int have_spacings);
void compute_allele_frequencies();
void read_spacings();

void join_cluster_files();

/* Auxiliarry functions */

char* fgets_adj(char **s, int *s_len, FILE *in);
void file_exists(char *name);
int read_config_file(char *name);

/* Helper functions to allocate 1d and 2d arrays */
int* new_int1d(int x);
double* new_double1d(int x);
char* new_char1d(int x);
double** new_double2d(int x,int y);
char** new_char2d(int x,int y);
int** new_int2d(int x,int y);
geno_t** new_genot2d(int x,int y);


/*
 * Functions to handle strings in input file.
 * stringmap_lookup(s): returns a unique integer for every string in input file
 * stringmap_get(a): get the string associated with integer a
 * s() an alias of stringmap_get
 */

void stringmap_init();
int stringmap_lookup(const char *s);
const char * stringmap_get(int a);
void stringmap_free();
int stringmap_max();
const char *s(int a) { return stringmap_get(a); }
