

// Map function if there is no genetic map available

double
bp_to_cm(int chr, double bp)
{
  return bp / 7.5e5;
}


char *pos;

char *
read_line(FILE *fp)
{
  char *s = fgets_adj(&buf, &buf_len, fp);
  pos = buf;
  return s;
}

int
read_token(char *s)
{
  if (*pos == 0) return 0;
  while ((*pos != 0) && ((*pos) == ' ' || (*pos) == '\t')) pos++;
  if (*pos == 0) return 0;
  while ((*pos != 0) && ((*pos) != ' ' && (*pos) != '\t')) *s++ = *pos++;

  *s++ = 0;
  return 1;
}

extern int MAP;
// Return 1 if frequencies were present in file

int read_map_files()
{
  int i, j, k;
  char s[10][100];

  nomark[0] = 0;
  nochrom = 0;

  file_exists(chromnameidx[1]);
  FILE *fp = fopen(chromnameidx[1], "r");

  int line = 1;
  char error[100] = {0, 0};

  int have_freq = 0, have_cm = MAP;

  /* first pass: get nochrom and nomark */
  while (read_line(fp)) {

    int n = read_token(s[1]) + read_token(s[2]) + read_token(s[3]) + read_token(s[4]);

    if (n != 4) {
      strcpy(error, "Not enough items in map file.");
      break;
    }

    int chr = atoi(s[1]);
    double cm = atof(s[2]);

    if (chr == 0) {
      strcpy(error, "Chromosome has to be a number.");
      break;
    }

    if (strlen(s[2]) > max_mark) max_mark = strlen(s[2]);

    if (line == 1 && read_token(s[5])) have_freq = 1;

    if (chr > nochrom) nochrom = chr;
    nomark[chr]++;
    nomark[0]++;

    line++;
  }

  if (error[0] != 0) {
    printf("Error in file %s line %d: %s\n", chromnameidx[1], line, error);
    exit(1);
  }


  printf("%d chromosomes\n", nochrom);
  printf("%d total markers\n", nomark[0]);
  printf("Markers per chromosome: ");
  for (i = 1; i <= nochrom; i++) printf("%d ", nomark[i]);
  printf("\n\n");

  printf("Allele frequencies: ");
  if (have_freq) printf("specified in map file\n\n");
  else printf("estimated from genotypes\n\n");

  int n = nochrom;
  q = (double ***) malloc((size_t) ((n + 1) * sizeof (double**)));
  qcum = (double ***) malloc((size_t) ((n + 1) * sizeof (double**)));
  theta = (double **) malloc((size_t) ((n + 1) * sizeof (double*)));
  cenm = (double **) malloc((size_t) ((n + 1) * sizeof (double*)));
  noalle = (int **) malloc((size_t) ((n + 1) * sizeof (int*)));
  names = (char***) malloc((n + 1) * sizeof (char**));

  for (i = 1; i <= nochrom; i++) {
    if (have_freq) {
      q[i] = (double **) malloc((size_t) ((nomark[i] + 1) * sizeof (double*)));
      qcum[i] = (double **) malloc((size_t) ((nomark[i] + 1) * sizeof (double*)));
      noalle[i] = (int *) malloc((size_t) ((nomark[i] + 1) * sizeof (int)));
    }

    theta[i] = (double *) malloc((size_t) ((nomark[i] + 1) * sizeof (double)));
    cenm[i] = (double *) malloc((size_t) ((nomark[i] + 1) * sizeof (double)));
    cenm[i][0] = 0;
    cenm[i][1] = 0;

    names[i] = (char **) malloc((size_t) ((nomark[i] + 1) * sizeof (char*)));

    for (j = 0; j <= nomark[i]; j++) names[i][j] = (char *) malloc((max_mark + 2));
  }


  rewind(fp);
  line = 1;
  int chr_prev = -1;
  int markno = 0;

  double prev = 0;

  while (read_line(fp)) {
    char name[100];

    int n = read_token(s[1]) + read_token(s[2]) + read_token(s[3]) + read_token(s[4]);

    strcpy(name, s[2]);

    int chr = atoi(s[1]);
    double cm = atof(s[3]), bp = atof(s[4]);

    if (chr < chr_prev) {
      printf("Error: Markers should be sorted accorind to chr/bp\n");
      exit(1);
    }

    if (chr != chr_prev) markno = 1;

    strcpy(names[chr][markno], s[2]);

    //printf("%s %d %f %f\n",name, markno, cm, cm-prev);

    if (uniformSPC > 0) cm = bp_to_cm(chr, bp);

    {
      cenm[chr][markno] = cm - prev;
      prev = cm;
    }

    if (have_freq) {
      int n = read_token(s[1]);

      i = chr;
      j = markno;

      noalle[i][j] = atoi(s[1]);
      q[i][j] = (double *) malloc((size_t) ((noalle[i][j] + 1) * sizeof (double)));
      qcum[i][j] = (double *) malloc((size_t) ((noalle[i][j] + 1) * sizeof (double)));

      q[i][j][0] = qcum[i][j][0] = 0;

      for (k = 1; k <= noalle[i][j]; k++) {
        n = read_token(s[1]);

        q[i][j][k] = atof(s[1]);
        qcum[i][j][k] = qcum[i][j][k - 1] + q[i][j][k];

        if (n != 1 || q[i][j][k] <= 0 || q[i][j][k] > 1) {
          printf("Error reading allele frequency for marker %s in file %s.\n", name, chromnameidx[1]);
          exit(1);
        }
      }

      if (fabs(qcum[i][j][noalle[i][j]] - 1) > 0.01) {
        printf("Error in file %s: marker %s: ", chromnameidx[1], name);
        printf("sum of allele freq = %f\n\n", qcum[i][j][noalle[i][j]]);
        exit(1);
      }

      qcum[i][j][noalle[i][j]] = 1;

    }

    line++;
    markno++;
    chr_prev = chr;
  }

  /* Compute theta and chromlength */
  for (i = 0; i <= nochrom; i++) chromlength[i] = 0;

  for (int i = 1; i <= nochrom; i++) {

    chromlength[i] = 0;

    for (int j = 1; j <= nomark[i] - 1; j++) {
      chromlength[i] += cenm[i][j + 1];
      theta[i][j + 1] = mapcmtotheta(cenm[i][j + 1], maptype);
    }

    theta[i][1] = cenm[i][1] = 0;

    if (chromlength[i] < 0) chromlength[i] = 0;

    /* add 10 CM for later simulation */
    chromlength[0] += chromlength[i];
    chromlength[i] += 10;

    if (debug > 0) {
      printf("Chromosome %d genetic map:\n", i);
      for (int j = 1; j <= nomark[i]; j++) printf("%f, %f ", cenm[i][j], theta[i][j]);
      printf("\n");
    }
  }

  printf("Genetic map: ");
  if (uniformSPC > 0) printf("estimated as constant rate cm/bp.\n");
  else if (have_cm) printf("specified in map file.\n");
  else printf("no genetic map\n");

  if(have_cm || uniformSPC > 0) {
    printf("Total genetic map size: %f cm\n", chromlength[0]);
    printf("Genetic map size per chromosome: ");
  
    for (int i = 1; i <= nochrom; i++) printf("%.2f ", chromlength[i] - 10);
    printf("\n");
  }

  printf("\n");

  return have_freq;
}
