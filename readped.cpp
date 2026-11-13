#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <string>

using namespace std;

// Functions to speed up reading of files by indexing individuals & file positions

// For each individual and genotype file, this contains the file position of the genotypes
map<string,fpos_t> fposmap[30];

// Save position of each individual in genotype file for quick reading

int indexgenotypes()
{
  //for (int k = 1; k <= nochrom; k++) {
  int k = 1;
  {
    FILE *fp_ped = fopen(chromnameped[k], "r");
    if (!fp_ped) {
      printf("Error: cannot open file %s\n", chromnameped[k]);
      exit(1);
    }

    char str[10][1000];
    fpos_t pos;
    int line = 1, numread;

    while (!feof(fp_ped)) {
      fgetpos(fp_ped, &pos);

      numread = fscanf(fp_ped, "%s %s %s %s %s %s", str[1], str[2], str[3], str[4], str[5], str[6]);
      if (numread != 6) break;

      string s = (string) str[1] + "__" + (string) str[2];

      fposmap[k][s] = pos;
      fgets_adj(&buf, &buf_len, fp_ped);

      line++;
    }
  }
}

int findpedgeno(pedigree &ped, int* nomem, int* nomark, int** noalle, geno_t** geno)
{
  int linenum;
  int j, m, numread;
  int error = 0;
  char str[7][500];

  FILE *fp_ped;
  map<string, fpos_t>::iterator it;
  //   printf("XXX No chrom= %d\n",nochrom);

  fp_ped = fopen(chromnameped[1], "r");
  if (!fp_ped) {
    printf("Error: cannot open file %s\n", chromnameped[1]);
    exit(1);
  }

  /* for each individual */
  for (j = 1; j <= nomem[2]; j++) {

    /* id of the individual */
    geno[j][0] = ped[j][1];

    char *id, *family;
    family = (char*) s(nomem[1]);
    id = (char*) s(geno[j][0]);

    if (debug) printf("Looking up %s %s\n", family, id);
    fflush(stdout);

    string st = (string) family + "__" + (string) id;

    it = fposmap[1].find(st);

    if (it == fposmap[1].end()) {
      printf("Warning: filling zeros for missing individual: %s %s\n", s(nomem[1]), s(geno[j][0]));

      for (m = 1; m <= nomark[0]; m++) {
        geno[j][2 * m - 1] = 0;
        geno[j][2 * m] = 0;
      }

      continue;
    }

    fpos_t pos = fposmap[1][st];
    fsetpos(fp_ped, &pos);

    numread = fscanf(fp_ped, "%s %s %s %s %s %s", str[1], str[2], str[3], str[4], str[5], str[6]);
    if (numread != 6) {
      printf("Error reading genotype data in file %s for individual %s\n", chromnameped[1], id);
      exit(1);
    }

    //printf("individual %s at position %d : %s,%s,%s,%s\n", id, (int) (pos),str[1],str[2],str[3],str[4]);
    
    if (strcmp(str[1], family) != 0 || strcmp(str[2], id) != 0) {
      printf("Problem reading genotypes for %s:%s : indexing failed\n", family, id);
      exit(1);
    }


    int N = 1;
    /* for each chromosome */
    for (int k = 1; k <= nochrom; k++) for (m = 1; m <= nomark[k]; m++) {

        linenum = 1;
        //printf("individual %s at position %d ", id, (int) (pos));

        /* Main loop to read genotypes */
        int a1, a2, n;
        n = fscanf(fp_ped, "%d %d", &a1, &a2);

        if (n != 2) {
          printf("\nError reading genotype data for %s (pedigree %s) in file %s\n",
                  s(geno[j][0]), s(nomem[1]), chromnameped[k]);

          printf("Note: Prest only accepts markers coded numerically.\n");
          exit(1);
        }

        // Order alleles (APO)
        if(a1 > a2) { int temp = a2; a2 = a1; a1 = temp; }

        geno[j][2 * N - 1] = a1;
        geno[j][2 * N] = a2;

        if ((a1 > noalle[k][m] || a1 < 0) ||
                (a2 > noalle[k][m] || a2 < 0) ||
                (a1 == 0 && a2 != 0) ||
                (a1 != 0 && a2 == 0)) {
          fprintf(stdout, "Genotyping errors:\n");
          fprintf(stdout, "pedigree %s, at %s's  %dth marker\n",
                  s(nomem[1]), chromnameped[1], N);
          fprintf(stdout, "%s's genotype is %d %d\n",
                  s(ped[j][1]), geno[j][2 * N - 1], geno[j][2 * N]);
          fprintf(stdout, "no. of alleles of this marker is %d\n\n",
                  noalle[k][m]);
          ++error;
        }


        N++;
      }

  }

  fclose(fp_ped);

  return error;
}

