#include "arrayUtility.h"
#include <stdlib.h>

#define max(a,b) a>b ? a:b
/* Read array of double from file */
int arrayRead(FILE *fid, double *v, size_t length) {
  size_t i;
  int rval=0;
  int rvalScanf;
  
  for (i=0; i<length; i++) {
    rvalScanf=fscanf(fid,"%le",&v[i]);
    if (rvalScanf==EOF || rvalScanf<1) {
      break;
    }
    rval+=rvalScanf;
  }
  return rval;
}

/*Write array of double to a file */
void arrayWrite(FILE *fid, double *v, size_t length, size_t stride) {
  size_t i;

  for (i=0; i<length; i++) {
    fprintf(fid," %+le", v[i]);
    if (stride>0 && i>0 && (i+1)%stride==0) {
      fprintf(fid, "\n");
    }
  }
  fprintf(fid,"\n");
}

/*Write array of double to screen */
void arrayWriteScreenStride(double *v, size_t length, size_t stride) {
  arrayWrite(stdout, v, length, stride);
}

/*Write array of double to screen without stride*/
void arrayWriteScreen(double *v, size_t length) {
  arrayWrite(stdout, v, length, 0);
}

/*Compute the l-infinity distance between two vectors*/
double arrayInftyDist(double *v1, double *v2, size_t length) {
  size_t i;
  double res=0;
  for (i=0; i<length; i++) {
    res=max(res,abs(v1[i]-v2[i]));
  }
  return res;
}
    
/*Set all the elements in the array to zero*/
void arrayZeros(double *v, size_t length) {
  size_t i;
  for (i=0; i<length; i++) {
    v[i]=0;
  }
}

/*Create a pointer to a vector of pointers to provide a 2-D matrix view of a vector */
/*The vector of pointers needs to be released with a free */
double **arrayAllocate2DView(double *A, size_t r, size_t c) {
  size_t i;
  double **matA;
  
  matA=malloc(sizeof(double *)*c);
  for (i=0; i<c; ++i) {
    matA[i]=&A[r*i];
  }
  
  return matA;
}
