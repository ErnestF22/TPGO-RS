#include "arrayUtility.h"
#include <stdio.h>
/* Read array of double from file */
int arrayRead(FILE *fid, double *v, size_t length) {
  size_t i;
  int rval=0;
  int rvalScanf;
  
  for (i=0; i<length; i++) {
    rvalScanf=fscanf(fid,"%e",&v[i]);
    if (rvalScanf==EOF || rvalScanf<1) {
      break;
    }
    rval+=rvalScanf;
  }
  return rval;
}

/*Write array of double to a file */
void arrayWrite(FILE *fid, double *v, size_t length, size_t stride=0) {
  size_t i;

  for (i=0; i<length; i++) {
    fprintf(fid," %+e", v[i]);
    if (stride>0 && i>0 && (i+1)%stride==0) {
      fprintf(fid, "\n");
    }
  }
}

/*Write array of double to screen */
int arrayWriteScreen(FILE *fid, double *v, size_t length, size_t stride=0) {
  size_t i;

  arrayWrite(stdout, v, length, stride);
}


    
    
