#include <stdio.h>
#include <stdlib.h>

#include "arrayUtility.h"

#define L 2
int main(void) {
  double A[9*L];
  FILE *fid;
  char *fileName="testData/matrix1.txt";
  size_t length=9*L;
  double **matA;

  fid=fopen(fileName,"r");
  if (fid==NULL) {
    printf("Error: input test file file not found\n");
    return 1;
  }
  printf("Reading file %s\n", fileName);
  arrayRead(fid,A,length);
  printf("Displaying without stride\n");
  arrayWriteScreen(A,length);
  printf("Displaying with stride\n");
  arrayWriteScreenStride(A,length,9);

  matA=arrayAllocate2DView(A,9,L);
  printf("matA[0]=%p\n",matA[0]);
  printf("A[]=%p\n",A);
  printf("Element [0,0] of the View: %+e\n", matA[0][0]);
  printf("Element [0,1] of the View: %+e\n", matA[0][1]);
  printf("Element [1,0] of the View: %+e\n", matA[1][0]);
  free(matA);

  fclose(fid);
  return 0;
}
