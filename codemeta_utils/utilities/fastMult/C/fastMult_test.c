#include <stdio.h>
#include <stdlib.h>
#include "arrayUtility.h"
#include "fastMult.h"

#include <string.h>
#include <time.h>
#include <unistd.h>

#define DATADIR "../testData/"

#define measureTime(s) \
  start = clock(); \
  s; \
  end = clock(); \
  printf("\tclock(): S=%7lu E=%7lu D=%lu=%.3f ms\n", \
	 start,end,end-start,(end-start)/(double) CLOCKS_PER_SEC*1000.0);

#define runTest(testName, cmd) \
  printf(testName " Test\n"); \
  readFile(DATADIR testName ".txt", Rref, lengthMat); \
  printf("Exec"); \
  measureTime(cmd); \
  res=arrayInftyDist(R,Rref,lengthMat); \
  printf("\tL-infty distance between R and Rref: %le\n", res);
  /*arrayZeros(R,lengthMat);					\*/


void readFile(char *fileName, double *A, size_t length) {
  FILE *fidA;
  size_t readLength;

  printf("\t%s\n",fileName);
  fidA=fopen(fileName,"r");
  if (fidA==NULL) {
    printf("Error: file %s not found\n",fileName);
  }
  readLength=arrayRead(fidA,A,length);
  printf("\t\tRead %u doubles\n",readLength);
  fclose(fidA);
}

int main(void) {
  double *A, *B, *C, *v, *e1, *e2, *Rref, *R;
  FILE *fidA;
  unsigned long int L, L2, lengthMat, lengthVec, lengthMat2;
  double res;
  clock_t start, end;

  printf("Reading number of matrices and vectors\n");
  fidA=fopen(DATADIR "dim.txt","r");
  fscanf(fidA,"%lu",&L);
  fclose(fidA);
  printf("\tN=%lu\n", L);

  lengthMat=9*L;
  lengthVec=3*L;

  A=(double *) malloc(sizeof(double)*lengthMat);
  B=(double *) malloc(sizeof(double)*lengthMat);
  C=(double *) malloc(sizeof(double)*lengthMat);  
  v=(double *) malloc(sizeof(double)*lengthVec);
  R=(double *) malloc(sizeof(double)*lengthMat);  
  Rref=(double *) malloc(sizeof(double)*lengthMat);
  
  printf("Reading files with inputs\n");
  readFile(DATADIR "matrix1.txt", A, lengthMat);
  readFile(DATADIR "matrix2.txt", B, lengthMat);
  readFile(DATADIR "matrix3.txt", C, lengthMat);
  readFile(DATADIR "vector1.txt", v, lengthVec);

  runTest("MatMat",fastMult3x3MatMat(A,B,R,L));
  runTest("HatMat",fastMult3x3HatMat(v,A,R,L));
  runTest("MatHat",fastMult3x3MatHat(A,v,R,L));
  runTest("MatMatMat",fastMult3x3MatMatMat(A,B,C,R,L));
  runTest("MatMatMat",fastMult3x3MatMatMat_opt(A,B,C,R,L));
  runTest("HatNormSqMat",fastMult3x3HatNormSqMat(v,A,R,L));
  runTest("MatHatNormSq",fastMult3x3MatHatNormSq(A,v,R,L));
  
  free(A);
  free(B);
  free(C);
  free(v);

  e1=(double *) malloc(sizeof(double)*lengthVec);
  e2=(double *) malloc(sizeof(double)*lengthVec);

  printf("Reading second number of matrices and vectors\n");
  fidA=fopen(DATADIR "dim2.txt","r");
  fscanf(fidA,"%lu",&L2);
  fclose(fidA);
  printf("\tN2=%lu\n", L2);

  lengthMat2=9*L2;

  A=(double *) malloc(sizeof(double)*lengthMat2);
  B=(double *) malloc(sizeof(double)*lengthMat);
  C=(double *) malloc(sizeof(double)*lengthMat2);
  e1=(double *) malloc(sizeof(double)*L);
  e2=(double *) malloc(sizeof(double)*L);
  /* R and Rref have already been allocated before */

  printf("Reading files with inputs\n");
  readFile(DATADIR "matrixShort1.txt", A, lengthMat2);
  readFile(DATADIR "matrix2.txt", B, lengthMat);
  readFile(DATADIR "matrixShort2.txt", C, lengthMat2);
  readFile(DATADIR "idx1.txt", e1, L);
  readFile(DATADIR "idx2.txt", e2, L);

  runTest("MatIdxTranspMatMatIdx",fastMult3x3MatIdxTranspMatMatIdx_opt(A, e1, B, C, e2, R, L));

  free(A);
  free(B);
  free(C);
  free(e1);
  free(e2);
  free(R);
  free(Rref);
  
  return 0;
}
