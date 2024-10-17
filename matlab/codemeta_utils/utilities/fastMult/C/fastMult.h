#ifndef __fastMult__
#define __fastMult__

#include <stdio.h>

void fastMult3x3MatMat(double *A, double *B, double *C, size_t length);
void fastMult3x3MatMat_opt(double *A, double *B, double *C, size_t length);
void fastMult3x3HatMat(double *v, double *A, double *B, size_t length);
void fastMult3x3MatHat(double *A, double *v, double *B, size_t length);
void fastMult3x3MatMatMat(double *A, double *B, double *C, double *D, size_t length);
void fastMult3x3MatMatMat_opt(double *A, double *B, double *C, double *D, size_t length);
void fastMult3x3HatNormSqMat(double *v, double *A, double *B, size_t length);
void fastMult3x3MatHatNormSq(double *A, double *v, double *B, size_t length);
void fastMult3x3MatIdxTranspMatMatIdx(double *A, double *eA, double *B, double *C, double *eC, double *D, size_t length);

#endif
