#ifndef __arrayUtility__
#define __arrayUtility__

#include <stdio.h>

int arrayRead(FILE *fid, double *v, size_t length);
void arrayWrite(FILE *fid, double *v, size_t length, size_t stride);
void arrayWriteScreenStride(double *v, size_t length, size_t stride);
void arrayWriteScreen(double *v, size_t length);
double arrayInftyDist(double *v1, double *v2, size_t length);
void arrayZeros(double *v, size_t length);
double **arrayAllocate2DView(double *A, size_t r, size_t c);

#endif
