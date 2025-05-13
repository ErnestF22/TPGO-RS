#ifndef WAVELET_H
#define WAVELET_H

#include <cmath>
#include <iostream>
#ifndef MATLAB_MEX_FILE
#include <thirdparty/roptlib/cwrapper/blas/f2c.h>
#include "others_def.h"
#else
#include "blas.h"
#include "lapack.h"
#include "others_def.h"
#define integer ptrdiff_t
#endif

void haarFWT_1d(int n, realdpcomplex *v);
void haarFWT_1d_inverse(int n, realdpcomplex *v);

void haarFWT_2d(int n1, int n2, realdpcomplex *v);
void haarFWT_2d_inverse(int n, int n2, realdpcomplex *v);

#endif /* WAVELET_H */
