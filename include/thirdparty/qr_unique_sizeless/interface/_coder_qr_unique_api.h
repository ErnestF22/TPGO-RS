//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// _coder_qr_unique_api.h
//
// Code generation for function 'qr_unique'
//

#ifndef _CODER_QR_UNIQUE_API_H
#define _CODER_QR_UNIQUE_API_H

// Include files
#include "coder_array_mex.h"
#include "emlrt.h"
#include "mex.h"
#include "tmwtypes.h"
#include <algorithm>
#include <cstring>

// Variable Declarations
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

// Function Declarations
void qr_unique(coder::array<real_T, 2U> *A, coder::array<real_T, 2U> *Q,
               coder::array<real_T, 2U> *R);

void qr_unique_api(const mxArray *prhs, int32_T nlhs, const mxArray *plhs[2]);

void qr_unique_atexit();

void qr_unique_initialize();

void qr_unique_terminate();

void qr_unique_xil_shutdown();

void qr_unique_xil_terminate();

#endif
// End of code generation (_coder_qr_unique_api.h)
