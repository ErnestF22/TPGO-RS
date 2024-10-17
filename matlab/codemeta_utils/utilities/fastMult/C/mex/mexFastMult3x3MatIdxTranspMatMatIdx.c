#include "mex.h"
#include "mexCheckMacros.h"
#include "../arrayUtility.h"
#include "../fastMult.h"

/*
Similar to fastMult3x3MatIdxTranspMatMatIdx

function D=mexFastMult3x3MatIdxTranspMatMatIdx(A, eA, B, C, eC)
Inputs (5)
  [0] A          3 x 3 x N
  [1] eA         NEdges x 1, indeces for A
  [2] B          3 x 3 x NEdges
  [3] C          3 x 3 x N
  [4] eC         NEdges x 1, indeces for C
Outputs (1)
  [0] D          3 x 3 x NEdges
*/

/* Entry point */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  unsigned int i;
  double *A, *B, *C, *D, *eA, *eC;
  size_t NEdges;

  mexCheckNargin(5);
  mexCheckNargout(1);
    
  /* make sure the arguments are real */
  for (i=0; i<5; i++) {
    mexCheckArginReal(i);
  }

  /* check dimensions of the inputs */
  /* mexCheckArginDims(0,3,3); */
  mexCheckArginNCols(1,1);
  /* mexCheckArginDims(2,3,3); */
  /* mexCheckArginDims(3,3,3); */
  mexCheckArginNCols(4,1);
      
  /* check inputs have consistent length */
  mexCheckGeneric(9*mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]),\
		  "wrongNumEl", "Inputs 2 and 3 must have consistent dimensions.");
  mexCheckGeneric(9*mxGetNumberOfElements(prhs[4]) != mxGetNumberOfElements(prhs[2]),\
		  "wrongNumEl", "Inputs 3 and 5 must have consistent dimensions.");

  /* get sizes and prepare pointers*/
  NEdges=mxGetNumberOfElements(prhs[2])/9;

  A=mxGetPr(prhs[0]);
  eA=mxGetPr(prhs[1]);
  B=mxGetPr(prhs[2]);
  C=mxGetPr(prhs[3]);
  eC=mxGetPr(prhs[4]);

  /* create the output matrix */
  plhs[0] = mxCreateDoubleMatrix(9,(mwSize) NEdges, mxREAL);

  /* get a pointer to the real data in the output matrix */
  D = mxGetPr(plhs[0]);

  /* call the computational routine */
  fastMult3x3MatIdxTranspMatMatIdx_opt(A, eA, B, C, eC, D, NEdges);
}
