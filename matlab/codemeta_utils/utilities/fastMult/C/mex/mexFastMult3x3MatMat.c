#include "mex.h"
#include "../fastMult.h"
#include "mexCheckMacros.h"

/* Entry point */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C;
    size_t numEl;                   /* size of matrix */

    mexCheckNargin(2);
    mexCheckNargout(1);
    
    /* make sure the arguments are real */
    mexCheckArginReal(0);
    mexCheckArginReal(1);
      
    /* check inputs have the same length */
    mexCheckGeneric(mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]),\
		    "wrongNumEl", "Inputs must have the same number of elements.");
    
    /* check inputs have a number of elements divisible by 9 */
    mexCheckNElemMultiple(0,9);
    mexCheckNElemMultiple(1,9);
    

    /* create a pointer to the real data in the input matrix  */
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);    

    /* get dimensions of the input matrix */
    numEl = mxGetNumberOfElements(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(9,(mwSize)numEl/9,mxREAL);

    /* get a pointer to the real data in the output matrix */
    C = mxGetPr(plhs[0]);

    /* call the computational routine */
    fastMult3x3MatMat_opt(A,B,C,numEl/9);
}
