#include "mex.h"
#include "../fastMult.h"
#include "mexCheckMacros.h"

/* Entry point */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A, *B, *C, *D;
    size_t numElMat;                   /* size of matrix */

    mexCheckNargin(3);
    mexCheckNargout(1);
    
    /* make sure the arguments are real */
    mexCheckArginReal(0);
    mexCheckArginReal(1);
    mexCheckArginReal(2);
      
    /* check inputs have compatible lengths */
    mexCheckGeneric(mxGetNumberOfElements(prhs[0]) != mxGetNumberOfElements(prhs[1]),\
		    "wrongNumEl", "Inputs must have the same number of elements.");
    mexCheckGeneric(mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2]),\
		    "wrongNumEl", "Inputs must have the same number of elements.");
    
    /* check inputs have a number of elements divisible by 3 and 9 */
    mexCheckNElemMultiple(0,9);
    mexCheckNElemMultiple(1,9);
    mexCheckNElemMultiple(2,9);
    
    /* create a pointer to the real data in the input matrix  */
    A = mxGetPr(prhs[0]);    
    B = mxGetPr(prhs[1]);
    C = mxGetPr(prhs[2]);
    
    /* get dimensions of the input matrix */
    numElMat = mxGetNumberOfElements(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(9,(mwSize)numElMat/9,mxREAL);

    /* get a pointer to the real data in the output matrix */
    D = mxGetPr(plhs[0]);

    /* call the computational routine */
    fastMult3x3MatMatMat_opt(A,B,C,D,numElMat/9);
}
