#include "mex.h"

/*
 * xtimesy.c - example found in API guide
 *
 * multiplies an input scalar times an input matrix and outputs a
 * matrix 
 *
 * This is a MEX-file for MATLAB.
 * Copyright 1984-2011 The MathWorks, Inc.
 */


void qtrapzc(double *x, int dim, double *z, mwSize n1, mwSize n2, mwSize n3)
{
    mwSize i,j,k;
    double sum;
    
    for (k=0; k<n3; k++) {
        for (j=0; j<n2; j++) {            
            sum = 0.;            
            sum += 0.5 * *(x++);            
            for (i=1; i<n1-1; i++) {
                sum += *(x++);
            }            
            sum += 0.5 * *(x++);            
            *(z++) = sum;
        }
    }
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  const mwSize *dims;
  double *x,*z;
  int d;
  size_t n1,n2,n3;
  
  /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgIdAndTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgIdAndTxt is executed. (mexErrMsgIdAndTxt breaks you out of
     the MEX-file) */
  if(nrhs!=2) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumInputs",
            "Two inputs required.");
  if(nlhs!=1) 
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:invalidNumOutputs",
            "One output required.");
  
  /* check to make sure the first input argument is a scalar */
  /* if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      mxGetN(prhs[0])*mxGetM(prhs[0])!=1 ) {
    mexErrMsgIdAndTxt( "MATLAB:xtimesy:xNotScalar",
            "Input x must be a scalar.");
  } */
  
  
  /*  create a pointer to the input matrix x */
  x = mxGetPr(prhs[0]);

  /*  get the dimension input */
  d = (int) mxGetScalar(prhs[1]);
  
  /*  get the dimensions of the matrix input y */
  dims = mxGetDimensions(prhs[0]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix((mwSize) dims[1], (mwSize) dims[2], mxREAL);
  
  /*  create a C pointer to a copy of the output matrix */
  z = mxGetPr(plhs[0]);
  
  /*  call the C subroutine */
  qtrapzc(x,d,z,(mwSize)dims[0],(mwSize)dims[1],(mwSize)dims[2]);
  
}
