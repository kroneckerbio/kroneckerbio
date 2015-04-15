#include "mex.h"
#include <math.h>
#include "matrix.h"
#include <string.h>
#include <stdlib.h>

#define NZMAX %NZMAX%
#define ROWS %ROWS%
#define COLS %COLS%
#define XROWS %XROWS%
#define UROWS %UROWS%
#define KROWS %KROWS%

#if (defined(__GNUC__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
    #define inline __inline__
#elif defined(__clang__)
    #define inline inline
#elif defined(_MSC_VER)
    #define inline __inline
#else
    #define inline inline
#endif

/*
#ifndef __GNUC__
#define __inline__ __inline
#endif
 */

/*
 * Auto-generated derivative code for KroneckerBio
 * Auto code generator based originally on the timestwo.c file provided by Mathworks
 * Written by David Flowers
 * Tidor Lab, 2014
/* $Revision: 1.8.6.5 $ */

void %FUN%(double T[], double t[], double x[], double u[], double k[])
{
  /* Start modified ccode */
    %CSTR%
  /* End modified ccode */
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *t,*x,*u,*k,dsym[NZMAX];
  size_t trows,tcols,xrows,xcols,urows,ucols,krows,kcols;
  double *prptr;
  mwIndex *irptr, *jcptr;
  /* Initialize sparse matrix indexing */
  static const mwIndex ir[] = {%IRSTR%};
  static const mwIndex jc[] = {%JCSTR%};
  
  /* Check for proper number of arguments. */
  if(nrhs!=4) {
    mexErrMsgIdAndTxt( "MATLAB:%FUN%:invalidNumInputs",
            "Four inputs required.");
  } else if(nlhs>1) {
    mexErrMsgIdAndTxt( "MATLAB:%FUN%:maxlhs",
            "Too many output arguments.");
  }
  
  /* Check inputs' sizes and types.*/
  trows = mxGetM(prhs[0]);
  tcols = mxGetN(prhs[0]);
  xrows = mxGetM(prhs[1]);
  xcols = mxGetN(prhs[1]);
  urows = mxGetM(prhs[2]);
  ucols = mxGetN(prhs[2]);
  krows = mxGetM(prhs[3]);
  kcols = mxGetN(prhs[3]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(trows==1 && tcols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
   if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      !(xrows==XROWS && xcols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(urows==UROWS && ucols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      !(krows==KROWS && kcols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
 
  
  /* Assign pointers to each input and output. */
  t = mxGetPr(prhs[0]);
  x = mxGetPr(prhs[1]);
  u = mxGetPr(prhs[2]);
  k = mxGetPr(prhs[3]);
  
  /* Call the %FUN% subroutine to calculate pr, called dsym. */
  %FUN%(dsym,t,x,u,k);
      
  /* Initialize an empty sparse matrix */
  plhs[0] = mxCreateSparse((mwSize)ROWS, (mwSize)COLS, (mwSize)NZMAX, mxREAL);
  
  /* Get pointers to relevant arrays of the output */
  prptr = mxGetPr(plhs[0]);
  irptr = mxGetIr(plhs[0]);
  jcptr = mxGetJc(plhs[0]);
  
  /* Fill the matrix */
  memcpy((void*)prptr, (const void*)dsym, sizeof(dsym));
  memcpy((void*)irptr, (const void*)ir, sizeof(ir));
  memcpy((void*)jcptr, (const void*)jc, sizeof(jc));
}
