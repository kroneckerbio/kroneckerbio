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
 * Auto code generator based originally on the timestwo.c file provided by the Mathworks
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
  double *t,*x,*u,*k;
  size_t trows,tcols,xrows,xcols,urows,ucols,krows,kcols;
  int xcolsmatch,ucolsmatch,kcolsmatch;
  double *prptr;
  mwIndex *irptr, *jcptr;
  /* Initialize sparse matrix indexing */
  #if NZMAX != 0
  double dsym[NZMAX];
  static const mwIndex ir[] = {%IRSTR%};
  static const mwIndex jc[] = {%JCSTR%};
  #endif
  
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
  /* Allow number of columns to be zero or one if the number of rows is supposed to be zero. Otherwise expect one column.*/
  if( XROWS==0  && xcols==0 ) {
      xcolsmatch = 0;
  }
  else {
      xcolsmatch = 1;
  };
    if( UROWS==0  && ucols==0 ) {
      ucolsmatch = 0;
  }
  else {
      ucolsmatch = 1;
  };
  if( KROWS==0  && kcols==0 ) {
      kcolsmatch = 0;
  }
  else {
      kcolsmatch = 1;
  };
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
      !(trows==1 && tcols==1) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
   if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
      !(xrows==XROWS && xcols==xcolsmatch) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
    if( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ||
      !(urows==UROWS && ucols==ucolsmatch) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
    if( !mxIsDouble(prhs[3]) || mxIsComplex(prhs[3]) ||
      !(krows==KROWS && kcols==kcolsmatch) ) {
    mexErrMsgIdAndTxt( "MATLAB:timestwo:inputNotRealScalarDouble",
            "Inputs must be noncomplex double column vectors with sizes 1, nx, nu, and nk, respectively.");
  }
 
  
  /* Assign pointers to each input and output. */
  t = mxGetPr(prhs[0]);
  x = mxGetPr(prhs[1]);
  u = mxGetPr(prhs[2]);
  k = mxGetPr(prhs[3]);
      
  /* Initialize an empty sparse matrix */
  plhs[0] = mxCreateSparse((mwSize)ROWS, (mwSize)COLS, (mwSize)NZMAX, mxREAL);
  
  /* Get pointers to relevant arrays of the output */
  prptr = mxGetPr(plhs[0]);
  irptr = mxGetIr(plhs[0]);
  jcptr = mxGetJc(plhs[0]);
  
  #if NZMAX != 0
  
  /* Call the %FUN% subroutine to calculate pr, called dsym. */
  %FUN%(dsym,t,x,u,k);
  
  /* Fill the matrix */
  memcpy((void*)prptr, (const void*)dsym, sizeof(dsym));
  memcpy((void*)irptr, (const void*)ir, sizeof(ir));
  memcpy((void*)jcptr, (const void*)jc, sizeof(jc));
  
  #endif
}
