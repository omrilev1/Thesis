/*
 * _coder_minSum_mex.c
 *
 * Code generation for function '_coder_minSum_mex'
 *
 */

/* Include files */
#include "minSum.h"
#include "_coder_minSum_mex.h"
#include "minSum_terminate.h"
#include "_coder_minSum_api.h"
#include "minSum_initialize.h"
#include "minSum_data.h"

/* Function Declarations */
static void minSum_mexFunction(int32_T nlhs, mxArray *plhs[3], int32_T nrhs,
  const mxArray *prhs[5]);

/* Function Definitions */
static void minSum_mexFunction(int32_T nlhs, mxArray *plhs[3], int32_T nrhs,
  const mxArray *prhs[5])
{
  const mxArray *outputs[3];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 5) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 5, 4, 6,
                        "minSum");
  }

  if (nlhs > 3) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 6,
                        "minSum");
  }

  /* Call the function. */
  minSum_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  minSum_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(minSum_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  minSum_initialize();

  /* Dispatch the entry-point. */
  minSum_mexFunction(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_minSum_mex.c) */
