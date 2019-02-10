/*
 * _coder_statUpdate_mex.c
 *
 * Code generation for function '_coder_statUpdate_mex'
 *
 */

/* Include files */
#include "statUpdate.h"
#include "_coder_statUpdate_mex.h"
#include "statUpdate_terminate.h"
#include "_coder_statUpdate_api.h"
#include "statUpdate_initialize.h"
#include "statUpdate_data.h"

/* Function Declarations */
static void statUpdate_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
  const mxArray *prhs[4]);

/* Function Definitions */
static void statUpdate_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T nrhs,
  const mxArray *prhs[4])
{
  const mxArray *outputs[2];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 4, 4,
                        10, "statUpdate");
  }

  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 10,
                        "statUpdate");
  }

  /* Call the function. */
  statUpdate_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  statUpdate_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(statUpdate_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  statUpdate_initialize();

  /* Dispatch the entry-point. */
  statUpdate_mexFunction(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_statUpdate_mex.c) */
