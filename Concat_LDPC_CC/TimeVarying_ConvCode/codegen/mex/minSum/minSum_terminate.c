/*
 * minSum_terminate.c
 *
 * Code generation for function 'minSum_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "minSum.h"
#include "minSum_terminate.h"
#include "_coder_minSum_mex.h"
#include "minSum_data.h"

/* Function Definitions */
void minSum_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void minSum_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (minSum_terminate.c) */
