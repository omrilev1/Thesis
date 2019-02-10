/*
 * minSum_initialize.c
 *
 * Code generation for function 'minSum_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "minSum.h"
#include "minSum_initialize.h"
#include "_coder_minSum_mex.h"
#include "minSum_data.h"

/* Function Definitions */
void minSum_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (minSum_initialize.c) */
