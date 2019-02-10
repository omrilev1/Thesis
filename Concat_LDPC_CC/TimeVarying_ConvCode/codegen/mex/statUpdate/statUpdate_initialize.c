/*
 * statUpdate_initialize.c
 *
 * Code generation for function 'statUpdate_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "statUpdate.h"
#include "statUpdate_initialize.h"
#include "_coder_statUpdate_mex.h"
#include "statUpdate_data.h"

/* Function Definitions */
void statUpdate_initialize(void)
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

/* End of code generation (statUpdate_initialize.c) */
