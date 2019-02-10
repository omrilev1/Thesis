/*
 * statUpdate_terminate.c
 *
 * Code generation for function 'statUpdate_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "statUpdate.h"
#include "statUpdate_terminate.h"
#include "_coder_statUpdate_mex.h"
#include "statUpdate_data.h"

/* Function Definitions */
void statUpdate_atexit(void)
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

void statUpdate_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (statUpdate_terminate.c) */
