/*
 * sign.c
 *
 * Code generation for function 'sign'
 *
 */

/* Include files */
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "minSum.h"
#include "sign.h"
#include "eml_int_forloop_overflow_check.h"
#include "minSum_data.h"

/* Variable Definitions */
static emlrtRSInfo l_emlrtRSI = { 13,  /* lineNo */
  "sign",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elfun\\sign.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 31,  /* lineNo */
  "applyScalarFunctionInPlace",        /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunctionInPlace.m"/* pathName */
};

/* Function Definitions */
void b_sign(const emlrtStack *sp, emxArray_real_T *x)
{
  int32_T nx;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &l_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  nx = x->size[0] * x->size[1];
  b_st.site = &m_emlrtRSI;
  if ((!(1 > nx)) && (nx > 2147483646)) {
    c_st.site = &n_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (k = 0; k < nx; k++) {
    x->data[k] = muDoubleScalarSign(x->data[k]);
  }
}

/* End of code generation (sign.c) */
