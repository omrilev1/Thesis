/*
 * _coder_minSum_api.c
 *
 * Code generation for function '_coder_minSum_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "minSum.h"
#include "_coder_minSum_api.h"
#include "minSum_emxutil.h"
#include "minSum_data.h"

/* Variable Definitions */
static emlrtRTEInfo u_emlrtRTEI = { 1, /* lineNo */
  1,                                   /* colNo */
  "_coder_minSum_api",                 /* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *prevL, const
  char_T *identifier, emxArray_real_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Ms, const
  char_T *identifier);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *Lci, const
  char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);

/* Function Definitions */
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  g_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m1;
  static const int32_T iv1[1] = { 0 };

  y = NULL;
  m1 = emlrtCreateNumericArray(1, iv1, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m1, (void *)&u->data[0]);
  emlrtSetDimensions((mxArray *)m1, u->size, 1);
  emlrtAssign(&y, m1);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *prevL, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(prevL), &thisId, y);
  emlrtDestroyArray(&prevL);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  h_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *Ms, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(Ms), &thisId);
  emlrtDestroyArray(&Ms);
  return y;
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *Lci, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  b_emlrt_marshallIn(sp, emlrtAlias(Lci), &thisId, y);
  emlrtDestroyArray(&Lci);
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv0[2] = { 0, 0 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m0, (void *)&u->data[0]);
  emlrtSetDimensions((mxArray *)m0, u->size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = i_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[1] = { -1 };

  const boolean_T bv0[1] = { true };

  int32_T iv2[1];
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims, &bv0[0],
    iv2);
  ret->size[0] = iv2[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv1[2] = { true, true };

  int32_T iv3[2];
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv1[0],
    iv3);
  ret->size[0] = iv3[0];
  ret->size[1] = iv3[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

static real_T i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

void minSum_api(const mxArray * const prhs[5], int32_T nlhs, const mxArray *
                plhs[3])
{
  emxArray_real_T *Lci;
  emxArray_real_T *prevL;
  emxArray_real_T *Lqij;
  emxArray_real_T *HT;
  emxArray_real_T *currL;
  emxArray_real_T *newLqij;
  emxArray_real_T *LQi;
  const mxArray *prhs_copy_idx_2;
  real_T Ms;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T1(&st, &Lci, 1, &u_emlrtRTEI, true);
  emxInit_real_T(&st, &prevL, 2, &u_emlrtRTEI, true);
  emxInit_real_T(&st, &Lqij, 2, &u_emlrtRTEI, true);
  emxInit_real_T(&st, &HT, 2, &u_emlrtRTEI, true);
  emxInit_real_T(&st, &currL, 2, &u_emlrtRTEI, true);
  emxInit_real_T(&st, &newLqij, 2, &u_emlrtRTEI, true);
  emxInit_real_T1(&st, &LQi, 1, &u_emlrtRTEI, true);
  prhs_copy_idx_2 = emlrtProtectR2012b(prhs[2], 2, false, -1);

  /* Marshall function inputs */
  emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "Lci", Lci);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "prevL", prevL);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs_copy_idx_2), "Lqij", Lqij);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "HT", HT);
  Ms = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "Ms");

  /* Invoke the target function */
  minSum(&st, Lci, prevL, Lqij, HT, Ms, currL, newLqij, LQi);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(currL);
  currL->canFreeData = false;
  emxFree_real_T(&st, &currL);
  HT->canFreeData = false;
  emxFree_real_T(&st, &HT);
  Lqij->canFreeData = false;
  emxFree_real_T(&st, &Lqij);
  prevL->canFreeData = false;
  emxFree_real_T(&st, &prevL);
  Lci->canFreeData = false;
  emxFree_real_T(&st, &Lci);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(newLqij);
  }

  newLqij->canFreeData = false;
  emxFree_real_T(&st, &newLqij);
  if (nlhs > 2) {
    plhs[2] = b_emlrt_marshallOut(LQi);
  }

  LQi->canFreeData = false;
  emxFree_real_T(&st, &LQi);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_minSum_api.c) */
