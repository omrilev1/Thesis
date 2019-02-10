/*
 * _coder_statUpdate_api.c
 *
 * Code generation for function '_coder_statUpdate_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "statUpdate.h"
#include "_coder_statUpdate_api.h"
#include "statUpdate_emxutil.h"
#include "statUpdate_data.h"

/* Variable Definitions */
static emlrtRTEInfo l_emlrtRTEI = { 1, /* lineNo */
  1,                                   /* colNo */
  "_coder_statUpdate_api",             /* fName */
  ""                                   /* pName */
};

/* Function Declarations */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *newL, const
  char_T *identifier, emxArray_real_T *y);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *numOfChuncks,
  const char_T *identifier);
static const mxArray *emlrt_marshallOut(const emxArray_real_T *u);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);

/* Function Definitions */
static real_T b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = e_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *newL, const
  char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  d_emlrt_marshallIn(sp, emlrtAlias(newL), &thisId, y);
  emlrtDestroyArray(&newL);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  f_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T emlrt_marshallIn(const emlrtStack *sp, const mxArray *numOfChuncks,
  const char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(numOfChuncks), &thisId);
  emlrtDestroyArray(&numOfChuncks);
  return y;
}

static const mxArray *emlrt_marshallOut(const emxArray_real_T *u)
{
  const mxArray *y;
  const mxArray *m0;
  static const int32_T iv1[2] = { 0, 0 };

  y = NULL;
  m0 = emlrtCreateNumericArray(2, iv1, mxDOUBLE_CLASS, mxREAL);
  emlrtMxSetData((mxArray *)m0, (void *)&u->data[0]);
  emlrtSetDimensions((mxArray *)m0, u->size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  static const int32_T dims[2] = { -1, -1 };

  const boolean_T bv0[2] = { true, true };

  int32_T iv2[2];
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims, &bv0[0],
    iv2);
  ret->size[0] = iv2[0];
  ret->size[1] = iv2[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)emlrtMxGetData(src);
  ret->canFreeData = false;
  emlrtDestroyArray(&src);
}

void statUpdate_api(const mxArray * const prhs[4], int32_T nlhs, const mxArray
                    *plhs[2])
{
  emxArray_real_T *newL;
  emxArray_real_T *Lqij_init;
  emxArray_real_T *Lqij;
  emxArray_real_T *nextL;
  real_T numOfChuncks;
  real_T Ms;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtHeapReferenceStackEnterFcnR2012b(&st);
  emxInit_real_T(&st, &newL, 2, &l_emlrtRTEI, true);
  emxInit_real_T(&st, &Lqij_init, 2, &l_emlrtRTEI, true);
  emxInit_real_T(&st, &Lqij, 2, &l_emlrtRTEI, true);
  emxInit_real_T(&st, &nextL, 2, &l_emlrtRTEI, true);

  /* Marshall function inputs */
  numOfChuncks = emlrt_marshallIn(&st, emlrtAliasP(prhs[0]), "numOfChuncks");
  Ms = emlrt_marshallIn(&st, emlrtAliasP(prhs[1]), "Ms");
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[2]), "newL", newL);
  c_emlrt_marshallIn(&st, emlrtAlias(prhs[3]), "Lqij_init", Lqij_init);

  /* Invoke the target function */
  statUpdate(&st, numOfChuncks, Ms, newL, Lqij_init, Lqij, nextL);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(Lqij);
  Lqij->canFreeData = false;
  emxFree_real_T(&st, &Lqij);
  Lqij_init->canFreeData = false;
  emxFree_real_T(&st, &Lqij_init);
  newL->canFreeData = false;
  emxFree_real_T(&st, &newL);
  if (nlhs > 1) {
    plhs[1] = emlrt_marshallOut(nextL);
  }

  nextL->canFreeData = false;
  emxFree_real_T(&st, &nextL);
  emlrtHeapReferenceStackLeaveFcnR2012b(&st);
}

/* End of code generation (_coder_statUpdate_api.c) */
