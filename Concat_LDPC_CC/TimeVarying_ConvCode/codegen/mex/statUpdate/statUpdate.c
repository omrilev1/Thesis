/*
 * statUpdate.c
 *
 * Code generation for function 'statUpdate'
 *
 */

/* Include files */
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "statUpdate.h"
#include "statUpdate_emxutil.h"
#include "statUpdate_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 8,     /* lineNo */
  "statUpdate",                        /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 11,  /* lineNo */
  "statUpdate",                        /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 25,  /* lineNo */
  "cat",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 100, /* lineNo */
  "cat",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 2,   /* lineNo */
  1,                                   /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 3, /* lineNo */
  1,                                   /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 112,/* lineNo */
  9,                                   /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 8, /* lineNo */
  44,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo e_emlrtRTEI = { 8, /* lineNo */
  5,                                   /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 11,/* lineNo */
  10,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo g_emlrtRTEI = { 11,/* lineNo */
  13,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 11,/* lineNo */
  22,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 11,/* lineNo */
  21,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo j_emlrtRTEI = { 1, /* lineNo */
  25,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo k_emlrtRTEI = { 103,/* lineNo */
  1,                                   /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

static emlrtRTEInfo m_emlrtRTEI = { 5, /* lineNo */
  9,                                   /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtDCInfo emlrtDCI = { 8,     /* lineNo */
  20,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 8,   /* lineNo */
  20,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 8,   /* lineNo */
  36,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = { 8,   /* lineNo */
  36,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  8,                                   /* lineNo */
  49,                                  /* colNo */
  "newL",                              /* aName */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  8,                                   /* lineNo */
  51,                                  /* colNo */
  "newL",                              /* aName */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  8,                                   /* lineNo */
  60,                                  /* colNo */
  "newL",                              /* aName */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = { 8,   /* lineNo */
  62,                                  /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  8,                                   /* lineNo */
  62,                                  /* colNo */
  "newL",                              /* aName */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { 1,   /* iFirst */
  1,                                   /* iLast */
  8,                                   /* lineNo */
  66,                                  /* colNo */
  "newL",                              /* aName */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { 1,   /* iFirst */
  1,                                   /* iLast */
  11,                                  /* lineNo */
  33,                                  /* colNo */
  "newL",                              /* aName */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m",/* pName */
  0                                    /* checkKind */
};

static emlrtECInfo emlrtECI = { -1,    /* nDims */
  11,                                  /* lineNo */
  5,                                   /* colNo */
  "statUpdate",                        /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 281,/* lineNo */
  27,                                  /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

/* Function Definitions */
void statUpdate(const emlrtStack *sp, real_T numOfChuncks, real_T Ms, const
                emxArray_real_T *newL, const emxArray_real_T *Lqij_init,
                emxArray_real_T *Lqij, emxArray_real_T *nextL)
{
  int32_T i0;
  int32_T loop_ub;
  int32_T k;
  emxArray_int32_T *r0;
  emxArray_int32_T *r1;
  emxArray_real_T *r2;
  cell_wrap_0 reshapes[2];
  emxArray_real_T *b_newL;
  emxArray_real_T *c_newL;
  int32_T b_loop_ub;
  real_T d0;
  int32_T result;
  int32_T i1;
  boolean_T empty_non_axis_sizes;
  int32_T b_result;
  int32_T newL_idx_0;
  int32_T i2;
  int32_T c_result;
  int32_T iv0[2];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  i0 = Lqij->size[0] * Lqij->size[1];
  Lqij->size[0] = Lqij_init->size[0];
  Lqij->size[1] = Lqij_init->size[1];
  emxEnsureCapacity_real_T(sp, Lqij, i0, &emlrtRTEI);
  loop_ub = Lqij_init->size[0] * Lqij_init->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    Lqij->data[i0] = Lqij_init->data[i0];
  }

  i0 = nextL->size[0] * nextL->size[1];
  nextL->size[0] = 0;
  nextL->size[1] = 0;
  emxEnsureCapacity_real_T(sp, nextL, i0, &b_emlrtRTEI);
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, numOfChuncks - 1.0, mxDOUBLE_CLASS,
    (int32_T)(numOfChuncks - 1.0), &m_emlrtRTEI, sp);
  k = 0;
  emxInit_int32_T(sp, &r0, 1, &j_emlrtRTEI, true);
  emxInit_int32_T(sp, &r1, 1, &j_emlrtRTEI, true);
  emxInit_real_T(sp, &r2, 2, &j_emlrtRTEI, true);
  emxInitMatrix_cell_wrap_0(sp, reshapes, &k_emlrtRTEI, true);
  emxInit_real_T(sp, &b_newL, 2, &d_emlrtRTEI, true);
  emxInit_real_T(sp, &c_newL, 2, &h_emlrtRTEI, true);
  while (k <= (int32_T)(numOfChuncks - 1.0) - 1) {
    /*  Future statistics section */
    if (1 > newL->size[0] - 2) {
      loop_ub = 0;
    } else {
      i0 = newL->size[0];
      if (!(1 <= i0)) {
        emlrtDynamicBoundsCheckR2012b(1, 1, i0, &emlrtBCI, sp);
      }

      i0 = newL->size[0];
      loop_ub = newL->size[0] - 2;
      if (!((loop_ub >= 1) && (loop_ub <= i0))) {
        emlrtDynamicBoundsCheckR2012b(loop_ub, 1, i0, &b_emlrtBCI, sp);
      }
    }

    if (1.0 > Ms) {
      b_loop_ub = 0;
    } else {
      i0 = newL->size[1];
      if (!(1 <= i0)) {
        emlrtDynamicBoundsCheckR2012b(1, 1, i0, &c_emlrtBCI, sp);
      }

      if (Ms != (int32_T)muDoubleScalarFloor(Ms)) {
        emlrtIntegerCheckR2012b(Ms, &e_emlrtDCI, sp);
      }

      i0 = newL->size[1];
      b_loop_ub = (int32_T)Ms;
      if (!((b_loop_ub >= 1) && (b_loop_ub <= i0))) {
        emlrtDynamicBoundsCheckR2012b(b_loop_ub, 1, i0, &d_emlrtBCI, sp);
      }
    }

    i0 = (int32_T)((1.0 + (real_T)k) + 1.0);
    if (!((i0 >= 1) && (i0 <= 1))) {
      emlrtDynamicBoundsCheckR2012b(i0, 1, 1, &e_emlrtBCI, sp);
    }

    st.site = &emlrtRSI;
    d0 = 2.0 * (Ms + 1.0) - 2.0;
    if (!(d0 >= 0.0)) {
      emlrtNonNegativeCheckR2012b(d0, &emlrtDCI, &st);
    }

    if (d0 != (int32_T)muDoubleScalarFloor(d0)) {
      emlrtIntegerCheckR2012b(d0, &b_emlrtDCI, &st);
    }

    if (!(Ms + 1.0 >= 0.0)) {
      emlrtNonNegativeCheckR2012b(Ms + 1.0, &c_emlrtDCI, &st);
    }

    if (Ms + 1.0 != (int32_T)muDoubleScalarFloor(Ms + 1.0)) {
      emlrtIntegerCheckR2012b(Ms + 1.0, &d_emlrtDCI, &st);
    }

    b_st.site = &c_emlrtRSI;
    if (!(((int32_T)(2.0 * (Ms + 1.0) - 2.0) == 0) || ((int32_T)(Ms + 1.0) == 0)))
    {
      result = (int32_T)(2.0 * (Ms + 1.0) - 2.0);
    } else if (!((loop_ub == 0) || (b_loop_ub == 0))) {
      result = loop_ub;
    } else {
      i1 = (int32_T)(2.0 * (Ms + 1.0) - 2.0);
      result = muIntScalarMax_sint32(i1, 0);
      if (loop_ub > result) {
        result = loop_ub;
      }
    }

    c_st.site = &d_emlrtRSI;
    if (((int32_T)(2.0 * (Ms + 1.0) - 2.0) == result) || (((int32_T)(2.0 * (Ms +
            1.0) - 2.0) == 0) || ((int32_T)(Ms + 1.0) == 0))) {
      empty_non_axis_sizes = true;
    } else {
      empty_non_axis_sizes = false;
      emlrtErrorWithMessageIdR2018a(&c_st, &n_emlrtRTEI,
        "MATLAB:catenate:matrixDimensionMismatch",
        "MATLAB:catenate:matrixDimensionMismatch", 0);
    }

    if ((loop_ub == result) || ((loop_ub == 0) || (b_loop_ub == 0))) {
    } else {
      empty_non_axis_sizes = false;
    }

    if (!empty_non_axis_sizes) {
      emlrtErrorWithMessageIdR2018a(&c_st, &n_emlrtRTEI,
        "MATLAB:catenate:matrixDimensionMismatch",
        "MATLAB:catenate:matrixDimensionMismatch", 0);
    }

    empty_non_axis_sizes = (result == 0);
    if (empty_non_axis_sizes || (!(((int32_T)(2.0 * (Ms + 1.0) - 2.0) == 0) ||
          ((int32_T)(Ms + 1.0) == 0)))) {
      b_result = (int32_T)(Ms + 1.0);
    } else {
      b_result = 0;
    }

    i0 = reshapes[0].f1->size[0] * reshapes[0].f1->size[1];
    reshapes[0].f1->size[0] = result;
    reshapes[0].f1->size[1] = b_result;
    emxEnsureCapacity_real_T(&b_st, reshapes[0].f1, i0, &c_emlrtRTEI);
    b_result *= result;
    for (i0 = 0; i0 < b_result; i0++) {
      reshapes[0].f1->data[i0] = 0.0;
    }

    if (empty_non_axis_sizes || (!((loop_ub == 0) || (b_loop_ub == 0)))) {
      b_result = b_loop_ub;
    } else {
      b_result = 0;
    }

    newL_idx_0 = newL->size[0];
    i0 = b_newL->size[0] * b_newL->size[1];
    b_newL->size[0] = loop_ub;
    b_newL->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(&b_st, b_newL, i0, &d_emlrtRTEI);
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_newL->data[i2 + b_newL->size[0] * i0] = newL->data[i2 + newL_idx_0 *
          i0];
      }
    }

    i0 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
    reshapes[1].f1->size[0] = result;
    reshapes[1].f1->size[1] = b_result;
    emxEnsureCapacity_real_T(&b_st, reshapes[1].f1, i0, &c_emlrtRTEI);
    for (i0 = 0; i0 < b_result; i0++) {
      for (i2 = 0; i2 < result; i2++) {
        reshapes[1].f1->data[i2 + reshapes[1].f1->size[0] * i0] = b_newL->
          data[i2 + result * i0];
      }
    }

    i0 = nextL->size[0] * nextL->size[1];
    nextL->size[0] = reshapes[0].f1->size[0];
    nextL->size[1] = reshapes[0].f1->size[1] + reshapes[1].f1->size[1];
    emxEnsureCapacity_real_T(&b_st, nextL, i0, &e_emlrtRTEI);
    loop_ub = reshapes[0].f1->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_loop_ub = reshapes[0].f1->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        nextL->data[i2 + nextL->size[0] * i0] = reshapes[0].f1->data[i2 +
          reshapes[0].f1->size[0] * i0];
      }
    }

    loop_ub = reshapes[1].f1->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_loop_ub = reshapes[1].f1->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        nextL->data[i2 + nextL->size[0] * (i0 + reshapes[0].f1->size[1])] =
          reshapes[1].f1->data[i2 + reshapes[1].f1->size[0] * i0];
      }
    }

    /*  Current statistics */
    i0 = k + 1;
    if (!((i0 >= 1) && (i0 <= 1))) {
      emlrtDynamicBoundsCheckR2012b(i0, 1, 1, &f_emlrtBCI, sp);
    }

    loop_ub = Lqij->size[0];
    i0 = r0->size[0];
    r0->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(sp, r0, i0, &f_emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      r0->data[i0] = i0;
    }

    loop_ub = Lqij->size[1];
    i0 = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity_int32_T(sp, r1, i0, &g_emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      r1->data[i0] = i0;
    }

    st.site = &b_emlrtRSI;
    b_st.site = &c_emlrtRSI;
    i0 = newL->size[0];
    i2 = newL->size[1];
    if (!((i0 == 0) || (i2 == 0))) {
      result = newL->size[1];
    } else if (!((nextL->size[0] == 0) || (nextL->size[1] == 0))) {
      result = nextL->size[1];
    } else {
      i0 = newL->size[1];
      if (i0 > 0) {
        result = newL->size[1];
      } else {
        result = 0;
      }

      if (nextL->size[1] > result) {
        result = nextL->size[1];
      }
    }

    c_st.site = &d_emlrtRSI;
    i0 = newL->size[1];
    if (i0 == result) {
      empty_non_axis_sizes = true;
    } else {
      i0 = newL->size[0];
      i2 = newL->size[1];
      if ((i0 == 0) || (i2 == 0)) {
        empty_non_axis_sizes = true;
      } else {
        empty_non_axis_sizes = false;
        emlrtErrorWithMessageIdR2018a(&c_st, &n_emlrtRTEI,
          "MATLAB:catenate:matrixDimensionMismatch",
          "MATLAB:catenate:matrixDimensionMismatch", 0);
      }
    }

    if ((nextL->size[1] == result) || ((nextL->size[0] == 0) || (nextL->size[1] ==
          0))) {
    } else {
      empty_non_axis_sizes = false;
    }

    if (!empty_non_axis_sizes) {
      emlrtErrorWithMessageIdR2018a(&c_st, &n_emlrtRTEI,
        "MATLAB:catenate:matrixDimensionMismatch",
        "MATLAB:catenate:matrixDimensionMismatch", 0);
    }

    empty_non_axis_sizes = (result == 0);
    if (empty_non_axis_sizes) {
      b_result = newL->size[0];
    } else {
      i0 = newL->size[0];
      i2 = newL->size[1];
      if (!((i0 == 0) || (i2 == 0))) {
        b_result = newL->size[0];
      } else {
        b_result = 0;
      }
    }

    if (empty_non_axis_sizes || (!((nextL->size[0] == 0) || (nextL->size[1] == 0))))
    {
      c_result = nextL->size[0];
    } else {
      c_result = 0;
    }

    newL_idx_0 = newL->size[0];
    loop_ub = newL->size[0];
    b_loop_ub = newL->size[1];
    i0 = c_newL->size[0] * c_newL->size[1];
    c_newL->size[0] = loop_ub;
    c_newL->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(&b_st, c_newL, i0, &h_emlrtRTEI);
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      for (i2 = 0; i2 < loop_ub; i2++) {
        c_newL->data[i2 + c_newL->size[0] * i0] = newL->data[i2 + newL_idx_0 *
          i0];
      }
    }

    i0 = r2->size[0] * r2->size[1];
    r2->size[0] = b_result + c_result;
    r2->size[1] = result;
    emxEnsureCapacity_real_T(&b_st, r2, i0, &i_emlrtRTEI);
    for (i0 = 0; i0 < result; i0++) {
      for (i2 = 0; i2 < b_result; i2++) {
        r2->data[i2 + r2->size[0] * i0] = c_newL->data[i2 + b_result * i0];
      }
    }

    for (i0 = 0; i0 < result; i0++) {
      for (i2 = 0; i2 < c_result; i2++) {
        r2->data[(i2 + b_result) + r2->size[0] * i0] = nextL->data[i2 + c_result
          * i0];
      }
    }

    iv0[0] = r0->size[0];
    iv0[1] = r1->size[0];
    emlrtSubAssignSizeCheckR2012b(&iv0[0], 2, &(*(int32_T (*)[2])r2->size)[0], 2,
      &emlrtECI, sp);
    b_result = Lqij->size[0];
    loop_ub = r2->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_loop_ub = r2->size[0];
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        Lqij->data[r0->data[i2] + b_result * r1->data[i0]] = r2->data[i2 +
          r2->size[0] * i0];
      }
    }

    k++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(sp, &c_newL);
  emxFree_real_T(sp, &b_newL);
  emxFreeMatrix_cell_wrap_0(sp, reshapes);
  emxFree_real_T(sp, &r2);
  emxFree_int32_T(sp, &r1);
  emxFree_int32_T(sp, &r0);

  /*  for k */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (statUpdate.c) */
