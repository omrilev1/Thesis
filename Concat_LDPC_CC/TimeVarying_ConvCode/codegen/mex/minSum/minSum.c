/*
 * minSum.c
 *
 * Code generation for function 'minSum'
 *
 */

/* Include files */
#include "mwmathutil.h"
#include "rt_nonfinite.h"
#include "minSum.h"
#include "minSum_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "sum.h"
#include "sign.h"
#include "minSum_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 20,    /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 23,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 29,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 30,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 36,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 52,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 69,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 75,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 80,  /* lineNo */
  "minSum",                            /* fcnName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 25,  /* lineNo */
  "cat",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 100, /* lineNo */
  "cat",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 16,  /* lineNo */
  "abs",                               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elfun\\abs.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 74,  /* lineNo */
  "applyScalarFunction",               /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\applyScalarFunction.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 41,  /* lineNo */
  "find",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 153, /* lineNo */
  "find",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 377, /* lineNo */
  "find",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 397, /* lineNo */
  "find",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 18,  /* lineNo */
  "indexShapeCheck",                   /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\indexShapeCheck.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 10,  /* lineNo */
  "prod",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\datafun\\prod.m"/* pathName */
};

static emlrtRTEInfo emlrtRTEI = { 112, /* lineNo */
  9,                                   /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

static emlrtRTEInfo b_emlrtRTEI = { 20,/* lineNo */
  10,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo c_emlrtRTEI = { 20,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo d_emlrtRTEI = { 301,/* lineNo */
  14,                                  /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

static emlrtRTEInfo e_emlrtRTEI = { 23,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo f_emlrtRTEI = { 26,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo g_emlrtRTEI = { 29,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo h_emlrtRTEI = { 16,/* lineNo */
  5,                                   /* colNo */
  "abs",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elfun\\abs.m"/* pName */
};

static emlrtRTEInfo i_emlrtRTEI = { 62,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo j_emlrtRTEI = { 63,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo k_emlrtRTEI = { 153,/* lineNo */
  13,                                  /* colNo */
  "find",                              /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pName */
};

static emlrtRTEInfo l_emlrtRTEI = { 69,/* lineNo */
  14,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo m_emlrtRTEI = { 41,/* lineNo */
  5,                                   /* colNo */
  "find",                              /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pName */
};

static emlrtRTEInfo n_emlrtRTEI = { 36,/* lineNo */
  4,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo o_emlrtRTEI = { 69,/* lineNo */
  4,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo p_emlrtRTEI = { 75,/* lineNo */
  16,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo q_emlrtRTEI = { 80,/* lineNo */
  64,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo r_emlrtRTEI = { 30,/* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtRTEInfo s_emlrtRTEI = { 103,/* lineNo */
  1,                                   /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

static emlrtRTEInfo t_emlrtRTEI = { 33,/* lineNo */
  6,                                   /* colNo */
  "find",                              /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pName */
};

static emlrtRTEInfo w_emlrtRTEI = { 88,/* lineNo */
  9,                                   /* colNo */
  "indexShapeCheck",                   /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\indexShapeCheck.m"/* pName */
};

static emlrtRTEInfo x_emlrtRTEI = { 387,/* lineNo */
  1,                                   /* colNo */
  "find",                              /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\lib\\matlab\\elmat\\find.m"/* pName */
};

static emlrtRTEInfo y_emlrtRTEI = { 281,/* lineNo */
  27,                                  /* colNo */
  "cat",                               /* fName */
  "C:\\Program Files\\MATLAB\\R2018a\\toolbox\\eml\\eml\\+coder\\+internal\\cat.m"/* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  80,                                  /* lineNo */
  69,                                  /* colNo */
  "Lrji",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  75,                                  /* lineNo */
  21,                                  /* colNo */
  "Lrji",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  69,                                  /* lineNo */
  17,                                  /* colNo */
  "HT",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo emlrtDCI = { 69,    /* lineNo */
  17,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtRTEInfo ab_emlrtRTEI = { 66,/* lineNo */
  9,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m"/* pName */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  52,                                  /* lineNo */
  59,                                  /* colNo */
  "alphaij",                           /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo e_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  52,                                  /* lineNo */
  52,                                  /* colNo */
  "alphaij",                           /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo f_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  52,                                  /* lineNo */
  40,                                  /* colNo */
  "alphaij",                           /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo g_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  36,                                  /* lineNo */
  20,                                  /* colNo */
  "HT",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo b_emlrtDCI = { 20,  /* lineNo */
  70,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo c_emlrtDCI = { 20,  /* lineNo */
  70,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo d_emlrtDCI = { 20,  /* lineNo */
  53,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo e_emlrtDCI = { 20,  /* lineNo */
  53,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo h_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  42,                                  /* colNo */
  "Lqij",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo i_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  29,                                  /* colNo */
  "Lqij",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo f_emlrtDCI = { 20,  /* lineNo */
  29,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo j_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  17,                                  /* colNo */
  "Lqij",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo g_emlrtDCI = { 20,  /* lineNo */
  17,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo k_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  20,                                  /* lineNo */
  15,                                  /* colNo */
  "Lqij",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo h_emlrtDCI = { 63,  /* lineNo */
  17,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo i_emlrtDCI = { 63,  /* lineNo */
  17,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo j_emlrtDCI = { 63,  /* lineNo */
  28,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo k_emlrtDCI = { 63,  /* lineNo */
  28,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo l_emlrtDCI = { 62,  /* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo m_emlrtDCI = { 62,  /* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtDCInfo n_emlrtDCI = { 63,  /* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtDCInfo o_emlrtDCI = { 63,  /* lineNo */
  1,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  4                                    /* checkKind */
};

static emlrtBCInfo l_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  80,                                  /* lineNo */
  72,                                  /* colNo */
  "Lrji",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo m_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  80,                                  /* lineNo */
  32,                                  /* colNo */
  "Lci",                               /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo p_emlrtDCI = { 80,  /* lineNo */
  32,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo n_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  80,                                  /* lineNo */
  4,                                   /* colNo */
  "LQi",                               /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo q_emlrtDCI = { 80,  /* lineNo */
  4,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo o_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  75,                                  /* lineNo */
  24,                                  /* colNo */
  "Lrji",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo p_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  74,                                  /* lineNo */
  46,                                  /* colNo */
  "Lci",                               /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo r_emlrtDCI = { 74,  /* lineNo */
  46,                                  /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo q_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  75,                                  /* lineNo */
  31,                                  /* colNo */
  "Lrji",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo r_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  75,                                  /* lineNo */
  31,                                  /* colNo */
  "r1",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo s_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  74,                                  /* lineNo */
  7,                                   /* colNo */
  "newLqij",                           /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtDCInfo s_emlrtDCI = { 74,  /* lineNo */
  7,                                   /* colNo */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  1                                    /* checkKind */
};

static emlrtBCInfo t_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  74,                                  /* lineNo */
  7,                                   /* colNo */
  "r1",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo u_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  52,                                  /* lineNo */
  36,                                  /* colNo */
  "alphaij",                           /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo v_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  52,                                  /* lineNo */
  52,                                  /* colNo */
  "c1",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo w_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  55,                                  /* lineNo */
  7,                                   /* colNo */
  "Lrji",                              /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo x_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  55,                                  /* lineNo */
  7,                                   /* colNo */
  "c1",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo y_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  45,                                  /* lineNo */
  16,                                  /* colNo */
  "betaij",                            /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo ab_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  45,                                  /* lineNo */
  16,                                  /* colNo */
  "c1",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo bb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  46,                                  /* lineNo */
  30,                                  /* colNo */
  "betaij",                            /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo cb_emlrtBCI = { -1, /* iFirst */
  -1,                                  /* iLast */
  46,                                  /* lineNo */
  30,                                  /* colNo */
  "c1",                                /* aName */
  "minSum",                            /* fName */
  "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\minSum.m",/* pName */
  0                                    /* checkKind */
};

/* Function Definitions */
void minSum(const emlrtStack *sp, const emxArray_real_T *Lci, const
            emxArray_real_T *prevL, emxArray_real_T *Lqij, const emxArray_real_T
            *HT, real_T Ms, emxArray_real_T *currL, emxArray_real_T *newLqij,
            emxArray_real_T *LQi)
{
  real_T minOfbetaij;
  int32_T i0;
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  int32_T idx;
  int32_T i5;
  boolean_T empty_non_axis_sizes;
  int32_T result;
  int32_T nx;
  cell_wrap_0 reshapes[2];
  int32_T loop_ub;
  emxArray_real_T *betaij;
  emxArray_real_T *Lrji;
  emxArray_real_T *alphaij;
  uint32_T uv0[2];
  int32_T k;
  emxArray_int32_T *c1;
  emxArray_int32_T *ii;
  boolean_T exitg1;
  real_T y;
  emxArray_real_T *r1;
  emxArray_int32_T *b_ii;
  emxArray_real_T *b_Lrji;
  real_T j;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /*  Min-sum product algorithm LDPC decoder */
  /*  */
  /*   Lci      : Received log-likelihood vector (column vector) */
  /*   prevL    : Previous statistics matrix */
  /*   Lqij     : Statistics matrix */
  /*   HT       : LDPC-CC H transpose matrix */
  /*   Ms       : Convolutional code memory */
  /*  */
  /*   currL    : Current statistics matrix */
  /*   newLqij  : Updated statistics matrix */
  /*   LQi      : Final statistics vector (column vector) */
  /*  */
  /*  Bagawan S. Nugroho 2007 */
  /*  Get the H transpose dimension */
  /*  Save the current statistics value for next period decoding */
  minOfbetaij = 2.0 * (Ms + 1.0);
  if (3.0 > minOfbetaij) {
    i0 = 0;
    i1 = 0;
  } else {
    i0 = Lqij->size[0];
    if (!(3 <= i0)) {
      emlrtDynamicBoundsCheckR2012b(3, 1, i0, &k_emlrtBCI, sp);
    }

    i0 = 2;
    i2 = Lqij->size[0];
    if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
      emlrtIntegerCheckR2012b(minOfbetaij, &g_emlrtDCI, sp);
    }

    i1 = (int32_T)minOfbetaij;
    if (!((i1 >= 1) && (i1 <= i2))) {
      emlrtDynamicBoundsCheckR2012b(i1, 1, i2, &j_emlrtBCI, sp);
    }
  }

  if ((Ms + 1.0) + 1.0 > Lqij->size[1]) {
    i2 = 0;
    i3 = 0;
  } else {
    i2 = Lqij->size[1];
    minOfbetaij = (Ms + 1.0) + 1.0;
    if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
      emlrtIntegerCheckR2012b(minOfbetaij, &f_emlrtDCI, sp);
    }

    i4 = (int32_T)minOfbetaij;
    if (!((i4 >= 1) && (i4 <= i2))) {
      emlrtDynamicBoundsCheckR2012b(i4, 1, i2, &i_emlrtBCI, sp);
    }

    i2 = i4 - 1;
    i4 = Lqij->size[1];
    i3 = Lqij->size[1];
    if (!((i3 >= 1) && (i3 <= i4))) {
      emlrtDynamicBoundsCheckR2012b(i3, 1, i4, &h_emlrtBCI, sp);
    }
  }

  st.site = &emlrtRSI;
  minOfbetaij = 2.0 * (Ms + 1.0) - 2.0;
  if (!(minOfbetaij >= 0.0)) {
    emlrtNonNegativeCheckR2012b(minOfbetaij, &e_emlrtDCI, &st);
  }

  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &d_emlrtDCI, &st);
  }

  if (!(Ms + 1.0 >= 0.0)) {
    emlrtNonNegativeCheckR2012b(Ms + 1.0, &c_emlrtDCI, &st);
  }

  minOfbetaij = Ms + 1.0;
  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &b_emlrtDCI, &st);
  }

  b_st.site = &j_emlrtRSI;
  if (!((i1 - i0 == 0) || (i3 - i2 == 0))) {
    idx = i1 - i0;
  } else if (!(((int32_T)(2.0 * (Ms + 1.0) - 2.0) == 0) || ((int32_T)(Ms + 1.0) ==
    0))) {
    idx = (int32_T)(2.0 * (Ms + 1.0) - 2.0);
  } else {
    i5 = i1 - i0;
    idx = muIntScalarMax_sint32(i5, 0);
    if ((int32_T)(2.0 * (Ms + 1.0) - 2.0) > idx) {
      idx = (int32_T)(2.0 * (Ms + 1.0) - 2.0);
    }
  }

  c_st.site = &k_emlrtRSI;
  if ((i1 - i0 == idx) || ((i1 - i0 == 0) || (i3 - i2 == 0))) {
    empty_non_axis_sizes = true;
  } else {
    empty_non_axis_sizes = false;
    emlrtErrorWithMessageIdR2018a(&c_st, &y_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  if (((int32_T)(2.0 * (Ms + 1.0) - 2.0) == idx) || (((int32_T)(2.0 * (Ms + 1.0)
         - 2.0) == 0) || ((int32_T)(Ms + 1.0) == 0))) {
  } else {
    empty_non_axis_sizes = false;
  }

  if (!empty_non_axis_sizes) {
    emlrtErrorWithMessageIdR2018a(&c_st, &y_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  empty_non_axis_sizes = (idx == 0);
  if (empty_non_axis_sizes || (!((i1 - i0 == 0) || (i3 - i2 == 0)))) {
    result = (i3 - i2) - 1;
  } else {
    result = -1;
  }

  if (empty_non_axis_sizes || (!(((int32_T)(2.0 * (Ms + 1.0) - 2.0) == 0) ||
        ((int32_T)(Ms + 1.0) == 0)))) {
    nx = (int32_T)(Ms + 1.0);
  } else {
    nx = 0;
  }

  emxInitMatrix_cell_wrap_0(&b_st, reshapes, &s_emlrtRTEI, true);
  i4 = reshapes[1].f1->size[0] * reshapes[1].f1->size[1];
  reshapes[1].f1->size[0] = idx;
  reshapes[1].f1->size[1] = nx;
  emxEnsureCapacity_real_T(&b_st, reshapes[1].f1, i4, &emlrtRTEI);
  loop_ub = idx * nx;
  for (i4 = 0; i4 < loop_ub; i4++) {
    reshapes[1].f1->data[i4] = 0.0;
  }

  emxInit_real_T(&b_st, &betaij, 2, &r_emlrtRTEI, true);
  i4 = betaij->size[0] * betaij->size[1];
  betaij->size[0] = i1 - i0;
  betaij->size[1] = i3 - i2;
  emxEnsureCapacity_real_T(&b_st, betaij, i4, &b_emlrtRTEI);
  loop_ub = i3 - i2;
  for (i4 = 0; i4 < loop_ub; i4++) {
    nx = i1 - i0;
    for (i3 = 0; i3 < nx; i3++) {
      betaij->data[i3 + betaij->size[0] * i4] = Lqij->data[(i0 + i3) +
        Lqij->size[0] * (i2 + i4)];
    }
  }

  i0 = currL->size[0] * currL->size[1];
  currL->size[0] = idx;
  currL->size[1] = (result + reshapes[1].f1->size[1]) + 1;
  emxEnsureCapacity_real_T(&b_st, currL, i0, &c_emlrtRTEI);
  loop_ub = result + 1;
  for (i0 = 0; i0 < loop_ub; i0++) {
    for (i2 = 0; i2 < idx; i2++) {
      currL->data[i2 + currL->size[0] * i0] = betaij->data[i2 + idx * i0];
    }
  }

  loop_ub = reshapes[1].f1->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    nx = reshapes[1].f1->size[0];
    for (i2 = 0; i2 < nx; i2++) {
      currL->data[i2 + currL->size[0] * ((i0 + result) + 1)] = reshapes[1]
        .f1->data[i2 + reshapes[1].f1->size[0] * i0];
    }
  }

  emxFreeMatrix_cell_wrap_0(&b_st, reshapes);

  /*  Prior log-likelihood. */
  st.site = &b_emlrtRSI;
  b_st.site = &j_emlrtRSI;
  if (!((prevL->size[0] == 0) || (prevL->size[1] == 0))) {
    idx = prevL->size[1];
  } else if (!((Lqij->size[0] == 0) || (Lqij->size[1] == 0))) {
    idx = Lqij->size[1];
  } else {
    idx = muIntScalarMax_sint32(prevL->size[1], 0);
    if (Lqij->size[1] > idx) {
      idx = Lqij->size[1];
    }
  }

  c_st.site = &k_emlrtRSI;
  if ((prevL->size[1] == idx) || ((prevL->size[0] == 0) || (prevL->size[1] == 0)))
  {
    empty_non_axis_sizes = true;
  } else {
    empty_non_axis_sizes = false;
    emlrtErrorWithMessageIdR2018a(&c_st, &y_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  if ((Lqij->size[1] == idx) || ((Lqij->size[0] == 0) || (Lqij->size[1] == 0)))
  {
  } else {
    empty_non_axis_sizes = false;
  }

  if (!empty_non_axis_sizes) {
    emlrtErrorWithMessageIdR2018a(&c_st, &y_emlrtRTEI,
      "MATLAB:catenate:matrixDimensionMismatch",
      "MATLAB:catenate:matrixDimensionMismatch", 0);
  }

  empty_non_axis_sizes = (idx == 0);
  if (empty_non_axis_sizes || (!((prevL->size[0] == 0) || (prevL->size[1] == 0))))
  {
    result = prevL->size[0];
  } else {
    result = 0;
  }

  if (empty_non_axis_sizes || (!((Lqij->size[0] == 0) || (Lqij->size[1] == 0))))
  {
    nx = Lqij->size[0];
  } else {
    nx = 0;
  }

  i0 = betaij->size[0] * betaij->size[1];
  betaij->size[0] = result + nx;
  betaij->size[1] = idx;
  emxEnsureCapacity_real_T(&b_st, betaij, i0, &d_emlrtRTEI);
  for (i0 = 0; i0 < idx; i0++) {
    for (i2 = 0; i2 < result; i2++) {
      betaij->data[i2 + betaij->size[0] * i0] = prevL->data[i2 + result * i0];
    }
  }

  for (i0 = 0; i0 < idx; i0++) {
    for (i2 = 0; i2 < nx; i2++) {
      betaij->data[(i2 + result) + betaij->size[0] * i0] = Lqij->data[i2 + nx *
        i0];
    }
  }

  i0 = Lqij->size[0] * Lqij->size[1];
  Lqij->size[0] = betaij->size[0];
  Lqij->size[1] = betaij->size[1];
  emxEnsureCapacity_real_T(&b_st, Lqij, i0, &e_emlrtRTEI);
  loop_ub = betaij->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    nx = betaij->size[0];
    for (i2 = 0; i2 < nx; i2++) {
      Lqij->data[i2 + Lqij->size[0] * i0] = betaij->data[i2 + betaij->size[0] *
        i0];
    }
  }

  emxInit_real_T(&b_st, &Lrji, 2, &f_emlrtRTEI, true);

  /*  Update the statistics */
  i0 = Lrji->size[0] * Lrji->size[1];
  Lrji->size[0] = HT->size[0];
  Lrji->size[1] = HT->size[1];
  emxEnsureCapacity_real_T(sp, Lrji, i0, &f_emlrtRTEI);
  loop_ub = HT->size[0] * HT->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    Lrji->data[i0] = 0.0;
  }

  emxInit_real_T(sp, &alphaij, 2, &g_emlrtRTEI, true);

  /*  Get the sign and magnitude of L(qij)    */
  i0 = alphaij->size[0] * alphaij->size[1];
  alphaij->size[0] = Lqij->size[0];
  alphaij->size[1] = Lqij->size[1];
  emxEnsureCapacity_real_T(sp, alphaij, i0, &g_emlrtRTEI);
  loop_ub = Lqij->size[0] * Lqij->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    alphaij->data[i0] = Lqij->data[i0];
  }

  st.site = &c_emlrtRSI;
  b_sign(&st, alphaij);
  st.site = &d_emlrtRSI;
  b_st.site = &o_emlrtRSI;
  nx = Lqij->size[0] * Lqij->size[1];
  for (i0 = 0; i0 < 2; i0++) {
    uv0[i0] = (uint32_T)Lqij->size[i0];
  }

  i0 = betaij->size[0] * betaij->size[1];
  betaij->size[0] = (int32_T)uv0[0];
  betaij->size[1] = (int32_T)uv0[1];
  emxEnsureCapacity_real_T(&b_st, betaij, i0, &h_emlrtRTEI);
  c_st.site = &p_emlrtRSI;
  if ((!(1 > nx)) && (nx > 2147483646)) {
    d_st.site = &n_emlrtRSI;
    check_forloop_overflow_error(&d_st);
  }

  for (k = 0; k < nx; k++) {
    betaij->data[k] = muDoubleScalarAbs(Lqij->data[k]);
  }

  /*  ------ Vertical step ------ */
  result = 0;
  emxInit_int32_T(sp, &c1, 1, &n_emlrtRTEI, true);
  emxInit_int32_T(sp, &ii, 1, &t_emlrtRTEI, true);
  while (result <= HT->size[1] - 1) {
    /*  Find non-zeros in the column */
    st.site = &e_emlrtRSI;
    i0 = HT->size[1];
    i2 = result + 1;
    if (!((i2 >= 1) && (i2 <= i0))) {
      emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &g_emlrtBCI, &st);
    }

    b_st.site = &q_emlrtRSI;
    i0 = HT->size[0];
    c_st.site = &r_emlrtRSI;
    idx = 0;
    i2 = HT->size[0];
    i1 = ii->size[0];
    ii->size[0] = i2;
    emxEnsureCapacity_int32_T(&c_st, ii, i1, &k_emlrtRTEI);
    d_st.site = &s_emlrtRSI;
    i2 = HT->size[0];
    if (1 > i2) {
      empty_non_axis_sizes = false;
    } else {
      i2 = HT->size[0];
      empty_non_axis_sizes = (i2 > 2147483646);
    }

    if (empty_non_axis_sizes) {
      e_st.site = &n_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    nx = 1;
    exitg1 = false;
    while ((!exitg1) && (nx <= i0)) {
      if (HT->data[(nx + HT->size[0] * result) - 1] != 0.0) {
        idx++;
        ii->data[idx - 1] = nx;
        if (idx >= i0) {
          exitg1 = true;
        } else {
          nx++;
        }
      } else {
        nx++;
      }
    }

    i0 = HT->size[0];
    if (!(idx <= i0)) {
      emlrtErrorWithMessageIdR2018a(&c_st, &x_emlrtRTEI,
        "Coder:builtins:AssertionFailed", "Coder:builtins:AssertionFailed", 0);
    }

    i0 = HT->size[0];
    if (i0 == 1) {
      if (idx == 0) {
        i0 = ii->size[0];
        ii->size[0] = 0;
        emxEnsureCapacity_int32_T(&c_st, ii, i0, &m_emlrtRTEI);
      }
    } else {
      if (1 > idx) {
        i0 = 0;
      } else {
        i0 = idx;
      }

      d_st.site = &t_emlrtRSI;
      empty_non_axis_sizes = !(ii->size[0] != 1);
      if (empty_non_axis_sizes) {
        empty_non_axis_sizes = false;
        if (i0 != 1) {
          empty_non_axis_sizes = true;
        }

        if (empty_non_axis_sizes) {
          empty_non_axis_sizes = true;
        } else {
          empty_non_axis_sizes = false;
        }
      } else {
        empty_non_axis_sizes = false;
      }

      e_st.site = &u_emlrtRSI;
      if (empty_non_axis_sizes) {
        emlrtErrorWithMessageIdR2018a(&e_st, &w_emlrtRTEI,
          "Coder:FE:PotentialVectorVector", "Coder:FE:PotentialVectorVector", 0);
      }

      i2 = ii->size[0];
      ii->size[0] = i0;
      emxEnsureCapacity_int32_T(&c_st, ii, i2, &m_emlrtRTEI);
    }

    i0 = c1->size[0];
    c1->size[0] = ii->size[0];
    emxEnsureCapacity_int32_T(&st, c1, i0, &n_emlrtRTEI);
    loop_ub = ii->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c1->data[i0] = ii->data[i0];
    }

    /*  Get the minimum of betaij */
    k = 0;
    while (k <= c1->size[0] - 1) {
      /*  Minimum of betaij\c1(k) */
      minOfbetaij = 1.7976931348623157E+308;
      nx = 0;
      while (nx <= c1->size[0] - 1) {
        if (1.0 + (real_T)nx != 1.0 + (real_T)k) {
          i0 = betaij->size[0];
          i2 = c1->size[0];
          if (!((nx + 1 >= 1) && (nx + 1 <= i2))) {
            emlrtDynamicBoundsCheckR2012b(nx + 1, 1, i2, &ab_emlrtBCI, sp);
          }

          i2 = c1->data[nx];
          if (!((i2 >= 1) && (i2 <= i0))) {
            emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &y_emlrtBCI, sp);
          }

          i0 = betaij->size[1];
          i1 = 1 + result;
          if (!((i1 >= 1) && (i1 <= i0))) {
            emlrtDynamicBoundsCheckR2012b(i1, 1, i0, &y_emlrtBCI, sp);
          }

          if (betaij->data[(i2 + betaij->size[0] * (i1 - 1)) - 1] < minOfbetaij)
          {
            i0 = betaij->size[0];
            i2 = c1->size[0];
            if (!((nx + 1 >= 1) && (nx + 1 <= i2))) {
              emlrtDynamicBoundsCheckR2012b(nx + 1, 1, i2, &cb_emlrtBCI, sp);
            }

            i2 = c1->data[nx];
            if (!((i2 >= 1) && (i2 <= i0))) {
              emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &bb_emlrtBCI, sp);
            }

            i0 = betaij->size[1];
            i1 = 1 + result;
            if (!((i1 >= 1) && (i1 <= i0))) {
              emlrtDynamicBoundsCheckR2012b(i1, 1, i0, &bb_emlrtBCI, sp);
            }

            minOfbetaij = betaij->data[(i2 + betaij->size[0] * (i1 - 1)) - 1];
          }
        }

        nx++;
        if (*emlrtBreakCheckR2012bFlagVar != 0) {
          emlrtBreakCheckR2012b(sp);
        }
      }

      /*  for l */
      /*  Multiplication alphaij\c1(k) (use '*' since alphaij are -1/1s) */
      st.site = &f_emlrtRSI;
      nx = alphaij->size[0];
      loop_ub = c1->size[0];
      for (i0 = 0; i0 < loop_ub; i0++) {
        i2 = c1->data[i0];
        if (!((i2 >= 1) && (i2 <= nx))) {
          emlrtDynamicBoundsCheckR2012b(i2, 1, nx, &u_emlrtBCI, &st);
        }
      }

      i0 = alphaij->size[1];
      i2 = result + 1;
      if (!((i2 >= 1) && (i2 <= i0))) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &f_emlrtBCI, &st);
      }

      b_st.site = &v_emlrtRSI;
      if ((c1->size[0] == 1) || (c1->size[0] != 1)) {
      } else {
        emlrtErrorWithMessageIdR2018a(&b_st, &v_emlrtRTEI,
          "Coder:toolbox:autoDimIncompatibility",
          "Coder:toolbox:autoDimIncompatibility", 0);
      }

      c_st.site = &w_emlrtRSI;
      nx = c1->size[0];
      if (c1->size[0] == 0) {
        y = 1.0;
      } else {
        d_st.site = &x_emlrtRSI;
        y = alphaij->data[(c1->data[0] + alphaij->size[0] * result) - 1];
        e_st.site = &y_emlrtRSI;
        empty_non_axis_sizes = ((!(2 > c1->size[0])) && (c1->size[0] >
          2147483646));
        if (empty_non_axis_sizes) {
          f_st.site = &n_emlrtRSI;
          check_forloop_overflow_error(&f_st);
        }

        for (idx = 2; idx <= nx; idx++) {
          y *= alphaij->data[(c1->data[idx - 1] + alphaij->size[0] * result) - 1];
        }
      }

      i0 = alphaij->size[0];
      i2 = c1->size[0];
      if (!((k + 1 >= 1) && (k + 1 <= i2))) {
        emlrtDynamicBoundsCheckR2012b(k + 1, 1, i2, &v_emlrtBCI, sp);
      }

      i2 = c1->data[k];
      if (!((i2 >= 1) && (i2 <= i0))) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &e_emlrtBCI, sp);
      }

      i0 = alphaij->size[1];
      i2 = result + 1;
      if (!((i2 >= 1) && (i2 <= i0))) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &d_emlrtBCI, sp);
      }

      /*  Update L(rji) */
      i0 = Lrji->size[0];
      i2 = c1->size[0];
      if (!((k + 1 >= 1) && (k + 1 <= i2))) {
        emlrtDynamicBoundsCheckR2012b(k + 1, 1, i2, &x_emlrtBCI, sp);
      }

      i2 = c1->data[k];
      if (!((i2 >= 1) && (i2 <= i0))) {
        emlrtDynamicBoundsCheckR2012b(i2, 1, i0, &w_emlrtBCI, sp);
      }

      i0 = Lrji->size[1];
      i1 = 1 + result;
      if (!((i1 >= 1) && (i1 <= i0))) {
        emlrtDynamicBoundsCheckR2012b(i1, 1, i0, &w_emlrtBCI, sp);
      }

      Lrji->data[(i2 + Lrji->size[0] * (i1 - 1)) - 1] = y * alphaij->data
        [(c1->data[k] + alphaij->size[0] * result) - 1] * minOfbetaij;
      k++;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    /*  for k */
    result++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_int32_T(sp, &ii);
  emxFree_int32_T(sp, &c1);
  emxFree_real_T(sp, &betaij);
  emxFree_real_T(sp, &alphaij);

  /*  for i */
  /*  ----- Horizontal step ----- */
  i0 = LQi->size[0];
  minOfbetaij = 2.0 * (Ms + 1.0);
  if (!(minOfbetaij >= 0.0)) {
    emlrtNonNegativeCheckR2012b(minOfbetaij, &m_emlrtDCI, sp);
  }

  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &l_emlrtDCI, sp);
  }

  LQi->size[0] = (int32_T)minOfbetaij;
  emxEnsureCapacity_real_T1(sp, LQi, i0, &i_emlrtRTEI);
  minOfbetaij = 2.0 * (Ms + 1.0);
  if (!(minOfbetaij >= 0.0)) {
    emlrtNonNegativeCheckR2012b(minOfbetaij, &m_emlrtDCI, sp);
  }

  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &l_emlrtDCI, sp);
  }

  loop_ub = (int32_T)minOfbetaij;
  for (i0 = 0; i0 < loop_ub; i0++) {
    LQi->data[i0] = 0.0;
  }

  i0 = newLqij->size[0] * newLqij->size[1];
  minOfbetaij = 2.0 * (Ms + 1.0);
  if (!(minOfbetaij >= 0.0)) {
    emlrtNonNegativeCheckR2012b(minOfbetaij, &i_emlrtDCI, sp);
  }

  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &h_emlrtDCI, sp);
  }

  newLqij->size[0] = (int32_T)minOfbetaij;
  minOfbetaij = 2.0 * (Ms + 1.0) - 1.0;
  if (!(minOfbetaij >= 0.0)) {
    emlrtNonNegativeCheckR2012b(minOfbetaij, &k_emlrtDCI, sp);
  }

  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &j_emlrtDCI, sp);
  }

  newLqij->size[1] = (int32_T)minOfbetaij;
  emxEnsureCapacity_real_T(sp, newLqij, i0, &j_emlrtRTEI);
  minOfbetaij = 2.0 * (Ms + 1.0);
  if (!(minOfbetaij >= 0.0)) {
    emlrtNonNegativeCheckR2012b(minOfbetaij, &o_emlrtDCI, sp);
  }

  if (minOfbetaij != (int32_T)muDoubleScalarFloor(minOfbetaij)) {
    emlrtIntegerCheckR2012b(minOfbetaij, &n_emlrtDCI, sp);
  }

  y = 2.0 * (Ms + 1.0) - 1.0;
  if (!(y >= 0.0)) {
    emlrtNonNegativeCheckR2012b(y, &o_emlrtDCI, sp);
  }

  if (y != (int32_T)muDoubleScalarFloor(y)) {
    emlrtIntegerCheckR2012b(y, &n_emlrtDCI, sp);
  }

  loop_ub = (int32_T)minOfbetaij * (int32_T)y;
  for (i0 = 0; i0 < loop_ub; i0++) {
    newLqij->data[i0] = 0.0;
  }

  /*  Process the statistics within Ms + 1 period */
  minOfbetaij = 2.0 * (Ms + 1.0) - 1.0;
  y = 4.0 * (Ms + 1.0) - 2.0;
  i0 = (int32_T)(y + (1.0 - minOfbetaij));
  emlrtForLoopVectorCheckR2012b(minOfbetaij, 1.0, y, mxDOUBLE_CLASS, i0,
    &ab_emlrtRTEI, sp);
  result = 0;
  emxInit_real_T(sp, &r1, 2, &o_emlrtRTEI, true);
  emxInit_int32_T1(sp, &b_ii, 2, &t_emlrtRTEI, true);
  emxInit_real_T(sp, &b_Lrji, 2, &q_emlrtRTEI, true);
  while (result <= i0 - 1) {
    j = minOfbetaij + (real_T)result;

    /*  Find non-zero in the row */
    st.site = &g_emlrtRSI;
    loop_ub = HT->size[1];
    i2 = HT->size[0];
    if (j != (int32_T)muDoubleScalarFloor(j)) {
      emlrtIntegerCheckR2012b(j, &emlrtDCI, &st);
    }

    nx = (int32_T)j;
    if (!((nx >= 1) && (nx <= i2))) {
      emlrtDynamicBoundsCheckR2012b(nx, 1, i2, &c_emlrtBCI, &st);
    }

    i2 = r1->size[0] * r1->size[1];
    r1->size[0] = 1;
    r1->size[1] = loop_ub;
    emxEnsureCapacity_real_T(&st, r1, i2, &l_emlrtRTEI);
    for (i2 = 0; i2 < loop_ub; i2++) {
      r1->data[r1->size[0] * i2] = HT->data[(nx + HT->size[0] * i2) - 1];
    }

    b_st.site = &q_emlrtRSI;
    i2 = HT->size[1];
    c_st.site = &r_emlrtRSI;
    idx = 0;
    i1 = HT->size[1];
    i4 = b_ii->size[0] * b_ii->size[1];
    b_ii->size[0] = 1;
    b_ii->size[1] = i1;
    emxEnsureCapacity_int32_T1(&c_st, b_ii, i4, &k_emlrtRTEI);
    d_st.site = &s_emlrtRSI;
    i1 = HT->size[1];
    if (1 > i1) {
      empty_non_axis_sizes = false;
    } else {
      i1 = HT->size[1];
      empty_non_axis_sizes = (i1 > 2147483646);
    }

    if (empty_non_axis_sizes) {
      e_st.site = &n_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    nx = 1;
    exitg1 = false;
    while ((!exitg1) && (nx <= i2)) {
      if (r1->data[nx - 1] != 0.0) {
        idx++;
        b_ii->data[idx - 1] = nx;
        if (idx >= i2) {
          exitg1 = true;
        } else {
          nx++;
        }
      } else {
        nx++;
      }
    }

    i2 = HT->size[1];
    if (!(idx <= i2)) {
      emlrtErrorWithMessageIdR2018a(&c_st, &x_emlrtRTEI,
        "Coder:builtins:AssertionFailed", "Coder:builtins:AssertionFailed", 0);
    }

    i2 = HT->size[1];
    if (i2 == 1) {
      if (idx == 0) {
        i2 = b_ii->size[0] * b_ii->size[1];
        b_ii->size[0] = 1;
        b_ii->size[1] = 0;
        emxEnsureCapacity_int32_T1(&c_st, b_ii, i2, &m_emlrtRTEI);
      }
    } else {
      i2 = b_ii->size[0] * b_ii->size[1];
      if (1 > idx) {
        b_ii->size[1] = 0;
      } else {
        b_ii->size[1] = idx;
      }

      emxEnsureCapacity_int32_T1(&c_st, b_ii, i2, &m_emlrtRTEI);
    }

    i2 = r1->size[0] * r1->size[1];
    r1->size[0] = 1;
    r1->size[1] = b_ii->size[1];
    emxEnsureCapacity_real_T(&st, r1, i2, &o_emlrtRTEI);
    loop_ub = b_ii->size[0] * b_ii->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      r1->data[i2] = b_ii->data[i2];
    }

    k = 1;
    while (k - 1 <= r1->size[1] - 1) {
      /*  Update L(qij) by summation of L(rij)\r1(k) */
      idx = Lrji->size[1];
      i2 = Lrji->size[0];
      nx = (int32_T)j;
      if (!((nx >= 1) && (nx <= i2))) {
        emlrtDynamicBoundsCheckR2012b(nx, 1, i2, &b_emlrtBCI, sp);
      }

      i2 = b_Lrji->size[0] * b_Lrji->size[1];
      b_Lrji->size[0] = 1;
      b_Lrji->size[1] = r1->size[1];
      emxEnsureCapacity_real_T(sp, b_Lrji, i2, &p_emlrtRTEI);
      loop_ub = r1->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        i1 = (int32_T)r1->data[r1->size[0] * i2];
        if (!((i1 >= 1) && (i1 <= idx))) {
          emlrtDynamicBoundsCheckR2012b(i1, 1, idx, &o_emlrtBCI, sp);
        }

        b_Lrji->data[b_Lrji->size[0] * i2] = Lrji->data[(nx + Lrji->size[0] *
          (i1 - 1)) - 1];
      }

      i2 = Lci->size[0];
      y = j - (2.0 * (Ms + 1.0) - 2.0);
      if (y != (int32_T)muDoubleScalarFloor(y)) {
        emlrtIntegerCheckR2012b(y, &r_emlrtDCI, sp);
      }

      i1 = (int32_T)y;
      if (!((i1 >= 1) && (i1 <= i2))) {
        emlrtDynamicBoundsCheckR2012b(i1, 1, i2, &p_emlrtBCI, sp);
      }

      i2 = Lrji->size[0];
      i4 = (int32_T)j;
      if (!((i4 >= 1) && (i4 <= i2))) {
        emlrtDynamicBoundsCheckR2012b(i4, 1, i2, &q_emlrtBCI, sp);
      }

      i2 = Lrji->size[1];
      i3 = r1->size[1];
      if (!((k >= 1) && (k <= i3))) {
        emlrtDynamicBoundsCheckR2012b(k, 1, i3, &r_emlrtBCI, sp);
      }

      i3 = (int32_T)r1->data[k - 1];
      if (!((i3 >= 1) && (i3 <= i2))) {
        emlrtDynamicBoundsCheckR2012b(i3, 1, i2, &q_emlrtBCI, sp);
      }

      i2 = newLqij->size[0];
      y = j - (2.0 * (Ms + 1.0) - 2.0);
      if (y != (int32_T)muDoubleScalarFloor(y)) {
        emlrtIntegerCheckR2012b(y, &s_emlrtDCI, sp);
      }

      nx = (int32_T)y;
      if (!((nx >= 1) && (nx <= i2))) {
        emlrtDynamicBoundsCheckR2012b(nx, 1, i2, &s_emlrtBCI, sp);
      }

      i2 = newLqij->size[1];
      idx = r1->size[1];
      if (!((k >= 1) && (k <= idx))) {
        emlrtDynamicBoundsCheckR2012b(k, 1, idx, &t_emlrtBCI, sp);
      }

      idx = (int32_T)r1->data[k - 1];
      if (!((idx >= 1) && (idx <= i2))) {
        emlrtDynamicBoundsCheckR2012b(idx, 1, i2, &s_emlrtBCI, sp);
      }

      st.site = &h_emlrtRSI;
      newLqij->data[(nx + newLqij->size[0] * (idx - 1)) - 1] = (Lci->data[i1 - 1]
        + sum(&st, b_Lrji)) - Lrji->data[(i4 + Lrji->size[0] * (i3 - 1)) - 1];
      k++;
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }

    /*  for k */
    /*  Final statistics for decision */
    idx = Lrji->size[1];
    i2 = Lrji->size[0];
    nx = (int32_T)j;
    if (!((nx >= 1) && (nx <= i2))) {
      emlrtDynamicBoundsCheckR2012b(nx, 1, i2, &emlrtBCI, sp);
    }

    i2 = b_Lrji->size[0] * b_Lrji->size[1];
    b_Lrji->size[0] = 1;
    b_Lrji->size[1] = r1->size[1];
    emxEnsureCapacity_real_T(sp, b_Lrji, i2, &q_emlrtRTEI);
    loop_ub = r1->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      i1 = (int32_T)r1->data[r1->size[0] * i2];
      if (!((i1 >= 1) && (i1 <= idx))) {
        emlrtDynamicBoundsCheckR2012b(i1, 1, idx, &l_emlrtBCI, sp);
      }

      b_Lrji->data[b_Lrji->size[0] * i2] = Lrji->data[(nx + Lrji->size[0] * (i1
        - 1)) - 1];
    }

    i2 = Lci->size[0];
    y = j - (2.0 * (Ms + 1.0) - 2.0);
    if (y != (int32_T)muDoubleScalarFloor(y)) {
      emlrtIntegerCheckR2012b(y, &p_emlrtDCI, sp);
    }

    i1 = (int32_T)y;
    if (!((i1 >= 1) && (i1 <= i2))) {
      emlrtDynamicBoundsCheckR2012b(i1, 1, i2, &m_emlrtBCI, sp);
    }

    i2 = LQi->size[0];
    y = j - (2.0 * (Ms + 1.0) - 2.0);
    if (y != (int32_T)muDoubleScalarFloor(y)) {
      emlrtIntegerCheckR2012b(y, &q_emlrtDCI, sp);
    }

    i4 = (int32_T)y;
    if (!((i4 >= 1) && (i4 <= i2))) {
      emlrtDynamicBoundsCheckR2012b(i4, 1, i2, &n_emlrtBCI, sp);
    }

    st.site = &i_emlrtRSI;
    LQi->data[i4 - 1] = Lci->data[i1 - 1] + sum(&st, b_Lrji);
    result++;
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(sp);
    }
  }

  emxFree_real_T(sp, &b_Lrji);
  emxFree_int32_T(sp, &b_ii);
  emxFree_real_T(sp, &r1);
  emxFree_real_T(sp, &Lrji);

  /*  for j */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (minSum.c) */
