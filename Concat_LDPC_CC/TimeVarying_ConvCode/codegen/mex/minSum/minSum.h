/*
 * minSum.h
 *
 * Code generation for function 'minSum'
 *
 */

#ifndef MINSUM_H
#define MINSUM_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "minSum_types.h"

/* Function Declarations */
extern void minSum(const emlrtStack *sp, const emxArray_real_T *Lci, const
                   emxArray_real_T *prevL, emxArray_real_T *Lqij, const
                   emxArray_real_T *HT, real_T Ms, emxArray_real_T *currL,
                   emxArray_real_T *newLqij, emxArray_real_T *LQi);

#endif

/* End of code generation (minSum.h) */
