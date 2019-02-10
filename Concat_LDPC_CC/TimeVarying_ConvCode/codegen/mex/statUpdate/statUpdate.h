/*
 * statUpdate.h
 *
 * Code generation for function 'statUpdate'
 *
 */

#ifndef STATUPDATE_H
#define STATUPDATE_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "statUpdate_types.h"

/* Function Declarations */
extern void statUpdate(const emlrtStack *sp, real_T numOfChuncks, real_T Ms,
  const emxArray_real_T *newL, const emxArray_real_T *Lqij_init, emxArray_real_T
  *Lqij, emxArray_real_T *nextL);

#endif

/* End of code generation (statUpdate.h) */
