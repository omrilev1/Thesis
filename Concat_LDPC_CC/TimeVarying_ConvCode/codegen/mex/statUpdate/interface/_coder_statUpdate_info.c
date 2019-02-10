/*
 * _coder_statUpdate_info.c
 *
 * Code generation for function '_coder_statUpdate_info'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "statUpdate.h"
#include "_coder_statUpdate_info.h"

/* Function Definitions */
mxArray *emlrtMexFcnProperties(void)
{
  mxArray *xResult;
  mxArray *xEntryPoints;
  const char * fldNames[6] = { "Name", "NumberOfInputs", "NumberOfOutputs",
    "ConstantInputs", "FullPath", "TimeStamp" };

  mxArray *xInputs;
  const char * b_fldNames[4] = { "Version", "ResolvedFunctions", "EntryPoints",
    "CoverageInfo" };

  xEntryPoints = emlrtCreateStructMatrix(1, 1, 6, fldNames);
  xInputs = emlrtCreateLogicalMatrix(1, 4);
  emlrtSetField(xEntryPoints, 0, "Name", emlrtMxCreateString("statUpdate"));
  emlrtSetField(xEntryPoints, 0, "NumberOfInputs", emlrtMxCreateDoubleScalar(4.0));
  emlrtSetField(xEntryPoints, 0, "NumberOfOutputs", emlrtMxCreateDoubleScalar
                (2.0));
  emlrtSetField(xEntryPoints, 0, "ConstantInputs", xInputs);
  emlrtSetField(xEntryPoints, 0, "FullPath", emlrtMxCreateString(
    "C:\\Users\\OmriLev\\Documents\\GitHub\\Thesis\\Concat_LDPC_CC\\TimeVarying_ConvCode\\statUpdate.m"));
  emlrtSetField(xEntryPoints, 0, "TimeStamp", emlrtMxCreateDoubleScalar
                (737459.86625));
  xResult = emlrtCreateStructMatrix(1, 1, 4, b_fldNames);
  emlrtSetField(xResult, 0, "Version", emlrtMxCreateString(
    "9.4.0.813654 (R2018a)"));
  emlrtSetField(xResult, 0, "ResolvedFunctions", (mxArray *)
                emlrtMexFcnResolvedFunctionsInfo());
  emlrtSetField(xResult, 0, "EntryPoints", xEntryPoints);
  return xResult;
}

const mxArray *emlrtMexFcnResolvedFunctionsInfo(void)
{
  const mxArray *nameCaptureInfo;
  const char * data[5] = {
    "789ce554cb4e0231142dc6e742435cb866e5c684aa2be34a1c148c8312198d891a2cc355aad396b485a01b7fc4ad7b3f81aff01bfc05770e0c953a6182d1e8c6"
    "9b90dbc3197ace3db9034aed955208a10514556626eaf3039c1ef409f4b9e27c2aa19b9a42939f7e67f8c741f705d7d0d1110828878316ab810c01270c3eaea9",
    "0b4639e1dabb6b0292a044d0867a9fb9a2017894812b2c50a42160bb16f5017a54efec34c0bfadb418920d35b41bd80059f974adf9df1e86f34f8ec8c7e6cdbc"
    "e709f9a463fcd9ce85b3898f1548850f99a42eb4715ef82d065c2b5ca0bad8aa61af018a2aec08ee135d75f365a7ea38b837dd099177945f5743aaed883a60a5",
    "893e6ed689862cb3e7695a7e6d1fd323e6b179e37f22864dcdc59e8f2ab315f5e52da3df49b87f549ea3f41713f4d331de0f4390591a2e99e424c8867945df1b"
    "1f97dff411af241fa68cdef337f5ccfdfb63f40c7fb6e79ef657a92cc5b5242cd3db79854b39cfcd6de3a3f5d5b50d82b510414d7430b0a0ff59e9a785574c5c",
    "388c2b5a9bb1797d756f92fe27e6d0ece0d47d7de22fb9bfd38beabfe8fdf4bd5b4ad04bc7f8fb7ca9e6b9c57d7213acb755a5595457371b85a18ff2189d713e"
    "5002feedfbdf01df7d867a", "" };

  nameCaptureInfo = NULL;
  emlrtNameCaptureMxArrayR2016a(data, 1848U, &nameCaptureInfo);
  return nameCaptureInfo;
}

/* End of code generation (_coder_statUpdate_info.c) */
