/*
 * GenerateLDPC_mex.C
 * This is a mex function for the CPP 
 * function GenerateLDPC.C
 * Author: Hatef Monajemi (monajemi@stanford.edu)
 * Date: 2013 Feb 23
 */

#include "mex.h"  
#include "GenerateLDPC.h"
#include "matrix.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){

//declare variables
int M, N, d;
char* codeName;
int StrLength;
bool output;


//Get integer values of M, N and d;
M = (int) mxGetScalar(prhs[0]);
N = (int) mxGetScalar(prhs[1]);
d = (int) mxGetScalar(prhs[3]);

//Get the string codeName
StrLength=mxGetN(prhs[2])+1;
codeName=(char*) mxCalloc(StrLength,sizeof(char));
mxGetString(prhs[2], codeName, StrLength);

//Call GenerateLDPC.C
GenerateLDPC(M, N, codeName, d);


mxFree(codeName);

return;
}
