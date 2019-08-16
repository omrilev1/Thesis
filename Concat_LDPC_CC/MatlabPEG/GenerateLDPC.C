

/* This is a function that calls 
 * the code of Xiaoyu Hu, et al
 * to generate LDPC codes using 
 * Progressive Edge Growth
 * Author: Hatef Monajemi (monajemi@stanford.edu)
 * Date : 11 Feb 2013 
 * PS: We are interested ONLY in regular graphs
 * therefore a degree value 'd' would suffice and
 * we do not require degFileName
*/


#include "GenerateLDPC.h"
#include <assert.h>

using namespace std;
/*
 * M = number of rows
 * N = number of columns
*/
bool GenerateLDPC(int M, int N, char* codeName, int d){

// Check the size of graph 
assert( M < N && "M must be smaller than N");

//cout << "M(#rows):" << M << endl;
//cout << "N(#cols):" << N << endl;

int *degSeq ;
degSeq = new int[N];

// set all the enteries of degSeq to d
for(int i=0;i<N; ++i){
degSeq[i] = d;
}

int sglConcent=1;  // default to non-strictly concentrated parity-check distribution
int targetGirth=100000; // default to greedy PEG version 
  
BigGirth *bigGirth = new BigGirth(M, N, degSeq, codeName, sglConcent, targetGirth);

(*bigGirth).writeToFile_Hcompressed();

delete [] degSeq;

return true;
}
