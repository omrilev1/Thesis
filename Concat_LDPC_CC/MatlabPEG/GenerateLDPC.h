#ifndef GENERATELDPC
#define GENERATELDPC



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "BigGirth.h"
#include "Random.h"
#include "CyclesOfGraph.h"

using namespace std;
/*
 * M = number of rows
 * N = number of columns
 * d = degree of regular graph
*/
bool GenerateLDPC(int M, int N, char *  codeName, int d);


#endif
