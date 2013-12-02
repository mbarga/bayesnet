#ifndef _MAIN_H_
#define _MAIN_H_

#include "node.h"

//TODO make some constant?
typedef struct params {
	double *X;	// n x m dataset
	NODE *Y;	// parent edges for each node
	int p;		// number of nodes
	int n;		// number of samples per node
	int r;		// number of discrete bins
	int m;		// max candidate parent set size
	int max_parents;// maximum number of parents per node
} PARAMS;

#endif
