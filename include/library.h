#ifndef _SUPP_H
#define _SUPP_H

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <math.h>

int * randperm(int, int);
void one_to_one(GHashTable **, double **);
double score(NODE *, NODE *);

#endif

