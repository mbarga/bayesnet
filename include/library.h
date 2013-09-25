#ifndef _LIBRARY_H_
#define _LIBRARY_H_

void
randperm(int *, int, int);
void
scramble(int **, int);
double*
max(double *, double*, double*);
void
one_to_one(double *, int, double **, int, int);
double
score(NODE *, NODE *);
int
generate_local_probability();
//void
//print_nodes(NODE *, int);
void
errlog(char *);

#endif

