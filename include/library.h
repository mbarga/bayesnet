#ifndef _LIBRARY_H_
#define _LIBRARY_H_

void randperm(int*, int, int);
void scramble(int**, int);
double * max(double *, double*, double*);
void one_to_one(NODE **, int, double **);
double score(NODE *, NODE *);
int generate_local_probability();
void free_nodes(NODE *, int);
void print_nodes(NODE *, int);
void print_matrix(int *, int);
void errlog(char *);

#endif

