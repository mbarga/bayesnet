#include "util.h"
//#include "score.h"
#include "BDE.h"

#include <stdlib.h>
#include <stdio.h>

double get_score(void *buff, int child, int *parents, int size)
{
	//return (double)bde_score(buff, parents, size);
	return (double)BDE_score(buff, child, parents, size);
}

void
util_print_dmatrix(double *c, int sizen, int sizem)
{
  for (int i = 0; i < sizen; ++i) {
    for (int j = 0; j < sizem; ++j)
      printf("%f, ", c[i * sizem + j]);
    printf("\n");
  }
}

void
util_print_imatrix(int *c, int size)
{
  for (int i = 0; i < size; ++i) {
    for (int j = 0; j < size; ++j)
      printf("%d, ", c[i * size + j]);
    printf("\n");
  }
}

void
util_errlog(char *s)
{
    fprintf(stderr, "%s", s);
    //syslog(LOG_INFO, "%s", "randperm(): m was found NULL\n");
}
