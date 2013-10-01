#include <stdio.h>
#include "bdescore.h"

int
main(int argc, char * argv[])
{
  /* (6,3)-matrix */
  double X[] =
    {   0, 1, 0, 1, 0, 1,
        0, 0, 1, 1, 0, 1,
        0, 1, 1, 0, 1, 0 };
  int p = 3;
  int n = 6;
  int parents[] =
    { 1, 2 };
  int c = 2;
  int m = 2;

  void * buff = bde_init(X, p, n, c, m);

  double s = get_score(buff, parents, 2);
  printf("score = %f\n", s);

  return 0;
}

