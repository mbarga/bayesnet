/*
  BDE score function

  Copyright (C) 2009, 2013 Yoshinori Tamada
*/

#ifndef __BDE_H
#define __BDE_H

void *      BDE_init(double * X, double * Y, int p, int n, int c, int m);
double      BDE_score(void * buff, const int j, const int * parents, const int q);
void        BDE_finalize(void * buff);

#endif /* __BDE_H */

/* End of file.  Copyright (C) 2009 HGC, IMS, Univeristy of Tokyo */
