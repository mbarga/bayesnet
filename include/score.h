#ifndef _SCORE_H_
#define _SCORE_H_

//TODO use of CONST INTS
void * 	bde_init(double*, int, int, int, int);
void    bde_destroy_buff(void *);
double 	bde_score(void*, int*, int);

#endif
