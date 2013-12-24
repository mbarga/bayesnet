#ifndef _SCORE_H_
#define _SCORE_H_

//TODO use of CONST INTS
void *	score_init(double*, int, int, int, int);
void    score_destroy_buff(void *);
double	score_calc(void *, const int, const int *, const int);

#endif
