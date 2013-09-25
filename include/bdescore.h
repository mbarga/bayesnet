#ifndef _BDESCORE_H_
#define _BDESCORE_H_

//TODO use of CONST INTS
void * 	bde_init(double*, int, int, int, int);
void    bde_destroy_buff(void *);
double 	get_score(void*, int*, int);

#endif
