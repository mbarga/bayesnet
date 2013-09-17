#ifndef _BDESCORE_H_
#define _BDESCORE_H_

typedef struct buffer
{
	int node_count;
	int sample_count;
	double *data;

	int prior;
	int categories;

	int *n_ij;
	int *n_ijk;

	double fixed_alpha;

	int max_parents;
} BUFFER;

//TODO use of CONST INTS
void* 	bde_init(double, int, int, int);
double 	get_score(void*, int, int*, int);
int 	count_nijk(BUFFER*, int*, int);
double 	calc_bde(int*, double, int*, int, int, int);

#endif
