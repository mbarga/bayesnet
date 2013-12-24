/*
   BDE score function

   Copyright (C) 2009, 2013 Yoshinori Tamada
   */

/*!\defgroup BDE
 *\brief BDE score function
 *
 *\author Yoshinori Tamada &lt;tamada at is.s.u-tokyo.ac.jp&gt;
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#define BDE_PRIOR_TYPE_FIXED 0

typedef struct
{
	int p; /* The number of nodes */
	int n; /* The number of samples */
	double * X; /* (n,p)-fortran style input data matrix for parents */
	double * Y; /* (n,p)-fortran style input data matrix for children */

	/* NOTE: X[i + j * n] corresponds to the (i,j) element in the matrix */

	int prior; /* Prior type. Currently, only BDE_PRIOR_TYPE_FIXED is supported. */
	int * r; /* r[i] = # of categories for the i-th variable (node) */

	int * n_ij; /* working area */
	int * n_ijk; /* working area */

	double fixed_alpha; /* fixed alpha parameter */

	int max_parents; /* The number of max parents */
	int verbose; /* verbose mode if non zero. */

} BDEBuff;

/* Prototype declaration of the internal functions */
BDEBuff *
BDEBuff_new();
int
BDE_count_n_ijk(BDEBuff * bf, const int * parents, int q, int i);
double
BDE_calc_bde(int * n_ijk, double fixed_alpha, int i, const int * parents, int q,
		int count, int * r);

/*****************************************************************************/
/* M A I N   R O U T I N E S                                                 */
/*****************************************************************************/

/*!\brief Initializes the BDE score function.
 *
 *\param X  (\a n, \a p)-fortran-style input data matrix for parents.
 *\param Y  (\a n, \a p)-fortran-style input data matrix for children.
 *\param p  The number of nodes.
 *\param n  The number of samples.
 *\param c  The number of categories (types of discretized values) in data.
 *\param m  The number of maximum parents.
 *
 *\return  Pointer to the score working buffer.
 *\ingroup BDE */
void *
BDE_init(double * X, double * Y, int p, int n, int c, int m)
{
	BDEBuff * bf = BDEBuff_new();
	size_t mem = sizeof(BDEBuff);

	bf->p = p;
	bf->n = n;
	bf->X = X;
	bf->Y = Y;

	bf->max_parents = m;
	bf->verbose = -1;

	if (bf->verbose >= 0)
	{
		fprintf(stderr, "BDE: Initializing the BDe score function.\n");
	}

	/* ------------------------------------------------------------- */
	/* Categories : bf->r_i                                          */
	/* ------------------------------------------------------------- */

	mem = sizeof(int) * bf->p;
	bf->r = (int *) malloc(mem);

	int max_r_i = 0;

	/* The program supports different numbers of categories for
	   variables.  Here the same number of categories is assumed.
	   You can change this by yourself. */
	for (int i = 0; i < bf->p; i++)
	{
		bf->r[i] = c;
		if (bf->r[i] > max_r_i)
			max_r_i = bf->r[i];
	}/* for i < bf->n */

	/* ------------------------------------------------------------- */
	/* n_ij & n_ijk                                                  */
	/* ------------------------------------------------------------- */
	/*
Note: The max num of patterns to count is equal to the number of
samples because if the patterns in all the samples are different,
then all the parent status patterns are different and they are
the entire possible patterns from the data.
*/
	mem = sizeof(int) * bf->max_parents * bf->n;
	bf->n_ij = (int *) malloc(mem);
	mem = sizeof(int) * bf->n * max_r_i;
	bf->n_ijk = (int *) malloc(mem);

	/* Output summary */
	if (bf->verbose >= 0)
	{
		fprintf(stderr, "BDE: Total memory usage: %zu [bytes]", mem);
		fprintf(stderr, " (%zu [MiB])\n", (mem + 524288) / 1048576);
		fprintf(stderr, "BDE: Finished BDE initialization.\n");
		fflush(stderr);
	}

	return bf;
}
/* BDE_init() ============================================================== */

double
BDE_score(void * buff, const int j, const int * parents, const int q)
{
	BDEBuff * bf = (BDEBuff *) buff;

	/* likelihood */
	int count = BDE_count_n_ijk(bf, parents, q, j);
#ifdef SIGN_DEBUG
	fprintf(SIGN_ERR, "BDE_score(): count = %d\n", count);
#endif
	double score = BDE_calc_bde(bf->n_ijk, bf->fixed_alpha, j, parents, q, count,
			bf->r);

	/* prior */
	switch (bf->prior)
	{
	case BDE_PRIOR_TYPE_FIXED:
		score = -2 * score;
		break;

	default:
		fprintf(stderr, "FATAL ERROR: BDE: Unknown prior type specified.\n");
		exit(1);
	}

	return score;
}
/* BDE_score() ============================================================= */

/*!\brief Releases memory allocated for the working buffer.
 *
 *\param buff Pointer to the working buffer returned by BDE_init().
 */
void
BDE_finalize(void * buff)
{
	BDEBuff * bf = (BDEBuff *) buff;

	free(bf->n_ij);
	free(bf->n_ijk);
	free(bf->r);

	free(buff);
}
/* BDE_finalize() ========================================================== */

/*****************************************************************************/
/* S U B   R O U T I N E S                                                   */
/*****************************************************************************/

BDEBuff *
BDEBuff_new()
{
	BDEBuff * bf = (BDEBuff *) malloc(sizeof(BDEBuff));

	bf->prior = BDE_PRIOR_TYPE_FIXED;
	bf->fixed_alpha = 1.0;

	return bf;
}
/* BDEBuff_new() =========================================================== */

/*
parents : parent ID vector

return: number of patterns of 'j' (parent status)

Using bf->n_ij & bf->n_jik.
*/
int
BDE_count_n_ijk(BDEBuff * bf, const int * parents, int q, int i)
{

	int n_ij_count = 0;

	for (int l = 0; l < bf->n; l++)
	{
		/* j[m] : parent stat vector */
		int j[q];
		for (int m = 0; m < q; m++)  {
			j[m] = bf->X[l + parents[m] * bf->n];
		}

		/* Find 'stat vector' from n_ij */
		int s = 0; /* stat vector position in n_ij */
		while (s < n_ij_count)
		{
			const int * vec = bf->n_ij + s * bf->max_parents;
			int flag2 = 0;
			for (int ii = 0; ii < q; ii++)
			{
				if (j[ii] != vec[ii])
				{
					flag2 = -1;
					break;
				}
			}
			if (flag2 == 0)
				break; /* found */

			s++;
		}/* while(s < n_ij_count) */

		if (s == n_ij_count)
		{
			/* Not found.  Creats a new element */
			memcpy(bf->n_ij + s * bf->max_parents, j, sizeof(int) * q);
			n_ij_count++;

			bzero(bf->n_ijk + s * bf->r[i], sizeof(int) * bf->r[i]);
		}

		int xi = (int) bf->Y[l + i * bf->n];
		bf->n_ijk[xi + s * bf->r[i]]++;
	}/* for l < n */

	return n_ij_count;
}
/* BDE_count_n_ijk() ======================================================= */

double
BDE_calc_bde(int * n_ijk, double fixed_alpha, int i, const int * parents, int q,
		int count, int * r)
{
	double score = 0.0;

	for (int j = 0; j < count; j++)
	{
		double aij = fixed_alpha * r[i];
		double nij = 0;

		for (int k = 0; k < r[i]; k++)
			nij += n_ijk[r[i] * j + k];

		double g_aij = lgamma(aij);
		double g_aijnij = lgamma(nij + aij);
		score += (g_aij - g_aijnij);

		for (int k = 0; k < r[i]; k++)
		{
			double g_a = lgamma(fixed_alpha + n_ijk[k + j * r[i]]);
			double g_b = lgamma(fixed_alpha);
			score += (g_a - g_b);
		}

	}/* for j < count */

	return score;
}
/* BDE_calc_bde() ========================================================== */

/* End of file.  Copyright (C) 2009, 2013 Yoshinori Tamada. */
