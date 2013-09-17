#include "main.h"
#include "bdescore.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

/* Function
 * -------------------
 *
 *
 */
void *bde_init(double *data, int node_count, int sample_count, int categories,
		int max_parents)
{
	BUFFER *buff;
	buff = (BUFFER *) malloc(siseof(BUFFER));

	buff->prior = 0;
	buff->fixed_alpha = 1.0;
	buff->node_count = node_count;
	buff->sample_count = sample_count;
	buff->categories = categories;
	buff->max_parents = max_parents;
	buff->data = data;

	size_t mem;
	mem = sizeof(int) * buff->max_parents * buff->sample_count;
	buff->n_ij = (int *) malloc(mem);
	mem = sizeof(int) * buff->sample_count * categories;
	buff->n_ijk = (int *) malloc(mem);
}

/* Function
 * -------------------
 *
 *
 */
double get_score(void *buff, int j, int *parents, int q)
{
	BUFFER *bf = (BUFFER *) buff;

	/* likelihood */
	int count = count_n_ijk(bf, parents, q, j);

	double score = calc_bde(bf->n_ijk, bf->fixed_alpha, j, parents, q, count,
			bf->categories);

	/* prior */
	score = -2 * score;

	return score;
}

/* Function
 * -------------------
 *
 *
 */
int count_nijk(BUFFER *buff, int *parents, int parents_count)
{
	int n_ij_count = 0;

	for (int l = 0; l < buff->sample_count; ++l)
	{
		int j[parents_count];
		for (int m = 0; m < parents_count; ++m)
		{
			//TODO rearrange indices??
			j[m] = buff->data[l + parents[m] * buff->sample_count];
		}

		/* Find 'stat vector' from n_ij */
		int s = 0; /* stat vector position in n_ij */
		while (s < n_ij_count)
		{
			const int * vec = buff->n_ij + s * buff->max_parents;
			int flag2 = 0;
			for (int ii = 0; ii < parents_count; ++ii)
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
			memcpy(buff->n_ij + s * buff->max_parents, j,
					sizeof(int) * parents_count);
			n_ij_count++;

			bzero(buff->n_ijk + s * buff->categories,
					sizeof(int) * buff->categories);
		}

		/*TODO properly implement this
		 int xi = (int) bf->Y[l + i * bf->n];
		 bf->n_ijk[xi + s * bf->r[i]]++;
		 */
	}/* for l < n */

	return n_ij_count;
}

/* Function
 * -------------------
 *
 *
 */
double calc_bde(int *n_ijk, double fixed_alpha, int *parents, int q, int count,
		int categories)
{
	double score = 0.0;

	for (int j = 0; j < count; ++j)
	{
		double aij = fixed_alpha * categories;
		double nij = 0;
		for (int k = 0; k < categories; k++)
		{
			nij += n_ijk[categories * j + k];
		}

		double g_aij = lgamma(aij);
		double g_aijnij = lgamma(nij + aij);

		score += (g_aij - g_aijnij);

		for (int k = 0; k < categories; k++)
		{
			double g_a = lgamma(fixed_alpha + n_ijk[k + j * categories]);
			double g_b = lgamma(fixed_alpha);
			score += (g_a - g_b);
		}
	}/* for j < count */

	return score;
}
