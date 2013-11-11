#include "score.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

//#define SIGN_DEBUG 1

typedef struct buffer
{
    int node_count;	// p
    int sample_count;	// n
    double *data;	// X
    int max_parents;	// m

    int prior;
    int multinomial;	// # categories TODO variable for each x_i?

    int *n_ij;		// used in BDe score
    int *n_ijk;		//

    double fixed_alpha;
} BUFFER;

// internal function declarations
int
count_nijk(BUFFER*, int*, int);
double
calc_bde(int*, double, int*, int, int, int);

/**Function
 * -------------------
 *
 *
 */
void*
bde_init(double *data, int node_count, int sample_count, int categories,
        int max_parents)
{
    BUFFER *buffer = Malloc(BUFFER, 1);

    buffer->data = data;
    buffer->node_count = node_count;
    buffer->sample_count = sample_count;
    buffer->multinomial = categories;
    buffer->max_parents = max_parents;

    buffer->prior = 0;
    buffer->fixed_alpha = 1.0;

    // |n_ij| = m * n
    buffer->n_ij = Malloc(int, (max_parents * sample_count));

    // |n_ijk| = n * r
    buffer->n_ijk = Malloc(int, (sample_count * categories));

    return buffer;
}

/**Function
 * -------------------
 *
 *
 */
void
bde_destroy_buff(void * buff)
{
    BUFFER * bf = (BUFFER *) buff;

    free(bf->n_ij);
    free(bf->n_ijk);

    free(buff);
}

/** Function
 * -------------------
 *
 *
 */
double
get_score(void *buffer, int *parents, int q)
{
    BUFFER *bf = (BUFFER *) buffer;

    /* likelihood */
    int count = count_nijk(bf, parents, q);

#ifdef SIGN_DEBUG
    fprintf(stderr, "BDE_score(): count = %d\n", count);
#endif

    double score = calc_bde(bf->n_ijk, bf->fixed_alpha, parents, q, count,
            bf->multinomial);

    /* prior */
    score = -2 * score;

    return score;
}

/** Function
 * -------------------
 *
 *
 */
int
count_nijk(BUFFER *buff, int *parents, int parents_count)
{
    int n_ij_count = 0;

    for (int i = 0; i < buff->sample_count; ++i)
    {
        int j[parents_count];
        for (int m = 0; m < parents_count; ++m)
        {
            //TODO rearrange indices??
            j[m] = buff->data[i + parents[m] * buff->sample_count];
        }
        //printf("DONE\n");

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
        }

        if (s == n_ij_count)
        {
            /* Not found.  Creats a new element */
            memcpy(buff->n_ij + s * buff->max_parents, j,
                    sizeof(int) * parents_count);
            n_ij_count++;

            bzero(buff->n_ijk + s * buff->multinomial,
                    sizeof(int) * buff->multinomial);
        }

        //TODO is this correct? buff->data[i + j(as in r[j]) * samplecount]
        int xi = (int) buff->data[i];
        //TODO remove
        //printf("reading index %d spots for n_ijk\n", xi+s*buff->multinomial);
        buff->n_ijk[xi + s * buff->multinomial]++;
    }

    return n_ij_count;
}

/** Function
 * -------------------
 *
 *
 */
double
calc_bde(int* n_ijk, double fixed_alpha, int* parents, int q, int count,
        int categories)
{
    double score = 0.0;

    // forall N_ij states:
    for (int j = 0; j < count; ++j)
    {
        double aij = fixed_alpha * categories;
        double nij = 0;

        // N_ij === SUM_{r_i}(N_ijk)
        for (int k = 0; k < categories; ++k)
            nij += n_ijk[categories * j + k];

        double g_aij = lgamma(aij);
        double g_aijnij = lgamma(nij + aij);
        score += (g_aij - g_aijnij);

        for (int k = 0; k < categories; ++k)
        {
            double g_a = lgamma(fixed_alpha + n_ijk[k + j * categories]);
            double g_b = lgamma(fixed_alpha);
            score += (g_a - g_b);
        }
    }/* for j < count */

    return score;
}
