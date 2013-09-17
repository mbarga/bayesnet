#include "main.h"
#include "bdescore.h"

/* Function
 * -------------------
 *
 *
 */
//double score(NODE *x, int *parents, int m_size)
double calc_bdescore() 
{
    double score = rand() % 100;
    BDEBuff * bf = (BDEBuff *) buff;

    /* likelihood */`
    int count = BDE_count_n_ijk(bf, parents, q, j);
#ifdef SIGN_DEBUG
    fprintf(SIGN_ERR, "BDE_score(): count = %d\n", count);
#endif
    double score = BDE_calc_bde(bf->n_ijk, bf->fixed_alpha, j, parents, q,
                                count, bf->r);

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

// (Buff *bf, "", "")
int count_nijk(const int *parents, int q, int i)
{

    int n_ij_count = 0;

    for (int l = 0; l < bf->n; l++)
    {
        /* j[m] : parent stat vector */
        int j[q];
        for (int m = 0; m < q; m++)
            j[m] = bf->X[l + parents[m] * bf->n];

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

double BDE_calc_bde(int * n_ijk, double fixed_alpha, int i, const int * parents,
                    int q, int count, int * r)
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
