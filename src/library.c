#include "library.h"
#include "score.h"

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define OMP_DEBUG = 1

/* ##############################################################################*/
/**
 * @brief  calculate all 1-to-1 (DIRECTED EDGE) scores
 *
 * @param X
 * @param p
 * @param local_scores
 * @param n
 * @param r
 */
/* ##############################################################################*/
void
one_to_one(double   *X, 
           int       p, 
           double  **local_scores, 
           int       n, 
           int       r)
{
    double *scores = Malloc(double, p*p);
    memset(scores, 0, sizeof(scores[0]) * p * p);

#pragma omp parallel shared(X, p, n, r, scores)
    {
        void *buff = bde_init(X, p, n, r, 1);
        //double t_start = omp_get_wtime();

#pragma omp for
        // loop over each possible child-pair edge
        for (int u = 0; u < p; ++u)
            for (int v = 0; v < p; ++v)
            {
                if (u == v)
                    continue;

                double score = get_score(buff, &v, 1);
                scores[v * p + u] = score;
            }
        //printf("local_score :: thread%d time on wall %f\n", omp_get_thread_num(), omp_get_wtime() - t_start);
        bde_destroy_buff(buff);
    } // end pragma omp parallel

    *local_scores = scores;
    scores = NULL;
}

/* ##############################################################################*/
/**
 * @brief 
 *
 * @param s
 */
/* ##############################################################################*/
void
errlog(char *s)
{
    fprintf(stderr, "%s", s);
    //syslog(LOG_INFO, "%s", "randperm(): m was found NULL\n");
}
