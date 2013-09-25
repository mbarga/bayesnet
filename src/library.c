#include "main.h"
#include "library.h"
#include "bdescore.h"

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <math.h>
#include <glib.h>

#define PROG_TAG "LIBRARY"

/* Function
 * -------------------
 *
 * TODO: return only array of size k OR optimize to shuffle only for k elements
 * Returns array of size N of random numbers in the range [1,N]
 *
 */
void
randperm(int *m, int k, int n)
{
    int i, j, t;

    /*
     for (i = 0; i < n; i++)
     {
     m[i] = i;
     }
     */
    if (m == NULL)
    {
        errlog("randperm(): m was found NULL\n");
        return;
    }

    for (i = 0; i < n; i++)
    {
        j = rand() % (n - i) + i;
        t = m[j];
        m[j] = m[i];
        m[i] = t;
    }
}

/* Function
 * -------------------
 *
 *
 */
void
scramble(int **indices, int size)
{
    int i, j, t;

    for (i = 0; i < size; ++i)
    {
        j = rand() % (size - i) + i;
        t = *indices[j];
        *indices[j] = *indices[i];
        *indices[i] = t;
    }
}

double *
max(double *a, double *b, double *c)
{
    double temp = 0;
    char action = 0;

    if (a[0] > temp)
    {
        temp = a[0];
        action = 1;
    }
    if (b[0] > temp)
    {
        temp = b[0];
        action = 2;
    }
    if (c[0] > temp)
    {
        temp = c[0];
        action = 3;
    }

    if (action == 1)
        return a;
    else if (action == 2)
        return b;
    else if (action == 3)
        return c;
    else
        return NULL;
}

/* Function
 * -------------------
 *
 * TODO use array of nodes instead
 * calculate all 1-to-1 scores
 *
 */
void
one_to_one(double *x, int x_size, double **local_scores, int sample_size,
        int categories)
{
    double *scores = Malloc(double, x_size*x_size);
    void *buff = bde_init(x, x_size, sample_size, categories, 1);

    // loop over each possible child-pair edge
    for (int u = 0; u < x_size; ++u)
    {
        for (int v = 0; v < x_size; ++v)
        {
            if (u == v)
                continue;

            /* score edge(u,v) */
            int parents[] =
                { v };
            double score = get_score(buff, parents, 1);
            scores[u * x_size + v] = score;
            //printf("local score (%d,%d) was: %f\n",u,parents[0],score);
        }
    }

    bde_destroy_buff(buff);

    *local_scores = scores;
    scores = NULL;
}

/* Function
 * -------------------
 *
 *TODO this cant be a negative # or throws off n_ijk count
 */
int
generate_local_probability()
{
    //TODO the resolution of the time seed is not high enough
    int seed = time(NULL);
    srand(seed);
    // {-1, 0, 1}
    return (int) (((rand() % 3) - 1));
}

/*
void
print_nodes(NODE *x, int x_size)
{
    int i;
    GSList * edge = NULL;
    NODE* test;

    for (i = 0; i < x_size; ++i)
    {
        printf("node :: %s with index %d, edges: ", x[i].name, x[i].index);
        edge = x[i].parents;
        while (edge != NULL)
        {
            test = edge->data;
            printf("%s, ", test->name);
            edge = edge->next;
        }
        printf("\n");
    }

    g_slist_free(edge);
    edge = NULL;
}
*/

void
errlog(char *s)
{
    fprintf(stderr, "%s", s);
    //syslog(LOG_INFO, "%s", "randperm(): m was found NULL\n");
}
