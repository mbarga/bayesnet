#include "main.h"
#include "library.h"
#include "score.h"

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <math.h>
#include <glib.h>

#define PROG_TAG "LIBRARY"

/* Function
 * -------------------
 * calculate all 1-to-1 scores
 * TODO use array of nodes instead
 * TODO use bidirectional edge? or directed edge scores?
 * TODO move to score file?
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
            //TODO fix this
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

/*
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
*/
