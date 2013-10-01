#include "search.h"
#include "main.h"
#include "score.h"
#include "library.h"

#include <stdio.h>
#include <stdlib.h>

void
randperm(int *, int, int);
double*
max(double *, double*, double*);

//TODO this should only traverse the group of highest scoring parent candidates
//TODO permute the candidate learning orders as well
/*
 * consideration (TRUE HC):
 * (1) find set of max socring parent candidates and hold them as a list
 * (2a) since no edges connected, try adding edges for each candidate and choose the
 *      best scoring
 * (2b) for each candidate try the other two operations that are not already employed
 *      (for empty graph, this will be reverse or add)
 * (3) when visiting the same node again, when visiting parent candidates, employ
 *      operations that are not already employed
 * (4) set maximum limit for parents when adding them
 *
 * another consideration:
 * (1) find set of maximum scoring parent candidates and add edges to all
 *      of thm for each node.
 * (2) select one gene at a time and then one parent edge at a time (in random order)
 *      and try reversing or adding
 *
 *      NOTES:
 *              when reversing edge, keep that parent candidate in the candidate history?
 *              (this candidate will now actually be a child)
 */
void
estimate_dag(double *X, NODE *Y, int n, int m, int max_parents,
        int max_candidates, int categories, int **G, int **C)
{
    double current_score = 0;
    double max_diff_a[3] =
        { 0, 1, 0 }; // {difference, action = 1, candidate parent index}
    double max_diff_d[3] =
        { 0, 2, 0 }; //
    double max_diff_r[3] =
        { 0, 3, 0 }; //
    char improvement = 1;

    /**
     * repeated random sampling on subsets of X
     */
    //TODO make static or something? need const length at compile time
    int candidates[n]; // array of permuted indices for the nodes X
    for (int i = 0; i < n; ++i)
        candidates[i] = i; // TODO what about not doing this every time?

    // initial set of random candidates
    randperm(candidates, max_candidates, n);

    /**
     * c_(u,v) <- c_(u,v) + 1 AND c_(v,u) <- c_(u,v)
     *      for all (X_u, X_v) in candidate set
     */
    int u, v;
    for (int i = 0; i < max_candidates; ++i)
        for (int j = i + 1; j < max_candidates; ++j)
        {
            u = candidates[i];
            v = candidates[j];
            C[u * n + v]++;
            C[v * n + u]++;
        }

    int i = 0;
    while ((improvement == 1) && (i < MAX_ITER))
    {
        improvement = 0;

        // TODO implement parent sets properly
        int parents[2];
        int parent_size = 2;

        void * buffer = bde_init(X, n, m, categories, max_candidates);
        current_score = get_score(buffer, parents, parent_size);

        // re-randomize only the candidate set indices
        randperm(candidates, max_candidates, max_candidates);

        /* for all children u */
        // only need to compute f() on the FAMILY SUBGRAPH containing child + parents
        for (int u = 0; u < max_candidates; ++u)
        {
            //TODO collapse all of this?
            //
            // test edge addition
            max_diff_a[0] = 0;
            for (int v = 0; v < max_candidates; ++v)
            {
                if (u == v)
                    continue; // node cannot be its own parent

                if (*G[u * n + v] == 0)
                {
                    // calculate gap = f(G + {x_v -> x_u)) - f(G)
                    double new_score = get_score(buffer, parents, parent_size);

                    double diff = new_score - current_score;

                    if (diff > max_diff_a[0])
                    {
                        // store the difference and the candidate parent index
                        max_diff_a[0] = diff;
                        max_diff_a[2] = v;
                        improvement = 1;
                    }
                }
            }

            // test edge deletion
            max_diff_d[0] = 0;
            for (int v = 0; v < max_candidates; ++v)
            {
                if (u == v)
                    continue;

                if (*G[u * n + v] == 1)
                {
                    // calculate gap = f(G / {x_v -> x_u)) - f(G)
                    double new_score = get_score(buffer, parents, parent_size);

                    double diff = new_score - current_score;

                    if (diff > max_diff_d[0])
                    {
                        // store the difference and the candidate parent index
                        max_diff_d[0] = diff;
                        max_diff_d[2] = v;
                        improvement = 1;
                    }
                }
            }

            // test edge reversal
            max_diff_r[0] = 0;
            for (int v = 0; v < max_candidates; ++v)
            {
                if (u == v)
                    continue;

                if (*G[u * n + v] == 1)
                {
                    // calculate gap = d1 + d2
                    double new_score_d1 = get_score(buffer, parents,
                            parent_size);
                    double d1 = new_score_d1 - current_score; // f(G / {x_v -> x_u}) - f(G)

                    double new_score_d2 = get_score(buffer, parents,
                            parent_size);
                    double d2 = new_score_d2 - current_score; // f(G + {x_u -> x_v}) - f(G)

                    double diff = d1 + d2;

                    if (diff > max_diff_r[0])
                    {
                        // store the difference and the candidate parent index
                        max_diff_r[0] = diff;
                        max_diff_r[2] = v;
                        improvement = 1;
                    }
                }
            }

            // choose max of three operations and apply the corresponding action to the graph
            double *action = max(max_diff_a, max_diff_d, max_diff_r);
            if (action != NULL)
            {
                int parent = action[2];
                if (action[1] == 1)
                {
                    *G[u * n + parent] = 1;
                }
                if (action[1] == 2)
                {
                    *G[u * n + parent] = 0;
                }
                if (action[1] == 3)
                {
                    *G[u * n + parent] = 0;
                    *G[parent * n + u] = 1;
                }
            }

        }
        //TODO ensure no cycles?
        //TODO what if the score improvement is the same between two actions?
        //TODO for each selected directed edge in the G_M DAG, add the edge to the adjacency matrix G
        //       how to encode?

        ++i;
    }
}

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
