#include "main.h"
#include "library.h"
#include "readfile.h"
#include "bdescore.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <syslog.h>
#include <string.h>

#define PROG_TAG "MAIN"

// private function declarations
int
populate_nodes(double **, int, int, int);
void
estimate_dag(double *X, NODE *, int, int, int, int *, int, int, int **);
void
init_edges(NODE *, int);
void
destroy_edges(NODE *, int);

void
print_matrix(int *c, int size)
{
    int i, j;
    for (i = 0; i < size; ++i)
    {
        for (j = 0; j < size; ++j)
        {
            printf("%d, ", c[i * size + j]);
        }
        printf("\n");
    }
}

int
main(int argc, char *argv[])
{
    // TODO pass in some of these as arguments to the program?
    //char input_file_name[1024] =
    //        "/home/mbarga/Dropbox/git/bayesnet/data/toy2.txt";

    double *X = NULL;           // n x m dataset
    //NODE *Y = NULL;             // contains parent sets for X
    int n = 10;                 // number of nodes
    int m = 2;                  // number of samples per node
    int max_candidates = 3;     // max candidate parent set size
    int categories = 3;         // number of discrete bins
    int max_parents = 5;        // maximum number of parents per node

    /* initialize logging */
    openlog(PROG_TAG, 0, LOG_USER);

    /* temporary measure to generate a test dataset */
    int status = populate_nodes(&X, n, m, categories);
    /* read in data samples from file and check for possible errors */
    //int status = read_problem(input_file_name, &hash, &x, &x_size);
    
    assert(status != 0);
    assert(X != NULL);
    assert(max_candidates < n);

    /* initialize parent sets for each node */
    NODE Y[n];
    init_edges(Y, n);

    /* allocate space for and initialize c and g matrices <- 0 */
    int *G = Malloc(int, n*n); // adjacency matrix of the network
    memset(G, 0, n * n);
    int *C = Malloc(int, n*n); // node selection count
    memset(C, 0, n * n);

    /**
     * calculate one to one local scores
     *  (i,j)th element is the local score() of graph gene_i -> gene_j
     */
    double *local_scores = NULL;
    one_to_one(X, n, &local_scores, m, categories);
    assert(local_scores != NULL);

    /**
     * repeated random sampling on subsets of X
     */
    int *candidates = Malloc(int, n); // array of permuted indices for the nodes X
    for (int i = 0; i < n; ++i)
        candidates[i] = i;

    int u, v;
    int k = 0;
    while (k < MAX_SUBSETS)
    {
        randperm(candidates, max_candidates, n);

        /**
         * c_(u,v) <- c_(u,v) + 1 AND c_(v,u) <- c_(u,v) for all (X_u, X_v) in M
         */
        for (int i = 0; i < max_candidates; ++i)
            for (int j = i + 1; j < max_candidates; ++j)
            {
                u = candidates[i];
                v = candidates[j];
                C[u * n + v]++;
                C[v * n + u]++;
            }

        estimate_dag(X, Y, n, m, max_parents, candidates, max_candidates, categories, &G);
        k++;
    }

    /**TODO
     * traverse the nodes and pick edges to add to the final graph G
     */

    free(local_scores);
    free(candidates);
    free(G);
    free(C);
    free(X);
    destroy_edges(Y, n);

    //syslog(LOG_INFO, "%s", "exiting cleanly\n");
    closelog();
    return 0;
}

//TODO this should only traverse the group of highest scoring parent candidates
//TODO permute the candidate learning orders as well
/*
 * consideration (TRUE HC):
 * (1) find set of max socring parent candidates and hold them as a list
 * (2a) since no edges connected, try adding edges for each candidate and choose the
 * 	best scoring
 * (2b) for each candidate try the other two operations that are not already employed
 * 	(for empty graph, this will be reverse or add)
 * (3) when visiting the same node again, when visiting parent candidates, employ
 * 	operations that are not already employed
 * (4) set maximum limit for parents when adding them
 *
 * another consideration:
 * (1) find set of maximum scoring parent candidates and add edges to all
 * 	of thm for each node.
 * (2) select one gene at a time and then one parent edge at a time (in random order)
 * 	and try reversing or adding
 *
 * 	NOTES:
 * 		when reversing edge, keep that parent candidate in the candidate history?
 * 		(this candidate will now actually be a child)
 */
void
estimate_dag(double *X, NODE *Y, int n, int m, int max_parents, int *candidates,
        int max_candidates, int categories, int **G)
{
    double current_score = 0;
    double max_diff_a[3] =
        { 0, 1, 0 }; // {difference, action = 1, candidate parent index}
    double max_diff_d[3] =
        { 0, 2, 0 }; //
    double max_diff_r[3] =
        { 0, 3, 0 }; //
    char improvement = 1;

    int i = 0;
    while (improvement == 1 && i < MAX_ITER)
    {
        improvement = 0;

        // TODO implement parent sets properly
        int parents[2];
        int parent_size = 2;
        void * buffer = bde_init(X, n, m, categories, max_candidates);
        current_score = get_score(buffer, parents, parent_size);

        randperm(candidates, max_candidates, max_candidates); // rerandomize the indices to traverse

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
                    double new_score_d1 = get_score(buffer, parents, parent_size);
                    double d1 = new_score_d1 - current_score; // f(G / {x_v -> x_u}) - f(G)

                    double new_score_d2 = get_score(buffer, parents, parent_size);
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

int
populate_nodes(double **X, int n, int m, int c)
{
    double *nodes = Malloc(double, n*m);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            nodes[i + j * m] = 1; //generate_local_probability();
            printf("%f, ", nodes[i + j * m]);
        }
        printf("\n");
    }

    *X = nodes;

    return 1;
}

void
init_edges(NODE *Y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        //&Y[i] = Malloc(NODE, 1);
        Y[i].parents = NULL;
    }
}

void
destroy_edges(NODE *Y, int n)
{
    for (int i = 0; i < n; ++i)
    {
        g_slist_free(Y[i].parents);
        //free(Y[i]);
    }
}
