#include "library.h"
#include "main.h"
#include "probability.h"
#include "readfile.h"
#include "score.h"
#include "search.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <syslog.h>
#include <string.h>

#define PROG_TAG "MAIN"

//TODO remove all redundant libs

// private function declarations
int
populate_nodes(double **, int, int, int);
void
init_edges(NODE *, int);
void
destroy_edges(NODE *, int);
void
print_dmatrix(double *c, int sizen, int sizem)
{
    int i, j;
    for (i = 0; i < sizen; ++i)
    {
        for (j = 0; j < sizem; ++j)
        {
            printf("%f, ", c[i * sizem + j]);
        }
        printf("\n");
    }
}
void
print_imatrix(int *c, int size)
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
    char input_file_name[1024] =
            "/home/mbarga/Dropbox/a_code/bayesnet/data/toy2.txt";

    double *X = NULL;           // n x m dataset
    int n = 0;                  // number of nodes
    int m = 0;                  // number of samples per node
    int max_candidates = 3;     // max candidate parent set size
    int categories = 3;         // number of discrete bins
    int max_parents = 5;        // maximum number of parents per node

    /* initialize logging */
    openlog(PROG_TAG, 0, LOG_USER);

    /* temporary measure to generate a test dataset */
    //int status = populate_nodes(&X, n, m, categories);
    /* read in data samples from file and check for possible errors */
    int status = read_problem(input_file_name, &X, &n, &m);

    assert(status == 0);
    assert(X != NULL);
    assert(max_candidates < n);

    print_dmatrix(X, n, m);

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

    int k = 0;
    while (k < MAX_SUBSETS)
    {
        estimate_dag(X, Y, n, m, max_parents, max_candidates, categories, &G,
                &C);
        k++;
    }

    /**TODO
     * traverse the nodes and pick edges to add to the final graph G
     */

    free(local_scores);
    free(G);
    free(C);
    free(X);
    destroy_edges(Y, n);

    //syslog(LOG_INFO, "%s", "exiting cleanly\n");
    closelog();
    return 0;
}

int
populate_nodes(double **X, int n, int m, int c)
{
    double *nodes = Malloc(double, n*m);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            nodes[i + j * m] = get_expression(NULL, 0);
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
