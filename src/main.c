#include "library.h"
#include "main.h"
#include "node.h"
#include "probability.h"
#include "readfile.h"
#include "score.h"
#include "search.h"
#include "globals.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <syslog.h>
#include <string.h>
#include <time.h>

#define PROG_TAG "MAIN"

//TODO remove all redundant libs
//TODO pass in some of these as parms to the program?
#define NUM_REPETITIONS 10
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

const char filename[1024] = "./data/toy2.txt";
int seed;

// private functions
int
populate_nodes(double **, int, int, int);
void
init_edges(NODE **, int, int);
void
destroy_edges(NODE *, int);
void
init_parents(double *, int, NODE *);
void
print_dmatrix(double *c, int sizen, int sizem);
void
print_imatrix(int *c, int size);

int
main(int argc, char *argv[])
{
    //TODO X matrix is accessed in (i, j) X[i * n + m]
    double *X = NULL;           // n x m dataset
    NODE *Y = NULL;             // parent edges for each node
    int p = 0;                  // number of nodes
    int n = 0;                  // number of samples per node
    int r = 3;                  // number of discrete bins
    int m = 3;                  // max candidate parent set size
    int max_parents = 3;        // maximum number of parents per node

    double *local_scores = NULL;

    // initialize logging
    openlog(PROG_TAG, 0, LOG_USER);

    // read in data samples from file and check for possible errors
    int s = read_problem(filename, &X, &p, &n);
    //int status = populate_nodes(&X, n, m, categories);
    assert(s == 0);
    assert(X != NULL);
    assert(m < p);
    assert(max_parents < p);

    // allocate and initialize parent edges for each node
    init_edges(&Y, p, max_parents);
    assert(Y != NULL);

    /**calculate one to one local scores
     *  (i,j)th element is the local score() of graph gene_i -> gene_j
     */
    one_to_one(X, p, &local_scores, n, r);
    assert(local_scores != NULL);

    // initialize c and g matrices
    int *G = Malloc(int, p*p); // adjacency matrix of the network
    memset(G, 0, sizeof(G[0]) * p * p);
    int *C = Malloc(int, p*p); // node selection frequency matrix
    memset(C, 0, sizeof(C[0]) * p * p);

    // find top candidate parents for each node and assign them
    init_parents(local_scores, p, Y);

    // reflect parents of node into the adjacency matrix
    for (int i = 0; i < p; ++i)
        for (int j = 0; j < Y[i].num_parents; ++j)
            if (Y[i].parents[j] != -1)
                G[Y[i].parents[j] * p + i] = 1;

    //TODO no parents with high indeces?
    //print_imatrix(G,p);

    /**
     * init seed to the time ONLY ONCE at the start -
     *  seed is increased by one every time a random number is drawn
     *  if the seed is init'd to the time every time, sometimes the time
     *  does not update between random number draws
     *  TODO remove globals.h
     *  TODO think of a way to do this better (include random routines in main?)
     *          do seed++ after every call of randinter?
     */
    seed = (int) time(NULL); // init seed only once at the start

    /**TODO partition into random sets to run randomly
     * repeatedly run the HC search & score routine
     */
//    printf("starting HC routine");
    for (int i = 0; i < NUM_REPETITIONS; ++i)
    {
//        if ((i % 100000) == 0)
//        {
//            printf(".");
//            fflush(stdout);
//        }
        estimate_dag(X, Y, p, n, max_parents, m, r, G, C);
    }
//    printf("\n\n");

//TODO traverse the nodes and pick edges to add to the final graph G based on C

    printf("G was: \n");
    print_imatrix(G, p);

    for(int i = 0; i < p; ++i)
      for(int j = 0; j < Y[i].num_parents; ++j)
        if (G[Y[i].parents[j] * p + i] != 1) 
          printf("parents dont match adj matrix\n");

    free(X);
    destroy_edges(Y, p);
    free(local_scores);
    free(G);
    free(C);

//syslog(LOG_INFO, "%s", "exiting cleanly\n");
    closelog();
    return 0;
}

/* ##############################################################################*/
/**
 * @brief
 *
 * @param X
 * @param n
 * @param m
 * @param c
 *
 * @return
 */
/* ##############################################################################*/
int
populate_nodes(double **X, int n, int m, int c)
{
    double *nodes = Malloc(double, n*m);

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            nodes[i * m + j] = get_expression(NULL, 0);
            printf("%f, ", nodes[i * m + j]);
        }
        printf("\n");
    }
    *X = nodes;

    return 1;
}

/* ##############################################################################*/
/**
 * @brief
 *
 * @param Y
 * @param n
 * @param max_parents
 */
/* ##############################################################################*/
void
init_edges(NODE **Y, int n, int max_parents)
{
    NODE *p = Malloc(NODE, n);
    for (int i = 0; i < n; ++i)
    {
        p[i].parents = Malloc(int, max_parents);
        // initialize all parents to -1
        memset(p[i].parents, -1, sizeof(p[i].parents[0]) * max_parents);
        p[i].num_parents = 0;
    }
    *Y = p;
}

/* ##############################################################################*/
/**
 * @brief
 *
 * @param Y
 * @param n
 */
/* ##############################################################################*/
void
destroy_edges(NODE *Y, int n)
{
    for (int i = 0; i < n; ++i)
        free(Y[i].parents);
    //g_slist_free(Y[i].parents);

    free(Y);
}

/* ##############################################################################*/
/**
 * @brief
 *
 * @param local_scores
 * @param p
 * @param Y
 */
/* ##############################################################################*/
void
init_parents(double *local_scores, int p, NODE *Y)
{
    int max_parents = 3;
    /**TODO case of score tie? (there could be many ties)
     * assign best candidate parents to each node
     */
    for (int i = 0; i < p; ++i)
    {
        //TODO this should be dynamic for max_parents
        double top[max_parents];
        memset(top, 0, sizeof(top[0]) * max_parents);

        for (int j = 0; j < p; ++j)
        {
            double score = local_scores[i * p + j];
            int index = -1;

            if (score > top[0])
                index = 0;
            if (score > top[1])
                index = 1;
            if (score > top[2])
                index = 2;

            if (index > -1)
            {
                // shift down elements below where new number is to be inserted
                for (int k = 0; k <= (index - 1); ++k)
                {
                    top[k] = top[k+1];
                    Y[i].parents[k] = Y[i].parents[k+1];
                }
                // insert new max into proper index
                top[index] = score;
                Y[i].parents[index] = j;
                if (Y[i].num_parents < max_parents)
                    Y[i].num_parents++;
            }
        }
    }
}

void
print_dmatrix(double *c, int sizen, int sizem)
{
    for (int i = 0; i < sizen; ++i)
    {
        for (int j = 0; j < sizem; ++j)
            printf("%f, ", c[i * sizem + j]);
        printf("\n");
    }
}

void
print_imatrix(int *c, int size)
{
    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
            printf("%d, ", c[i * size + j]);
        printf("\n");
    }
}
