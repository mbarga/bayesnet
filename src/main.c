#include "main.h"
#include "library.h"
#include "readfile.h"
#include "bdescore.h"

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <string.h>
#include <glib.h> //TODO will this be supported on all machines?

#define PROG_TAG "MAIN"

int main(int argc, char *argv[])
{
    // TODO pass in some of these as arguments to the program?
    char input_file_name[1024] =
        "/home/mbarga/Dropbox/git/bayesnet/data/toy2.txt";

    int x_size = 0;
    int m_size = 3;

    GHashTable* hash = NULL;
    NODE *x = 	NULL;

    /* initialize logging */
    openlog(PROG_TAG, 0, LOG_USER);

    /* read in data samples from file and check for possible errors */
    int status = read_problem(input_file_name, &hash, &x, &x_size);
    g_hash_table_destroy(hash); //TODO remove

    if (status != 0)
    {
        errlog("main(): failed to read input file\n");
        return 1;
    }

    if (x == NULL)
    {
        errlog("main(): x was returned NULL\n");
        return 1;
    }

    if (m_size > x_size)
    {
        errlog("main(): |M| was larger than |X|\n");
        free_nodes(x, x_size);
        return 1;
    }

    /* allocate space for and initialize c and g matrices <- 0 */
    int *g = Malloc(int, x_size*x_size); // adjacency matrix of the network
    int *c = Malloc(int, x_size*x_size); // node selection count

    for (int i = 0; i < x_size; ++i)
    {
        for (int j = 0; j < x_size; ++j)
        {
            g[i * x_size + j] = 0;
            c[i * x_size + j] = 0;
        }
    }

    /*
     * calculate one to one local scores
     *  (i,j)th element is the local score() of graph gene_i -> gene_j
     */
    double *local_scores = NULL;
    one_to_one(&x, x_size, &local_scores);

    /*
     * repeated random sampling on subsets of X
     */
    int *m = Malloc(int, x_size); // array of permuted indices for the nodes X
    for (int i=0; i < x_size; ++i) m[i] = i;

    int n = 0;
    int u, v;
    while(n < MAX_SUBSETS)
    {
        randperm(m, m_size, x_size);

        /*
         * c_(u,v) <- c_(u,v) + 1 AND c_(v,u) <- c_(u,v) for all (X_u, X_v) in M
         */
        for (int i = 0; i < m_size; ++i)
        {
            for (int j = i+1; j < m_size; ++j)
            {
                u = m[i];
                v = m[j];
                ++c[u * x_size + v];
                ++c[v * x_size + u];
            }
        }

        /*
         * Estimate DAG on G_M using HC algorithm
         */
        //estimate_dag(&hash, m, m_size);
        ++n;
    }

    /*TODO
     * traverse the nodes and pick edges to add to the final graph G
     */

    /*
     * cleanup before exiting the algorithm
     */
    //print_matrix(c, x_size);
    //print_matrix(g,x_size);
    free(m);
    free(g);
    free(c);
    free(local_scores);
    free_nodes(x, x_size);

    syslog(LOG_INFO, "%s", "exiting cleanly\n");
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
void estimate_dag(NODE *x, int x_size, int *m, int m_size, int **g)
{
    double current_score = 0;
		double max_diff_a[3] = {0, 1, 0}; // {difference, action = 1, candidate parent index} 
		double max_diff_d[3] = {0, 2, 0}; //
		double max_diff_r[3] = {0, 3, 0}; //
    char improvement = 1;

    int i = 0;
    while (improvement == 1 && i < MAX_ITER)
    {
        improvement = 0;
        current_score = calc_bdescore();

        randperm(m, m_size, m_size); // rerandomize the indices to traverse

        // for all children u
				// only need to compute f() on the FAMILY SUBGRAPH containing child + parents
        for (int u = 0; u < m_size; ++u)
        {
            //TODO collapse all of this?
            //
            // test edge addition
						max_diff_a[0] = 0;
            for (int v = 0; v < m_size; ++v)
            {
                if (u == v) continue; // node cannot be its own parent

                if ( *g[u*x_size + v] == 0 )
                {
                    // calculate score diff = f(G + {x_v -> x_u)) - f(G)
										double diff = calc_bdescore() - current_score;

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
            for (int v = 0; v < m_size; ++v)
            {
                if (u == v) continue;

                if ( *g[u*x_size + v] == 1 )
                {
                    // calculate score diff = f(G / {x_v -> x_u)) - f(G)
										double diff = calc_bdescore() - current_score;

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
            for (int v = 0; v < m_size; ++v)
            {
                if (u == v) continue;

                if ( *g[u*x_size + v] == 1 )
                {
                    // calculate score diff = d1 + d2
										double d1 = calc_bdescore() - current_score; // f(G / {x_v -> x_u}) - f(G)
										double d2 = calc_bdescore() - current_score; // f(G + {x_u -> x_v}) - f(G)
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
										if (action[1] == 1) {*g[u*x_size + parent] = 1;}
										if (action[1] == 2) {*g[u*x_size + parent] = 0;}
										if (action[1] == 3) {*g[u*x_size + parent] = 0; *g[parent*x_size + u] = 1;}
						}

        }
        //TODO ensure no cycles?
				//TODO what if the score improvement is the same between two actions?
        //TODO for each selected directed edge in the G_M DAG, add the edge to the adjacency matrix g
        //       how to encode?

        ++i;
    }
}
