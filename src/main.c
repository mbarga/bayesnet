#include "util.h"
#include "main.h"
#include "node.h"
#include "readfile.h"
#include "score.h"
#include "BDE.h" // TODO REMOVE
#include "search.h"
#include "globals.h"

#include <stdlib.h>
#include <stdio.h>
#include <syslog.h>
#include <assert.h>
#include <string.h>
#include <time.h>

//TODO remove all redundant libs
#define NUM_REPETITIONS 40
#define PROG_TAG "BAYES_HC"

const char filename[1024] = "./data/1000sim";
int seed;

// private functions
static void finalize_params(PARAMS params);
static int init_parents(NODE **, const int, const int);
static void destroy_parents(NODE *, const int);
//TODO make double * const * const * ...
static void populate_candidate_parents(NODE *, const int, double *);
static int one_to_one(const PARAMS, double **);
//static int populate_nodes(double **, int, int, int);

int main(int argc, char *argv[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	// initialize logging
	//setlogmask (LOG_UPTO (LOG_NOTICE));
	openlog(PROG_TAG, 0, LOG_USER);

	//TODO set these in arguments?
	PARAMS params = {NULL, NULL, 0, 0, 3, 0, 7}; // X, Y, p, n, r, m, max_parents

	// read in data samples from file and record size of data
	status = read_problem(filename, &params.X, &params.p, &params.n);
	params.m = params.p; // TODO change this later?
	assert(status == E_SUCCESS);
	assert(params.max_parents < params.p);
	assert(params.m <= params.p);

#ifdef DEBUG
	printf("%d nodes, %d samples, %d categories, %d max candidate parents, %d max parents\n", params.p, params.n, params.r, params.m, params.max_parents);
#endif

	// allocate and initialize parent edges for each node
	status = init_parents(&params.Y, params.p, params.max_parents);
	assert(status == E_SUCCESS);

	//calculate 1x1 directed edge scores
	//  (i,j)th element is the local score() of graph gene_j -> gene_i
	double *local_scores = NULL;
	status = one_to_one(params, &local_scores);
	assert(status == E_SUCCESS);
	//util_print_dmatrix(local_scores, params.p, params.p);

	// find top candidate parents for each node and assign them
	// TODO update this to choose lowest scoring parents
	//populate_candidate_parents(params.Y, params.p, local_scores);

	// adjacency matrix of the network
	int *G = Calloc(int, params.p * params.p);

	// reflect parents of nodes into the adjacency matrix
	for (int i = 0; i < params.p; ++i)
		for (int j = 0; j < params.Y[i].num_parents; ++j)
			if (params.Y[i].parents[j] != -1)
				matrix(G, params.p, params.Y[i].parents[j], i) = 1;

	//printf("G matrix after initial candidate parents: \n");
	//util_print_imatrix(G, params.p);
	//printf("\n");

#ifdef DEBUG
	//void *buff = score_init(params.X, params.p, params.n, params.r, params.max_parents);
	void *buff = BDE_init(params.X, params.X, params.p, params.n, params.r, params.max_parents);
	util_print_score_table(buff);
	//score_destroy_buff(buff);
#endif

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

	//TODO partition into random sets to run randomly
	// repeatedly run the HC search & score routine
	printf("\n## Starting HC routine\n");
	for (int i = 0; i < NUM_REPETITIONS; ++i) {
		/*
		if ((i % 10) == 0) {
			printf(".");
			fflush(stdout);
		}
		*/
		estimate_dag(params, G);
	}
	//printf("\n---- FINAL OUTPUT GRAPH -> G[] ----\n");
	util_print_imatrix(G, params.p);

	for (int i = 0; i < params.p; ++i)
		for (int j = 0; j < params.Y[i].num_parents; ++j)
			if (G[params.Y[i].parents[j] * params.p + i] != 1)
				util_errlog("PARENTS DONT MATCH ADJ MATRIX");

	free(G);
	free(local_scores);
	finalize_params(params);

	closelog();
	printf("\n## Exiting Normally\n");
	return status;
}

/*
##############################################################################*/
/**
 * @brief frees up all memory allocated for problem data defined in PARAMS p
 *
 * @param p PARAMS struct target
 */
/*
##############################################################################*/
void finalize_params(PARAMS p)
{
	free(p.X);
	destroy_parents(p.Y, p.p);

	return;
}

/* ##############################################################################*/
/**
 * @brief allocates space for Y; allocates arrays of size max_parents for each 
 * node in Y and inits values to -1
 *
 * @param Y target array of NODEs to have parent arrays initialized
 * @param p total number of NODEs in Y
 * @param max_parents size of parent array to allocate
 */
/* ##############################################################################*/
int init_parents(NODE **Y, const int p, const int max_parents)
{
	_CONFIG_ERROR status = E_SUCCESS;
	NODE *a = Calloc(NODE, p);

	for (int i = 0; i < p; ++i) {
		a[i].parents = Calloc(int, max_parents);

		// initialize all parents to -1
		memset(a[i].parents, -1, sizeof(a[i].parents[0]) * max_parents);
		a[i].num_parents = 0;

		// remember index for the node
		a[i].index = i;
	}

	*Y = a;

	if (a == NULL)
		status = E_NULL_POINTER;

	return status;
}

/* ##############################################################################*/
/**
 * @brief frees memory allocated for parent arrays of each NODE
 *
 * @param Y target array of NODEs to have parent elements destroyed
 * @param p total number of NODEs in Y
 */
/* ##############################################################################*/
void destroy_parents(NODE *Y, const int p)
{
	for (int i = 0; i < p; ++i)
		free(Y[i].parents);

	free(Y);
}

/* ##############################################################################*/
/**
 * @brief initialize parent candidate sets based on highest 1x1 scores (assumes parents
 * are empty)
 *
 * @param p struct containing node data
 * @param local_scores local score matrix (non-NULL)
 */
/* ##############################################################################*/
void populate_candidate_parents(NODE *Y, const int p, double *local_scores)
{
	// holds highest scoring candidate scores for a node (descending order)
	//double *min_score_buff = Calloc(double, p.max_parents);
	double score = INF;
	double min_score = INF;
	int min_parent = -1;
	//int parent_slot = -1;

	for (int i = 0; i < p; ++i) {

		min_parent = -1;
		min_score = INF;

		printf("new line\n");
		for (int j = 0; j < p; ++j) {
			if (j == i)
				continue;

			score = matrix(local_scores, p, j, i);
			printf("score was %f\n", score);
			if (score < min_score) {
				min_score = score;
				min_parent = j;
			}

			/* GATHER LOWEST SCORING 1-to-1 NODES 
			// (i,j) represents score of i for parent set {j}
			score = matrix(local_scores, p.p, i, j);

			parent_slot = -1;

			//TODO case of a tie?
			// see if score is higher than previous j 
			for (int k = (p.max_parents-1); k >= 0; --k)
				if (score < min_score_buff[k])
					parent_slot = k;

			// if a new max score was found
			if (parent_slot > -1) {
				// shift up elements above where new number is to be inserted
				for (int k = (p.max_parents-1); k > parent_slot; --k) {
					min_score_buff[k] = min_score_buff[k-1];
					p.Y[i].parents[k] = p.Y[i].parents[k-1];
				}

				// insert new max score at parent_slot
				min_score_buff[parent_slot] = score;
				p.Y[i].parents[parent_slot] = j;

				// if parent set full, do not increase count
				if (p.Y[i].num_parents < p.max_parents)
					p.Y[i].num_parents++;
			}
			***********************************/
		}

		if (min_parent > -1) {
			printf("ADDING CANDIDATE: ");
			Y->parents[0] = min_parent;
			printf("Y[0] is: %d\n",Y->parents[0]);
		}

		//reset buff to all 0's
		//memset(min_score_buff, 0, sizeof(min_score_buff[0]) * p.max_parents);
	}

	//TODO invalid size for free() below?
	//free(min_score_buff);
	return;
}

/* ##############################################################################*/
/**
 * @brief calculate all 1-to-1 (DIRECTED EDGE) scores
 *
 * @param p struct containing node data 
 * @param local_scores 1x1 scores for network (typically NULL)
 */
/* ##############################################################################*/
int one_to_one(const PARAMS p, double **local_scores)
{
	_CONFIG_ERROR status = E_SUCCESS;
	double *scores = Calloc(double, p.p * p.p);

//#pragma omp parallel shared(p.X, p.p, p.n, p.r, scores)
	{
		//void *buff = score_init(p.X, p.p, p.n, p.r, 1);
		void *buff = BDE_init(p.X, p.X, p.p, p.n, p.r, 1); //TODO remove tamada
		//double t_start = omp_get_wtime();

//#pragma omp for
		// loop over each possible child-pair edge
		for (int u = 0; u < p.p; ++u) {
			for (int v = 0; v < p.p; ++v) {
				if (u == v)
					continue;

				double score = get_score(buff, u, &v, 1);
				//double score = BDE_score(buff, 0, &v, 1); //tamada

				// mark parent scores {v} for each row u
				matrix(scores, p.p, u, v)= score;
			}
		}

		//printf("local_score :: thread%d time on wall %f\n", omp_get_thread_num(), omp_get_wtime() - t_start);
		//score_destroy_buff(buff);
	} // end pragma omp parallel

	*local_scores = scores;

	if (scores == NULL)
		status = E_NULL_POINTER;

	return status;
}

//TODO REMOVE
/*
   int populate_nodes(double **X, int n, int m, int c)
   {
   double *nodes = Calloc(double, n*m);
   _CONFIG_ERROR status = E_SUCCESS;

   for (int i = 0; i < n; ++i) {
   for (int j = 0; j < m; ++j) {
   nodes[i * m + j] = get_expression(NULL, 0);
   printf("%f, ", nodes[i * m + j]);
   }
   printf("\n");
   }

 *X = nodes;

 return status;
 }
 */
