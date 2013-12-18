#include "util.h"
#include "main.h"
#include "node.h"
#include "readfile.h"
//#include "score.h"
#include "search.h"
#include "globals.h"
#include "BDE.h" // TODO REMOVE

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <time.h>

//TODO remove all redundant libs
#define NUM_REPETITIONS 10

const char filename[1024] = "./5sim";
int seed;

// private functions
static void finalize_params(PARAMS params);
static int init_parents(NODE **, const int, const int);
static void destroy_parents(NODE *, const int);
//TODO make double * const * const * ...
static void populate_candidate_parents(PARAMS, double *);
static int one_to_one(const PARAMS, double **);
//static int populate_nodes(double **, int, int, int);

int main(int argc, char *argv[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	//TODO set these in arguments?
	PARAMS params = { NULL, NULL, 0, 0, 3, 0, 3 }; // X, Y, p, n, r, m, max_parents
	double *local_scores = NULL;

	// read in data samples from file and record size of data
	status = read_problem(filename, &params.X, &params.p, &params.n);
	// TODO change this later?
	params.m = params.p;
	assert(status == E_SUCCESS);
	assert(params.max_parents < params.p);
	assert(params.m <= params.p);

	printf("%d nodes, %d samples, %d categories, %d max candidate parents, %d max parents\n", params.p, params.n, params.r, params.m, params.max_parents);
	/*
	for (int i = 0; i < params.p; ++i) {
		for (int j = 0; j < params.n; ++j)
			//printf("%d, ", (int)params.X[i * params.n + j]);
			printf("%d, ", (int)params.X[j * params.p + i]);
		printf("\n\n");
	}
	*/
	for (int i = 0; i < params.n*params.p; ++i) {
			//printf("%d, ", (int)params.X[i * params.n + j]);
		//	printf("%d, ", (int)params.X[i]);
	}

	// allocate and initialize parent edges for each node
	status = init_parents(&params.Y, params.p, params.max_parents);
	assert(status == E_SUCCESS);

	//calculate 1x1 directed edge scores
	//  (i,j)th element is the local score() of graph gene_j -> gene_i
	status = one_to_one(params, &local_scores);
	assert(status == E_SUCCESS);
	//util_print_dmatrix(local_scores, params.p, params.p);

	// find top candidate parents for each node and assign them
	// TODO update this to choose lowest scoring parents
	//populate_candidate_parents(params, local_scores);

	// adjacency matrix of the network
	int *G = Calloc(int, params.p * params.p);
	// node selection frequency matrix
	int *C = Calloc(int, params.p * params.p);

	// reflect parents of nodes into the adjacency matrixS
	for (int i = 0; i < params.p; ++i)
		for (int j = 0; j < params.Y[i].num_parents; ++j)
			if (params.Y[i].parents[j] != -1)
				matrix(G, params.p, params.Y[i].parents[j], i) = 1;

	//printf("G matrix after initial candidate parents: \n");
	//util_print_imatrix(G, params.p);
	//printf("\n\n");

	//void *buff = bde_init(params.X, params.p, params.n, params.r, 1);
	void *buff = BDE_init(params.X, params.X, params.p, params.n, params.r, params.max_parents);

	int parents[5] = {0,1,2,3,4};
	int four[4] = {0,0,0,0};
	int three[3] = {0,0,0};
	int two[2] = {0,0};
	double score;
	// repeat 5 times;
	for(int i = 0; i < 5; ++i) {
		printf("child: %d =======\n",parents[0]);

		// four set
		//printf("|4|: ===\n");
		for (int j = 0; j < 4; ++j)
			four[j] = parents[j+1];
		printf("parent set %d, {%d, %d, %d, %d} -> ",parents[0], four[0],four[1],four[2],four[3]);
		score = get_score(buff, parents[0], four, 4);
		printf("%f\n",score);

		// three set
		//printf("|3|: ===\n");
		int tmp[4];
		for (int p=0; p<4; ++p) tmp[p] = parents[p+1];
		for (int j = 0; j < 4; ++j) {
			for (int k=0; k<3; ++k) three[k] = tmp[k];
			score = get_score(buff, parents[0], three, 3);
			printf("parent set %d, {%d, %d, %d} -> ",parents[0], three[0],three[1],three[2]);
			printf("%f\n",score);

			int temp = tmp[3];
			for(int m=3; m>0; --m) tmp[m] = tmp[m-1];
			tmp[0] = temp;
		}

		// two set
		//printf("|2|: ===\n");
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[0],tmp[1]);
		two[0] = tmp[0];
		two[1] = tmp[1];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		//printf("parent set {%d, %d} -> ",tmp[1],tmp[0]);
		//two[0] = tmp[1];
		//two[1] = tmp[0];
		//score = get_score(buff, two, 2);
		//printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[1],tmp[2]);
		two[0] = tmp[1];
		two[1] = tmp[2];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		//printf("parent set {%d, %d} -> ",tmp[2],tmp[1]);
		//two[0] = tmp[2];
		//two[1] = tmp[1];
		//score = get_score(buff, two, 2);
		//printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[2],tmp[3]);
		two[0] = tmp[2];
		two[1] = tmp[3];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[0],tmp[2]);
		two[0] = tmp[0];
		two[1] = tmp[2];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[1],tmp[3]);
		two[0] = tmp[1];
		two[1] = tmp[3];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);
		printf("parent set %d, {%d, %d} -> ",parents[0],tmp[0],tmp[3]);
		two[0] = tmp[0];
		two[1] = tmp[3];
		score = get_score(buff, parents[0], two, 2);
		printf("%f\n",score);

		// one set
		//printf("|1|: ===\n");
		int test[] = { 0 };
		for (int j = 0; j<4; ++j) {
			//printf("parent set {%d} -> ",tmp[j]);
			//score = get_score(buff, &tmp[j], 1);
			test[0] = tmp[j];
			printf("parent set %d, {%d} -> ",parents[0], test[0]);
			score = get_score(buff, parents[0], test, 1);
			printf("%f\n",score);
		}

		//printf("|0|: ===\n");
		test[0] = i;
		printf("parent set %d, {} -> ", parents[0]);
		score = get_score(buff, parents[0], test, 0);
		printf("%f\n",score);

		//rotate the parent array
		int temp = parents[4];
		for(int j = 4; j > 0; --j) parents[j] = parents[j-1];
		parents[0] = temp;
	}


	//bde_destroy_buff(buff);

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
	//printf("##\nStarting HC routine");
	for (int i = 0; i < NUM_REPETITIONS; ++i) {
		/*
		   if ((i % 1000) == 0) {
		   printf(".");
		   fflush(stdout);
		   }
		   */
		//estimate_dag(X, Y, p, n, max_parents, m, r, G, C);
		//TODO how does parms struct impose parallellism
		estimate_dag(params, G, C);
	}
	printf("\n\n");

	//TODO traverse the nodes and pick edges to add to the final graph G based on C

	printf("G was: \n");
	util_print_imatrix(G, params.p);

	for (int i = 0; i < params.p; ++i)
		for (int j = 0; j < params.Y[i].num_parents; ++j)
			if (G[params.Y[i].parents[j] * params.p + i] != 1)
				printf("parents dont match adj matrix\n");

	//free(C);
	free(G);
	//free(local_scores);
	finalize_params(params);

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
void populate_candidate_parents(PARAMS p, double *local_scores)
{
	// holds highest scoring candidate scores for a node (descending order)
	double *max_score_buff = Calloc(double, p.max_parents);
	double score = 0;
	int parent_slot = -1;

	for (int i = 0; i < p.p; ++i) {
		for (int j = 0; j < p.p; ++j) {
			if (j == i)
				continue;

			// (i,j) represents score of i for parent set {j}
			score = matrix(local_scores, p.p, i, j);

			parent_slot = -1;

			//TODO case of a tie?
			// see if score is higher than previous j 
			for (int k = (p.max_parents-1); k >= 0; --k)
				if (score > max_score_buff[k])
					parent_slot = k;

			// if a new max score was found
			if (parent_slot > -1) {
				// shift up elements above where new number is to be inserted
				for (int k = (p.max_parents-1); k > parent_slot; --k) {
					max_score_buff[k] = max_score_buff[k-1];
					p.Y[i].parents[k] = p.Y[i].parents[k-1];
				}

				// insert new max score at parent_slot
				max_score_buff[parent_slot] = score;
				p.Y[i].parents[parent_slot] = j;

				// if parent set full, do not increase count
				if (p.Y[i].num_parents < p.max_parents)
					p.Y[i].num_parents++;
			}
		}

		//reset buff to all 0's
		memset(max_score_buff, 0, sizeof(max_score_buff[0]) * p.max_parents);
	}

	//TODO invalid size for free() below?
	free(max_score_buff);
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

#pragma omp parallel shared(p.X, p.p, p.n, p.r, scores)
	{
		//void *buff = bde_init(p.X, p.p, p.n, p.r, 1);
		void *buff = BDE_init(p.X, p.X, p.p, p.n, p.r, 1); //TODO remove tamada
		//double t_start = omp_get_wtime();

#pragma omp for
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
		//bde_destroy_buff(buff);
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
