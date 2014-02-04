//TODO this should only traverse the group of highest scoring parent candidates
//TODO tabu restart, etc.
//TODO ensure no cycles?
//TODO what if the score improvement is the same between two actions?
//       how to encode?
/*
 * repeated random sampling on subsets of X
 *
 * consideration (TRUE HC):
 * (1) find set of max scoring parent candidates and hold them as a list
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
 *      of them for each node.
 * (2) select one gene at a time and then one parent edge at a time (in random order)
 *      and try reversing or adding
 *
 *      NOTES:
 *              when reversing edge, keep that parent candidate in the candidate history?
 *              (this candidate will now actually be a child)
 */
#include "search.h"
#include "node.h"
#include "ran2.h"
#include "globals.h"
#include "score.h"
#include "BDE.h"
#include "util.h"
#include "main.h"

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <string.h>

// local constants
#define TOL 0.0001
#define MAX_ITER 40

int seed;

// local functions
void error_check(int *, int, NODE *, int[]);
void op_addition(int, int, void *, NODE, double, double *);
int op_deletion(int, int, void *, NODE, double, double *);
int op_reversal(int, int, int, void *, NODE, NODE, double, double *);
void random_permute(int *, int);
double * max_reduction(double *, double *, double *);
int apply_action(double[], NODE *, int *, int, int, int[], int);
void add_parent(int, NODE *);
int remove_parent(int, int, NODE *);

/* ##############################################################################*/
/** estimate_dag()
 * @brief
 *
 * @param p 
 * @param G
 * @param C
 */
/* ##############################################################################*/
//void estimate_dag(double *X, NODE *Y, int p, int n, int max_parents, int m, int r,
void estimate_dag(PARAMS parms, int *G)
{
	_CONFIG_ERROR status = E_SUCCESS;

// READ ONLY DATA
	// make local copies of data used from PARAMS due to multiple threads
	// running this routine
	const int p = parms.p;
	const int n = parms.n;
	const int m = parms.p; // FIXME clean this up later?
	const int r = parms.r;
	const int max_parents = parms.max_parents;
	double *X = parms.X;

// WRITE ONCE DATA
	NODE *Y = parms.Y;

// NEED PRIVATE COPIES OF THIS DATA 
	// {delta score, action, candidate parent index}
	double add_op[3] = { TOL, 1, -1 };
	double del_op[3] = { TOL, 2, -1 };
	double rev_op[3] = { TOL, 3, -1 };

	int u = -1; // child node index in global adjacency matrix G
	int j = -1; // child node index in local adjacency matrix G_M
	int v = -1; // parent node index in global adjacency matrix G
	int k = -1; // parent node index in local adjacency matrix G_M

	int no_improvement_cnt = 0; //


	// sample input data
	//printf("%d %d %d %d %d %d %f\n",p, n, m, r, max_parents, Y[0].num_parents, X[0]);

	// initialize candidate indices and shuffle them
	int candidates[p];
	for (int i = 0; i < p; ++i)
		candidates[i] = i;
	random_permute(candidates, p);

	// local adjacency matrix
	int *G_M = Calloc(int, m * m);
// END PRIVATE COPIES

	// initialize local adjacency matrix G_M
	//TODO is this right?
#pragma omp parallel shared(G_M, p) private(u,v)
#pragma omp for
	for (j = 0; j < m; ++j) {
		for (k = 0; k < m; ++k) {
			u = candidates[j];
			v = candidates[k];
			matrix(G_M, m, k, j) = matrix(G, p, v, u);
		}
	}
// end pragma omp parallel

	//void *buffer = score_init(X, p, n, r, max_parents);
	void *buffer = BDE_init(X, X, p, n, r, max_parents);

	// repeatedly apply HC on the candidate set until no no_improvement_cnt
	int i = 0;
	while ((no_improvement_cnt < max_parents) && (i < MAX_ITER)) {
		// chooses candidate at random from candidate set
		j = randinter(0, m);
		u = candidates[j];

		// initial score of candidate
		//printf("Y[%d] parents(cnt:%d) were: ", u, Y[u].num_parents);
		//for (int z = 0; z < max_parents; ++z)
		//    printf("%d, ", Y[u].parents[z]);
		//printf("\n");

		//TODO code this to a memory matrix to hold prev. calculated score
		double score = get_score(buffer, u, Y[u].parents, Y[u].num_parents);

#ifdef DEBUG
		printf("========================================\n");
		printf("CURRENT SCORE (of %d, #parents %d): %f\n",u, Y[u].num_parents, score);
#endif

		// reset all of the operation buffers
		add_op[0] = TOL;
		del_op[0] = TOL;
		rev_op[0] = TOL;
		add_op[2] = -1;
		del_op[2] = -1;
		rev_op[2] = -1;

		for (k = 0; k < m; ++k) {
			v = candidates[k];
			if (u == v)
				continue;

			// test hill climbing operations
			// if edge doesnt exist and parent set is not full
			if ((matrix(G_M,m,k,j)== 0) & (Y[u].num_parents < max_parents)) {
				op_addition(k, v, buffer, Y[u], score, add_op);
			}

			// if edge exists and parent set is not empty
			//TODO redundant?? edge exists and still check parent set?
			if ((matrix(G_M,m,k,j)== 1) & (Y[u].num_parents > 0)) {
				status = op_deletion(k, v, buffer, Y[u], score, del_op);

				if (Y[v].num_parents < max_parents) {
					status = op_reversal(j, k, v, buffer, Y[u], Y[v], score, rev_op);
				}
			}
		}

		// choose max of three operations and apply the corresponding action to the graph
		double *action = max_reduction(add_op, del_op, rev_op);

		// no improvement
		if (action == NULL) {
#ifdef DEBUG
			printf("NO ACTION TAKEN\n");
#endif
			no_improvement_cnt++;
		} else if (action[2] != -1) {
			status = apply_action(action, Y, G_M, max_parents, j, candidates,
					m);
			if (status != E_SUCCESS)
				util_errlog("ERROR IN apply_action()");

#ifdef DEBUG
			score = get_score(buffer, u, Y[u].parents, Y[u].num_parents);
			printf("SCORE AFTER: %f (diff: %f)\n", score, action[0]);
#endif

			no_improvement_cnt = 0;
		}


		error_check(G_M, m, Y, candidates);

		++i;
	} // while improvement

	//score_destroy_buff(buffer);

	/**
	 * reflect changes on G_M back into the global adjacency matrix G
	 */
#pragma omp parallel shared(G_M, p) private(u,v)
#pragma omp for
	for (int ii = 0; ii < m; ++ii)
		for (int jj = 0; jj < m; ++jj) {
			u = candidates[ii];
			v = candidates[jj];
			matrix(G,p,v,u)= matrix(G_M,m,jj,ii);
		}
	// end parallel
	
	free(G_M);
	//TODO return G_M and C_
}

void error_check(int *G_M, int m, NODE *Y, int candidates[])
{
	for (int ii = 0; ii < m; ++ii)
		//for(int jj = 0; jj < Y[candidates[ii]].num_parents; ++jj)
		for (int jj = 0; jj < m; ++jj)
			if (matrix(G_M,m,jj,ii)== 1) {
				NODE u = Y[candidates[ii]];
				//NODE v = Y[candidates[jj]];

				int g_idx = -1;
				for(int kk=0; kk < u.num_parents; ++kk)
					if(u.parents[kk] == candidates[jj]) {
						g_idx = kk;
						break;
					}

				if (g_idx == -1)
					fprintf(stderr,"AFTER:: PARENT %d(%d) OF %d(%d) DOESNT MATCH G_M\n",jj,candidates[jj],ii,candidates[ii]);
			}

	for (int ii = 0; ii < m; ++ii) {
		NODE u = Y[candidates[ii]];
		for (int jj = 0; jj < u.num_parents; ++jj) {
			int idx = -1;
			for (int kk = 0; kk < m; ++kk)
				if (candidates[kk] == u.parents[jj])
					idx = kk;

			if (idx != -1)
				if (matrix(G_M,m,idx,ii)!= 1) {
					fprintf(stderr,"AFTER :: G_M DOESNT MATCH PARENT (%d)%d OF (%d)%d\n",idx,candidates[idx],ii,candidates[ii]);
					printf("GM::::::::::::\n");
					for (int pp = 0; pp < m; ++pp) {
						for (int ll = 0; ll < m; ++ll)
							printf("%d, ",matrix(G_M,m,pp,ll));
						printf("\n");
					}
				}
		}
	}

}

/* ##############################################################################*/
/** op_addition
 * @brief calculate gap = f(G + {x_v -> x_u}) - f(G)
 *
 * @param k
 * @param buffer
 * @param u
 * @param current_score
 * @param max_diff_a[]
 */
/* ##############################################################################*/
void op_addition(int k, int v, void *buffer, NODE u, double current_score,
		double max_diff_a[])
{
	u.parents[u.num_parents] = v;
	double new_score = get_score(buffer, u.index, u.parents, u.num_parents + 1);
	u.parents[u.num_parents] = -1;

	double diff = current_score - new_score;

	if (diff > max_diff_a[0]) {
#ifdef DEBUG
		printf("diff of %d(%f) > than current diff %f\n",v,diff,max_diff_a[0]);
#endif

		max_diff_a[0] = diff;
		max_diff_a[2] = k;
	}

	return;
}
//TODO rename the int v variable
/* ##############################################################################*/
/** op_deletion()
 * @brief calculate gap = f(G / {x_v -> x_u)) - f(G)
 *
 * @param k
 * @param buffer
 * @param u
 * @param current_score
 * @param max_diff_d[]
 */
/* ##############################################################################*/
int op_deletion(int k, int v, void *buffer, NODE u, double current_score,
		double max_diff_d[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	int *temp = Calloc(int, (u.num_parents-1));
	memset(temp, -1, sizeof(temp[0]) * (u.num_parents-1));
	int temp_iter = 0;

	for (int i = 0; i < u.num_parents; ++i) {
		if (u.parents[i] != v)
			temp[temp_iter++] = u.parents[i];
		else
			status = E_DUPLICATE;
	}

	if (status != E_SUCCESS)
		return status;

	//    printf("REMOVE [ - %d] testing parents size %d: ",k, u.num_parents-1);
	//    for (int z = 0; z < u.num_parents-1; ++z)
	//        printf("%d, ", temp[z]);
	//    printf("\n");

	double new_score = get_score(buffer, u.index, temp, u.num_parents - 1);

	double diff = current_score - new_score;

	if (diff > max_diff_d[0]) {
		max_diff_d[0] = diff;
		max_diff_d[2] = k;
	}

	free(temp);
	return status;
}
//TODO rename the int nv variable
/* ##############################################################################*/
/** op_reversal()
 * @brief calculate gap = d1 + d2 where
 *  d1 = f(G / {x_v -> x_u}) - f(G)
 *  d2 = f(G + {x_u -> x_v}) - f(G)
 *
 * @param j
 * @param k
 * @param buffer
 * @param u
 * @param v
 * @param current_score
 * @param max_diff_r[]
 */
/* ##############################################################################*/
int op_reversal(int j, int k, int nv, void *buffer, NODE u, NODE v,
		double current_score, double max_diff_r[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	int *temp = Calloc(int, u.num_parents - 1);
	memset(temp, -1, sizeof(temp[0]) * (u.num_parents-1));
	int temp_iter = 0;

	// fill up temp buffer; check if parent already exists
	for (int i = 0; i < u.num_parents; ++i) {
		if (u.parents[i] != nv)
			temp[temp_iter++] = u.parents[i];
		else
			status = E_DUPLICATE;
	}

	if (status != E_SUCCESS)
		return status;

	//    printf("REV [ - %d] testing parents size %d: ",k,u.num_parents-1);
	//    for (int z = 0; z < u.num_parents-1; ++z)
	//        printf("%d, ", temp[z]);
	//    printf("\n");

	double new_score_d1 = get_score(buffer, u.index, temp, u.num_parents - 1);
	double d1 = current_score - new_score_d1;

	v.parents[v.num_parents] = j;

	//    printf("REV [ + %d] testing parents size %d: ",j,v.num_parents+1);
	//    for (int z = 0; z < u.num_parents-1; ++z)
	//        printf("%d, ", temp[z]);
	//    printf("\n");

	double new_score_d2 = get_score(buffer, v.index, v.parents, v.num_parents + 1);
	v.parents[v.num_parents] = -1;
	double d2 = current_score - new_score_d2;

	double diff = d1 + d2;

	if (diff > max_diff_r[0]) {
		max_diff_r[0] = diff;
		max_diff_r[2] = k;
	}

	free(temp);
	return status;
}

/* ##############################################################################*/
/** random_permute()
 * @brief randomize the first n values of the members of input array m
 *
 * @param m
 * @param n
 */
/* ##############################################################################*/
void random_permute(int *m, int n) {
	//srandinter((int)time(NULL));
	srandinter(seed++);

	if (m == NULL) {
		util_errlog("random_permute(): m was found NULL\n");
		return;
	}

	int j = 0, t = 0;
	for (int i = 0; i < n; ++i) {
		//j = rand() % (n - i) + i;
		j = randinter(0, n);
		t = m[j];
		m[j] = m[i];
		m[i] = t;
	}
}

/* ##############################################################################*/
/** max_reduction()
 * @brief
 *
 * @param a
 * @param b
 * @param c
 *
 * @return
 */
/* ##############################################################################*/
double * max_reduction(double *a, double *b, double *c) {
	double temp = TOL;
	char action = 0;

	if (a[0] > temp) {
		temp = a[0];
		action = 1;
	}

	if (b[0] > temp) {
		temp = b[0];
		action = 2;
	}

	if (c[0] > temp) {
		temp = c[0];
		action = 3;
	}

	if (action == 1)
		return a;
	else if (action == 2)
		return b;
	else if (action == 3)
		return c;

	return NULL;
}

/* ##############################################################################*/
/** apply_action()
 * @brief
 *
 * @param action[3]
 * @param Y
 * @param G_M
 * @param max_parents
 * @param j
 * @param candidates[]
 * @param m
 *
 * @return
 */
/* ##############################################################################*/
int apply_action(double action[3], NODE *Y, int *G_M, int max_parents, int child,
		int candidates[], int m) {
	_CONFIG_ERROR status = E_SUCCESS;

	int parent = action[2];
	int action_code = (int) action[1];
	int u = candidates[child];
	int v = candidates[parent];

	switch (action_code) {
	case 1: // edge addition
		add_parent(v, &Y[u]);
		matrix(G_M, m, parent, child) = 1;
#ifdef DEBUG
		printf("ADDING PARENT %d to CHILD %d\n",v,u);
		printf("CHILD parents after: ");
		for (int z = 0; z < max_parents; ++z)
			printf("%d, ", Y[u].parents[z]);
		printf("\n");
#endif
		break;

	case 2: // edge removal
		//              printf("CHILD parents before(%d): ", Y[u].num_parents);
		//              for (int z = 0; z < max_parents; ++z)
		//                  printf("%d, ", Y[u].parents[z]);
		//              printf("\n");
		status = remove_parent(v, max_parents, &Y[u]);
		if (status == E_SUCCESS)
			matrix(G_M, m, parent, child) = 0;

#ifdef DEBUG
		printf("REMOVING PARENT %d FROM CHILD %d ", v, u);
		printf("CHILD parents after(%d): ", Y[u].num_parents);
		for (int z = 0; z < max_parents; ++z)
			printf("%d, ", Y[u].parents[z]);
		printf("\n");
#endif

		break;

	case 3: // edge reversal
		//TODO what to do for reversal if target parent list is already full??
		if (Y[v].num_parents < max_parents) {
			/*
			   printf("CHILD parents before(%d): ", Y[u].num_parents);
			   for (int z = 0; z < max_parents; ++z)
			   printf("%d, ", Y[u].parents[z]);
			   printf("\n");
			   printf("TARGET parents before(%d): ", Y[v].num_parents);
			   for (int z = 0; z < max_parents; ++z)
			   printf("%d, ", Y[v].parents[z]);
			   printf("\n");
			   */

			// remove parent from u parent list
			status = remove_parent(v, max_parents, &Y[u]);
			if (status == E_SUCCESS)
				matrix(G_M, m, parent, child) = 0;

			// add u parent parent list
			add_parent(u, &Y[v]);
			matrix(G_M, m, child, parent) = 1;

			//            printf("CHILD parents after(%d): ", Y[u].num_parents);
			//            for (int z = 0; z < max_parents; ++z)
			//                printf("%d, ", Y[u].parents[z]);
			//            printf("\n");
			//            printf("TARGET parents before(%d): ", Y[v].num_parents);
			//            for (int z = 0; z < max_parents; ++z)
			//                printf("%d, ", Y[v].parents[z]);
			//            printf("\n");
		}

#ifdef DEBUG
		printf("REVERSING PARENT %d AND CHILD %d\n",v, u);
#endif

		break;
	}

	return status;
}

/* ##############################################################################*/
/** add_parent()
 * @brief add parent v to child node Y
 *
 * @param v parent node to add
 * @param Y child node to recieve parent
 */
/* ##############################################################################*/
void add_parent(int v, NODE *Y)
{
	// dont add a duplicate parent
	for (int i = 0; i < Y->num_parents; ++i)
		if (Y->parents[i] == v)
			return;

	int tail = Y->num_parents;
	Y->parents[tail] = v;
	Y->num_parents++;

	return;
}

//TODO better way?
/* ##############################################################################*/
/** remove_parent()
 * @brief remove parent v from child u
 *
 * @param u child node to remove parent from
 * @param v parent node to remove
 * @param max_parents maximum number of parents (for node parent array)
 * @param Y list of child nodes
 *
 * @return
 */
/* ##############################################################################*/
int remove_parent(int v, int max_parents, NODE *Y)
{
	_CONFIG_ERROR status = E_NOT_FOUND;

	for (int i = 0; i < max_parents; ++i)
	{
		// find target parent in the list
		if (Y->parents[i] == v)
		{   // shift all parents above target index down (overwrite target)
			for (int k = i; k < (max_parents - 1); ++k)
				Y->parents[k] = Y->parents[k + 1];

			// set the last element to empty and decrease count
			Y->parents[(max_parents - 1)] = -1;
			Y->num_parents--;
			status = E_SUCCESS;
			break;
		}
	}

	return status;
}

