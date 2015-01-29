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
#include <assert.h>
#include <syslog.h>
#include <sys/time.h>
#include <time.h>

// local constants
#define TOL		0.0001
#define MAX_ITER	100
#define IMP_LIMIT	40
#define NTHREADS	1

struct timeval t1, t2;
int seed;
clock_t sbegin;
clock_t send;
double wbegin;
double wend;
double stime_spent = 0;
double wtime_spent = 0;

// local functions
void	error_check(int *, int, NODE *, int[]);
int	op_addition(int, int, void *, int *, NODE, double, double *);
int	op_deletion(int, int, void *, int *, NODE, double, double *);
int	op_reversal(int, int, void *, int *, NODE, NODE, double, double *);
void	random_permute(int *, int);
double *max_reduction(double *, double *, double *);
int	apply_action(double[], NODE *, int *, int, int, int[], int);
int	add_parent(int, NODE *);
int	remove_parent(int, int, NODE *);

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
void estimate_dag(PARAMS parms, int *G, int *iters)
{
	_CONFIG_ERROR status = E_SUCCESS;
	//struct timeval;
#ifdef PAR
	omp_set_num_threads(NTHREADS);
#endif

	// READ ONLY DATA ##########################################################
	// make local copies of parameter data
	const int p = parms.p;
	const int n = parms.n;
	const int m = parms.m;
	const int r = parms.r;
	const int max_parents = parms.max_parents;
	double *X = parms.X;
	NODE *Y = parms.Y;

	// create shuffled list of network node indices
	int candidates[p];
	for (int i = 0; i < p; ++i) {
		candidates[i] = i;
	}
	random_permute(candidates, p);

	// WRITE MULTIPLE ###########################################################
	int u = -1;	// child node index in global adjacency matrix G
	int u_M = -1;	// child node index in local adjacency matrix G_M
	int v = -1;	// parent node index in global adjacency matrix G
	int v_M = -1;	// parent node index in local adjacency matrix G_M

	int *G_M = Calloc(int, m * m); // local adjacency matrix
	// initialize local adjacency matrix G_M
	#pragma omp parallel for shared(G_M) private(u,v)
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			u = candidates[i];
			v = candidates[j];
			matrix(G_M, m, j, i) = matrix(G, p, v, u);
		}
	}

#ifdef PAR
	//double	*queue = Calloc(double, NTHREADS * 4);
	double	*queue[NTHREADS];
	void	*score_buffers[NTHREADS];
	int	*parent_buffers[NTHREADS];
	double	timing[NTHREADS];
	struct timeval btimers[NTHREADS];
	struct timeval atimers[NTHREADS];
	for (int i = 0; i < NTHREADS; ++i) {
		queue[i] = Calloc(double, 4);
		score_buffers[i] = BDE_init(X, X, p, n, r, max_parents);
		parent_buffers[i] = Calloc(int, max_parents);
		timing[i] = 0;
	}
#else
	//TODO could this be more elegant? (naming for default serial buffers)
	void	*score_buff = BDE_init(X, X, p, n, r, max_parents);
	int	*parent_buff = Calloc(int, max_parents);
#endif

	// ######## begin procedure #################################################
	/*
	 * repeatedly apply HC on the candidate set until no no_improvement_cnt
	 */
	int iteration = 0;
	int no_improvement_cnt = 0;
	//while (iteration < MAX_ITER) {
	//while ((no_improvement_cnt < IMP_LIMIT && (iteration < MAX_ITER))) {
	while ((no_improvement_cnt < (IMP_LIMIT*NTHREADS) && (iteration < (MAX_ITER/NTHREADS)))) {


		#pragma omp parallel \
		     firstprivate(u, u_M, v, v_M, status) \
		     shared(G_M, X, Y, candidates, queue, score_buffers, parent_buffers)
		{
		double wbegin = 0;
		double wend = 0;

		void    *score_buffer;
		int     *parent_buffer;
		double  score;

		// {delta score, action, candidate parent index, child index}
		double add_op[4] = { TOL, 1, -1, -1 };
		double del_op[4] = { TOL, 2, -1, -1 };
		double rev_op[4] = { TOL, 3, -1, -1 };

#ifdef PAR	/* INITIALIZE BUFFERS */
		const int tid = omp_get_thread_num();
		//printf(".");
		printf("t/all %d/%d, ", tid, omp_get_num_threads());
		score_buffer = score_buffers[tid];
		parent_buffer = parent_buffers[tid];
		u_M = (int)randinter(0, (m + 0.99));
		//u_M = (int)randinter(0, m);

		// chooses candidate at random from candidate set
		//TODO this is not evenly distributed
		//TODO error check rand within bounds
		//TODO fix the random generator
/*
		int rand = -1;
		while (!(rand >= 0 && rand <= m/NTHREADS)) {
			if (tid < (m % NTHREADS)) {
				rand = (int)randinter(0, (m/NTHREADS + 0.5));
				//printf("tid was %d, rand was %d MNTHR was %d\n", tid, rand, m/NTHREADS);
			} else {
				rand = (int)randinter(0, ((m/NTHREADS + 0.5) - 1));
				//printf("tid was %d, rand was %d MNTHR was %d\n", tid, rand, m/NTHREADS);
			}
		}
		//assert(rand >= 0 && rand <= m/NTHREADS);
		u_M = NTHREADS * rand + tid;
		//TODO remove this check
		if (!(rand >= 0 && rand < (m/NTHREADS+0.5)))
			syslog(LOG_INFO, "rand:%d, m/N:%d, tid:%d\n", rand, (int)(m/NTHREADS+0.5), tid);
		//assert(u_M >= 0 && u_M < m);
		//if (tid > 0 && u_M > 990) printf("TID WAS %d, u_M: %d\n", tid, u_M);
*/
#else		/* INITIALIZE SERIAL BUFFERS */
		score_buffer = score_buff;
		parent_buffer = parent_buff;
		u_M = randinter(0, (m + 0.99));
#endif

		u = candidates[u_M];
		assert(u >= 0 && u <= p);
		//if (!(u_M >= 0 && u_M <= p) || !(u >= 0 && u <= p))
		//	syslog(LOG_ERROR, "u is: %d, u_M is: %d; GETTING SCORE:", u, u_M);

		// initial score of candidate
		//printf("Y[%d] parents(cnt:%d) were: ", u, Y[u].num_parents);
		//for (int z = 0; z < max_parents; ++z)
		//    printf("%d, ", Y[u].parents[z]);
		//printf("\n");

		//TODO code this to a memory matrix to hold prev. calculated score
#ifdef PAR
		#ifdef DEBUG
		syslog(LOG_INFO, "%d :: u is: %d, u_M is: %d; GETTING SCORE:", tid, u, u_M);
		#endif
#endif
		//printf("buffer: %p, %p\n", score_buffer, score_buffers[tid]);
		score = get_score(score_buffer, u, Y[u].parents, Y[u].num_parents);
		//syslog(LOG_INFO, " %f\n", score);

		// reset all of the operation buffers
		add_op[0] = TOL;
		del_op[0] = TOL;
		rev_op[0] = TOL;
		add_op[2] = -1;
		del_op[2] = -1;
		rev_op[2] = -1;

//wbegin=omp_get_wtime();
//gettimeofday(&t1,0);
//gettimeofday(&btimers[tid],0);
wbegin=omp_get_wtime();
		for (v_M = 0; v_M < m; ++v_M) {
			v = candidates[v_M];
			assert(v >= 0 && v <= p);
			if (u == v) {
				continue;
			}

			/*
			int loop = 0;
			while (loop < 1000000) { loop++; }
			*/

			// if edge doesnt exist and parent set is not full
			if ((matrix(G_M, m, v_M, u_M) == 0) && \
			   (Y[u].num_parents < max_parents)) {
				status = op_addition(v_M, v, score_buffer, parent_buffer, Y[u], score, add_op);
				if (status != E_SUCCESS)
					syslog(LOG_PERROR, "error code %d in op_addition()", status);
			}

			// if edge exists and parent set is not empty
			//TODO redundant?? edge exists and still check parent set?
			if ((matrix(G_M, m, v_M, u_M) == 1) && \
			   (Y[u].num_parents > 0))
			{
				status = op_deletion(v_M, v, score_buffer, parent_buffer, Y[u], score, del_op);
				if (status != E_SUCCESS)
					syslog(LOG_PERROR, "error code %d in op_deletion()", status);
/*
				// if v has room for another parent and u isnt already a parent
				if ((Y[v].num_parents < max_parents) && \
				   (matrix(G_M, m, u_M, v_M) == 0))
				{
					status = op_reversal(u_M, v_M, score_buffer, parent_buffer, Y[u], Y[v], score, rev_op);
					if (status != E_SUCCESS)
						syslog(LOG_PERROR, "error code %d in op_reversal()", status);
				}
*/
			}

		}
wend=omp_get_wtime();
timing[tid] += (wend-wbegin);
//printf("writing time: %f to ID %d\n", (wend-wbegin), tid);
//wend=omp_get_wtime();
//gettimeofday(&t2,0);
//gettimeofday(&atimers,0);
//timing[tid] += t2.tv_sec+t2.tv_usec/1e6-(t1.tv_sec+t1.tv_usec/1e6);
//timing[tid] += atimers[tid].tv_sec+atimers[tid].tv_usec/1e6-(btimers[tid].tv_sec+atimers[tid].tv_usec/1e6);
//timing[tid] += (double)(wend-wbegin);

		// choose max of three operations and apply the corresponding action to the graph
		double *action = max_reduction(add_op, del_op, rev_op);

#ifndef PAR	/* SERIAL LOOP ACTION */
		/* apply the action that was chosen */
		if (action == NULL) {
			no_improvement_cnt++;

		} else if (action[2] != -1) {
			#ifdef DEBUG
			printf("\n========================================\n");
			printf("CURRENT SCORE (of %d, #parents %d): %f\n", u, Y[u].num_parents, score);
			if ((int)action[1] == 3) {
				int parent = candidates[(int)action[2]];
				score = get_score(score_buffer, parent, Y[parent].parents, Y[parent].num_parents);
				printf("CURRENT SCORE (of %d, #parents %d): %f\n", parent, Y[parent].num_parents, score);
			}
			#endif

			status = apply_action(action, Y, G_M, max_parents, u_M, candidates, m);
			if (status != E_SUCCESS) {
				printf("ABORTED\n");
				syslog(LOG_PERROR, "ERROR CODE %d IN apply_action(), action code: %d", status, (int)action[1]);
			}

			#ifdef DEBUG
			score = get_score(score_buffer, u, Y[u].parents, Y[u].num_parents);
			if ((int)action[1] == 3) {
				int parent = candidates[(int)action[2]];
				double score_other = get_score(score_buffer, parent, Y[parent].parents, Y[parent].num_parents);
				printf("SCORE AFTER: %f, %f (diff: %f)\n", score, score_other, action[0]);
			} else {
				printf("SCORE AFTER: %f (diff: %f)\n", score, action[0]);
			}
			#endif

			#ifdef DEBUG
			//printf("->%d", no_improvement_cnt);
			#endif
			no_improvement_cnt = 0;

		} else {
			// do nothing
		}
#endif

#ifdef PAR	/* PARALLEL LOOP ACTION */
		/* add the chosen action to the queue */
		if (action == NULL) {
			queue[tid][2] = -1;

		} else if (action[2] != -1) {
			action[3] = u_M;
			memcpy(queue[tid], action, sizeof(double)*4);

		} else {
			queue[tid][2] = -1;
		}

#endif
		} //END OMP PARALLEL
#ifdef PAR
		//printf("number of threads: %d\n", NTHREADS);
		//TODO check for redundant or conflicting thread action choices?
		//TODO sort by score first
		for(int i = 0; i < NTHREADS; ++i) {
			//if (queue[i] == NULL) {
			if (queue[i][2] == -1) {
				no_improvement_cnt++;
			} else if (queue[i][2] != -1) {
				//printf("idx: %d, :ITEMS %f, %f, %f, child: %f\n", h, queue[h][0], queue[h][1], queue[h][2], queue[h][3]);
				status = apply_action(queue[i], Y, G_M, max_parents, queue[i][3], candidates, m);
				if (status != E_SUCCESS) {
					syslog(LOG_PERROR, "ERROR CODE %d IN apply_action(), action code: %d, while: %d P:%dC:%d", status, (int)queue[i][1],iteration,(int)queue[i][2],(int)queue[i][3]);
					syslog(LOG_INFO, "idx: %d, :ITEMS %f, %f, child: %f", i, queue[i][0], queue[i][2], queue[i][3]);
				}
				#ifdef DEBUG
				printf("->%d", no_improvement_cnt);
				#endif
				no_improvement_cnt = 0;
			}
		}
#endif

//wend=omp_get_wtime();
//wtime_spent += (double)(wend-wbegin);


		/* END OF LOOP PREP */
		#ifdef DEBUG
		error_check(G_M, m, Y, candidates);
		#endif
		++iteration;
	} // while improvement

for (int i=0; i < NTHREADS; ++i) {
	printf("%f, ", timing[i]);
	wtime_spent += timing[i];
}
printf("\n");

printf("** WALL TIME ON WHILE LOOP: %f / %f(avg single iter) **\n", wtime_spent, wtime_spent/iteration);
//#ifdef DEBUG
printf("** FOR %d ITERATIONS **\n", iteration);
//#endif
	*iters = iteration;

#ifdef PAR
	for (int i = 0; i < NTHREADS; ++i) {
		BDE_finalize(score_buffers[i]);
		free(parent_buffers[i]);
		free(queue[i]);
	}
#else
	BDE_finalize(score_buff);
	free(parent_buff);
#endif

	// reflect changes on G_M back into the global adjacency matrix G
	#pragma omp parallel for shared(G_M) private(u,v)
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < m; ++j) {
			u = candidates[i];
			v = candidates[j];
			matrix(G, p, v, u) = matrix(G_M, m, j, i);
		}
	}

	free(G_M);
	//TODO return G_M and C
}

void error_check(int	*G_M,
		 int	m,
		 NODE	*Y,
		 int	candidates[])
{
	for (int i = 0; i < m; ++i) {
		//for(int jj = 0; jj < Y[candidates[ii]].num_parents; ++jj)
		for (int j = 0; j < m; ++j) {
			if (matrix(G_M, m, j, i) == 1) {
				NODE u = Y[candidates[i]];
				//NODE v = Y[candidates[jj]];

				int g_idx = -1;
				for(int k = 0; k < u.num_parents; ++k) {
					if(u.parents[k] == candidates[j]) {
						g_idx = k;
						break;
					}
				}

				if (g_idx == -1)
					syslog(LOG_PERROR, "AFTER :: PARENT %d(%d) OF %d(%d) DOESNT MATCH G_M\n",j,candidates[j],i,candidates[i]);
			}
		}
	}

	for (int i = 0; i < m; ++i) {
		NODE u = Y[candidates[i]];
		for (int j = 0; j < u.num_parents; ++j) {
			int idx = -1;
			for (int k = 0; k < m; ++k) {
				if (candidates[k] == u.parents[j]) {
					idx = k;
				}
			}

			if (idx != -1) {
				if (matrix(G_M, m, idx, i)!= 1) {
					syslog(LOG_PERROR, "AFTER :: G_M DOESNT MATCH PARENT (%d)%d OF (%d)%d\n",idx,candidates[idx],i,candidates[i]);
					/*
					printf("GM::::::::::::\n");
					for (int p = 0; p < m; ++p) {
						for (int l = 0; l < m; ++l)
							printf("%d, ",matrix(G_M,m,p,l));
						printf("\n");
					}
					*/
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
int op_addition(int	v_M,
		int	v,
		void	*buffer,
		int	*parent_buffer,
		NODE	u,
		double	current_score,
		double	max_diff_a[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	memcpy(parent_buffer, u.parents, sizeof(u.parents[0]) * (u.num_parents+1));
	parent_buffer[u.num_parents] = v;

	double new_score = get_score(buffer, u.index, parent_buffer, (u.num_parents+1));
	double diff = current_score - new_score;

	if (diff > max_diff_a[0]) {
		/* show the difference of new op vs current */
		//#ifdef DEBUG
		//printf("diff of %d(%f) > than current diff %f\n",v,diff,max_diff_a[0]);
		//#endif
		max_diff_a[0] = diff;
		max_diff_a[2] = v_M;
	}

	return status;
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
int op_deletion(int	v_M,
		int	v,
		void	*buffer,
		int	*parent_buffer,
		NODE	u,
		double	current_score,
		double	max_diff_d[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	memset(parent_buffer, -1, sizeof(u.parents[0]) * (u.num_parents-1));

	// create parent buffer with all of u's parents except for v
	int j = 0;
	for (int i = 0; i < u.num_parents; ++i) {
		if (u.parents[i] != v) {
			parent_buffer[j++] = u.parents[i];
		}
	}

	// if for some reason v was not found in u's parent list return error
	if (j > u.num_parents) {
		status = E_NOT_FOUND;
		return status;
	}
	/*
	printf("REMOVE [ - %d] testing parents size %d: ",k, u.num_parents-1);
	for (int z = 0; z < u.num_parents-1; ++z)
	    printf("%d, ", temp[z]);
	printf("\n");
	*/
	double new_score = get_score(buffer, u.index, parent_buffer, (u.num_parents-1));
	double diff = current_score - new_score;

	if (diff > max_diff_d[0]) {
		max_diff_d[0] = diff;
		max_diff_d[2] = v_M;
	}

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
int op_reversal(int	u_M,
		int	v_M,
		void	*buffer,
		int	*parent_buffer,
		NODE	u,
		NODE	v,
		double	current_score,
		double	max_diff_r[])
{
	_CONFIG_ERROR status = E_SUCCESS;

	memset(parent_buffer, -1, sizeof(u.parents[0]) * (u.num_parents-1));

	// create parent buffer with all of u's parents except for v
	int j = 0;
	for (int i = 0; i < u.num_parents; ++i) {
		if (u.parents[i] != v.index) {
			parent_buffer[j++] = u.parents[i];
		}
	}

	// if for some reason v was not found in u's parent list return error
	if (j > u.num_parents) {
		status = E_NOT_FOUND;
		return status;
	}

	double new_score_d1 = get_score(buffer, u.index, parent_buffer, (u.num_parents-1));
	double d1 = current_score - new_score_d1;

	memcpy(parent_buffer, v.parents, sizeof(v.parents[0]) * (v.num_parents+1));
	parent_buffer[v.num_parents] = u.index;

	/*
	#ifdef DEBUG
	printf("REV [ + %d ] testing parents size %d: ", v.index, v.num_parents+1);
	for (int z = 0; z < v.num_parents+1; ++z)
	    printf("*%d*::%d,  ", parent_buffer[z], v.parents[z]);
	printf("\n");
	#endif
	*/
	double current_rev_score = get_score(buffer, v.index, v.parents, v.num_parents);
	double new_score_d2 = get_score(buffer, v.index, parent_buffer, (v.num_parents+1));
	double d2 = current_rev_score - new_score_d2;

	double diff = d1 + d2;

	if (diff > max_diff_r[0]) {
		/* show diff of new op vs current */
		/*
		#ifdef DEBUG
		printf("current: %f, new_d2 %f === ", current_score, new_score_d2);
		printf("diff of d1 %d(%f) + d2 %d(%f) \n", u.index, d1, v.index, d2);
		#endif
		*/
		max_diff_r[0] = diff;
		max_diff_r[2] = v_M;
	}

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
void random_permute(int	*m,
		    int	n)
{
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
double * max_reduction(double *a,
		       double *b,
		       double *c)
{
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
int apply_action(double	action[3],
		 NODE	*Y,
		 int	*G_M,
		 int	max_parents,
		 int	child,
		 int	candidates[],
		 int	m)
{
	_CONFIG_ERROR status = E_SUCCESS;

	int parent = action[2];
	int action_code = (int)action[1];
	int u = candidates[child];
	int v = candidates[parent];

	switch (action_code) {
	case 1: // edge addition
		status = add_parent(v, &Y[u]);
		if (status != E_SUCCESS) {
			break;
		}
		matrix(G_M, m, parent, child) = 1;

		#ifdef DEBUG
		printf("ADDING PARENT %d to CHILD %d\n",v,u);
		printf("CHILD parents after: ");
		for (int z = 0; z < max_parents; ++z) {
			printf("%d, ", Y[u].parents[z]);
		}
		printf("\n");
		#endif

		break;

	case 2: // edge removal
		status = remove_parent(v, max_parents, &Y[u]);
		if (status != E_SUCCESS) {
			break;
		}
		matrix(G_M, m, parent, child) = 0;

		#ifdef DEBUG
		printf("REMOVING PARENT %d FROM CHILD %d ", v, u);
		printf("CHILD parents after(%d): ", Y[u].num_parents);
		for (int z = 0; z < max_parents; ++z) {
			printf("%d, ", Y[u].parents[z]);
		}
		printf("\n");
		#endif

		break;

	case 3: // edge reversal
		//TODO what to do for reversal if target parent list is already full??
		if (Y[v].num_parents < max_parents) {
			// remove parent from u parent list
			status = remove_parent(v, max_parents, &Y[u]);
			if (status != E_SUCCESS) {
				syslog(LOG_PERROR, "error in reverse reverse() remove parent");
				break;
			}

			// add u parent parent list
			status = add_parent(u, &Y[v]);
			if (status != E_SUCCESS) {
				syslog(LOG_PERROR, "TODO revert remove:: error in reverse() add parent");
				break;
			}

			matrix(G_M, m, parent, child) = 0;
			matrix(G_M, m, child, parent) = 1;

			#ifdef DEBUG
			printf("CHILD parents after(%d): ", Y[u].num_parents);
			for (int z = 0; z < max_parents; ++z) {
				printf("%d, ", Y[u].parents[z]);
			}
			printf("\n");
			printf("TARGET parents after(%d): ", Y[v].num_parents);
			for (int z = 0; z < max_parents; ++z) {
				printf("%d, ", Y[v].parents[z]);
			}
			printf("\n");
			#endif
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
int add_parent(int	v,
	       NODE	*Y)
{
	_CONFIG_ERROR status = E_SUCCESS;

	// dont add a duplicate parent
	for (int i = 0; i < Y->num_parents; ++i) {
		if (Y->parents[i] == v) {
			status = E_DUPLICATE;
			return status;
		}
	}

	int tail = Y->num_parents;
	Y->parents[tail] = v;
	Y->num_parents++;

	return status;
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
int remove_parent(int	v,
		  int	max_parents,
		  NODE	*Y)
{
	_CONFIG_ERROR status = E_NOT_FOUND;

	for (int i = 0; i < max_parents; ++i) {
		// find target parent in the list
		if (Y->parents[i] == v) {
			// shift all parents above target index down (overwrite target)
			for (int k = i; k < (max_parents - 1); ++k) {
				Y->parents[k] = Y->parents[k + 1];
			}

			// set the last element to empty and decrease count
			Y->parents[(max_parents - 1)] = -1;
			Y->num_parents--;
			status = E_SUCCESS;
			break;
		}
	}

	return status;
}

