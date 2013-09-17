#include "main.h"
#include "library.h"
#include "bdescore.h"

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <math.h>
#include <glib.h>

#define PROG_TAG "LIBRARY"

/* Function
 * -------------------
 *
 * TODO: return only array of size k OR optimize to shuffle only for k elements
 * Returns array of size N of random numbers in the range [1,N]
 *
 */
void randperm(int *m, int k, int n) {
	int i, j, t;

	/*
	 for (i = 0; i < n; i++)
	 {
	 m[i] = i;
	 }
	 */
	if (m == NULL) {
		errlog("randperm(): m was found NULL\n");
		return;
	}

	for (i = 0; i < n; i++) {
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
void scramble(int **indices, int size) {
	int i, j, t;

	for (i = 0; i < size; ++i) {
		j = rand() % (size - i) + i;
		t = *indices[j];
		*indices[j] = *indices[i];
		*indices[i] = t;
	}
}

double * max(double *a, double *b, double *c)
{
	double temp = 0;
	char action = 0;

	if (a[0] > temp) {temp = a[0]; action = 1;}
	if (b[0] > temp) {temp = b[0]; action = 2;}
	if (c[0] > temp) {temp = c[0]; action = 3;}

	if (action == 1) return a;
	else if (action == 2) return b;
	else if (action == 3) return c;
	else return NULL;
}

/* Function
 * -------------------
 *
 * TODO use array of nodes instead
 * calculate all 1-to-1 scores
 *
 */
void one_to_one(NODE **x, int x_size, double **local_scores) {
	int u, v;
	double *scores = scores = Malloc(double, x_size*x_size);

	// loop over each possible child-pair edge
	for (u = 0; u < x_size; ++u) {
		for (v = 0; v < x_size; ++v) {
			scores[u * x_size + v] = calc_bdescore();
		}
	}

	*local_scores = scores;
	scores = NULL;
}

/* Function
 * -------------------
 *
 *
 */
int generate_local_probability() {
	//TODO the resolution of the time seed is not high enough
	int seed = time(NULL);
	srand(seed);
	// {-1, 0, 1}
	return (int) (((rand() % 3) - 1));
}

/* Function
 * -------------------
 *
 *
 */
void free_nodes(NODE *x, int x_size) {
	NODE *nodes = (NODE *) x;
	int i;

	if (nodes != NULL) {
		for (i = 0; i < x_size; ++i) {
			//free(x[i].parents);
			g_slist_free(nodes[i].parents);
			free(x[i].name);
		}
	}

	free(nodes);
	x = NULL;
	nodes = NULL;
}

/* Function
 * -------------------
 * TODO remove?
 * The following are only useful for debugging purposes
 *
 */
void print_matrix(int *c, int size) {
	int i, j;
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			printf("%d, ", c[i * size + j]);
		}
		printf("\n");
	}
}

void print_nodes(NODE *x, int x_size) {
	int i;
	GSList * edge = NULL;
	NODE* test;

	for (i = 0; i < x_size; ++i) {
		printf("node :: %s with index %d, edges: ", x[i].name, x[i].index);
		edge = x[i].parents;
		while (edge != NULL) {
			test = edge->data;
			printf("%s, ", test->name);
			edge = edge->next;
		}
		printf("\n");
	}

	g_slist_free(edge);
	edge = NULL;
}

void errlog(char *s)
{
	fprintf(stderr, "%s", s);
	//syslog(LOG_INFO, "%s", "randperm(): m was found NULL\n");
}
