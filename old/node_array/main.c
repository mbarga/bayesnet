#include "main.h"
#include "library.h"
#include "readfile.h"

#define PROG_TAG "MAIN"

void print_matrix(int *c, int size)
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

void print_nodes(GHashTable *hash)
{
	GList * keys = NULL;
	GSList * edge = NULL;
	NODE * test = NULL;

	keys = g_hash_table_get_keys(hash);
	while (keys != NULL)
	{
		test = g_hash_table_lookup(hash, keys->data);
		if (test != NULL)
		{	
			printf("node :: %s with index %d, edges: ", test->name, test->index);
			edge = test->edges;
			while (edge != NULL) 
			{
				test = edge->data;
				printf("%s, ",test->name);
				edge = edge->next;
			}
			printf("\n");
		}
		keys = keys->next;
	}

	g_list_free(keys);
	keys = NULL;
	edge = NULL;
	test = NULL;
}

void free_nodes(GHashTable *hash)
{
	GList *keys = g_hash_table_get_keys(hash);
	NODE *test = NULL;

	while (keys != NULL)
	{
		test = g_hash_table_lookup(hash, keys->data);
		if (test != NULL)
		{
			g_slist_free(test->edges);
			g_strfreev(&test->name);
			free(test);
		}
		keys = keys->next;
	}

	test = NULL;
	g_hash_table_destroy(hash);
	hash = NULL;
	g_list_free(keys);
	keys = NULL;
}

int main(int argc, char *argv[])
{

	int status = 1;
	char input_file_name[1024] =
			"/home/mbarga/Dropbox/git/bayesnet/data/toy2.txt";

	int i, j, u, v, n;
	int x_size = 0;
	int m_size = 3;

	GHashTable* hash = NULL;
	int *m = NULL; // array of permuted indices for the nodes X
	int *g = NULL; // adjacency matrix of the network
	int *c = NULL; // node selection count matrix
	NODE *x = NULL; // array of node structs

	double *l_scores = NULL;

	// initialize logging
	openlog(PROG_TAG, 0, LOG_USER);

	/*
	 * read in data samples from file and check for possible errors
	 */ 
	status = read_problem(input_file_name, &hash, &x);

	if (status != 0)
	{
		syslog(LOG_INFO, "%s", "main(): failed to read input file\n");
		return 1;
	}

	if (hash == NULL)
	{
		syslog(LOG_INFO, "%s", "main(): x was returned NULL\n");
		return 1;
	}

	x_size = g_hash_table_size(hash);
	
	if (m_size > x_size)
	{
		syslog(LOG_INFO, "main(): |M|=%d was larger than |X|=%d", m_size,
				x_size);
		free_nodes(hash);
		free(x);
		return 1;
	}

	/*
	 * allocate space for and initialize c and g matrices <- 0
	 */
	g = Malloc(int, x_size*x_size);
	c = Malloc(int, x_size*x_size);

	for (i = 0; i < x_size; ++i)
	{
		for (j = 0; j < x_size; ++j)
		{
			g[i * x_size + j] = 0;
			c[i * x_size + j] = 0;
		}
	}

	/*
	 * calculate one to one local scores
	 *  (i,j)th element is the local score() of graph gene_i -> gene_j
	 */ 
	one_to_one(&hash, &l_scores);

	/*
	 * repeated random sampling on subsets of X
	 */ 
	m = Malloc(int, x_size);
	for (i=0; i < x_size; ++i) m[i] = i;

	n=0;
	while(n < MAX_SUBSETS)
	{	
		randperm(m, m_size, x_size);

		/* 
		 * c_(u,v) <- c_(u,v) + 1 AND c_(v,u) <- c_(u,v) for all (X_u, X_v) in M
		 */
		for (i = 0; i < m_size; ++i)
		{
			for (j = i+1; j < m_size; ++j)
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


	// cleanup before exiting the algorithm
	//print_matrix(c, x_size);
	//print_matrix(g,x_size);
	free_nodes(hash);
	free(m);
	free(g);
	free(c);
	free(x);
	free(l_scores);

	syslog(LOG_INFO, "%s", "exiting cleanly\n");
	closelog();
	return 0;
}

void estimate_dag(GHashTable **hash, int *m, int m_size)
{
	int i=0, j=0, k=0;
	double max_score = 0;
	char score_changed = 1;

	GList *keys = g_hash_table_get_keys(*hash);
	NODE *u, *v;

	while (score_changed == 1 && i < MAX_ITER)
	{
		// jumble the indices in m[] since we want to visit nodes in random order
		randperm(m, m_size, m_size);
		
		// for all members in the randomly chosen subset M
		for (j = 0; j < m_size; ++j)
		{
			u = g_hash_table_lookup(*hash, g_list_nth_data(keys, m[j]));
			
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
			
			for (k = 0; k < m_size; ++k)
			{
				if (j == k) continue; // node cannot be its own parent
		
				v = g_hash_table_lookup(*hash, g_list_nth_data(keys, m[k]));

				// test each operation (add, delete, reverse)
			//	if (no edge)
			//	{
					// calculate added edge score
					// calculate reversed edge score
			//	}
			//	else if (edge)
			//	{
					// calculate deleted edge score
			//	}

				
			//	{
			//		score_changed = 0;
			//	}
			}
		}

		++i;
	}

	g_list_free(keys);
	keys = NULL;
	u = NULL;
	v = NULL;
}
