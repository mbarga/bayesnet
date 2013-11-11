#include "main.h"
#include "library.h"

#define PROG_TAG "library"

/*
 * Returns array of size N of random numbers in the range [1,N]
 */ 
void randperm(int *m, int k, int n)
{
	int i, j, t;
//int *a = (int *) malloc(n * sizeof(int));
//int *a = Malloc(int, n);

/*
	for (i = 0; i < n; i++)
	{
		m[i] = i;
	}
*/
	if (m == NULL)
	{
		syslog(LOG_INFO, "%s", "randperm(): m was found NULL\n");
		return;
	}

	for (i = 0; i < n; i++)
	{
		j = rand() % (n - i) + i;
		t = m[j];
		m[j] = m[i];
		m[i] = t;
	}

	//return a;
}

void scramble(int **indices, int size)
{
	int i, j, t;

	for (i = 0; i < size; ++i)
	{
		j = rand() % (size -i) + i;
		t = *indices[j];
		*indices[j] = *indices[i];
		*indices[i] = t;
	}
}

// calculate all 1-to-1 scores
void one_to_one(GHashTable **hash, double **local_scores)
{	
	int u_index, v_index;
	NODE *u = NULL;
	NODE *v = NULL;

	GList *u_iter = g_hash_table_get_keys(*hash);
	GList *v_iter = g_hash_table_get_keys(*hash);

	int x_size = g_hash_table_size(*hash);
	double *l_scores = NULL;

	l_scores = Malloc(double, x_size*x_size);

	// loop over each possible child-pair edge
	while (u_iter != NULL)
	{
		u = g_hash_table_lookup(*hash, u_iter->data);

		if (u != NULL)
		{
			u_index = u->index;

			//small change

			while (v_iter != NULL)
			{
				v = g_hash_table_lookup(*hash, v_iter->data);
				
				if (v != NULL)
				{
					v_index = v->index;
				
					if (u_index != v_index)
					{
						l_scores[u_index * x_size + v_index] = score(u, v);
					}
				}

				v_iter = v_iter->next;
			}
			// restart parent iterator from the beginning of the list
			v_iter = g_hash_table_get_keys(*hash);
		}

		u_iter = u_iter->next;
	}

	g_list_free(u_iter);
	g_list_free(v_iter);
	u_iter 	= NULL;
	v_iter 	= NULL;
	u 		= NULL;
	v 		= NULL;

	*local_scores = l_scores;
	l_scores = NULL;
}

double score(NODE *u, NODE *v)
{
	double score = rand() % 100;

	return score;
}

int generate_random_p_x ()
{
	return (int)(rand()%2);
}
