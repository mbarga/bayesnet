#include "main.h"
#include "library.h"

#define PROG_TAG "library"

//TODO: setup to only return array of size k OR optimize to shuffle only for k elements
int * randperm(int k, int n)
{
	int i, j, t;
	int *a = (int *) malloc(n * sizeof(int));

	for (i = 0; i < n; i++)
	{
		a[i] = i;
	}

	for (i = 0; i < n; i++)
	{
		j = rand() % (n - i) + i;
		t = a[j];
		a[j] = a[i];
		a[i] = t;
	}

	return a;
}

// calculate all 1-to-1 scores
void one_to_one(GHashTable *hash, double **local_scores)
{	
	int u_index, v_index;
	NODE *u = NULL;
	NODE *v = NULL;

	GList *u_iter = g_hash_table_get_keys(hash);
	GList *v_iter = g_hash_table_get_keys(hash);

	int x_size = g_hash_table_size(hash);
	double *l_scores = NULL;

	l_scores = Malloc(double, x_size*x_size);

	// loop over each possible child-pair edge
	while (u_iter != NULL)
	{
		u = g_hash_table_lookup(hash, u_iter->data);

		if (u != NULL)
		{
			u_index = u->index;

			//small change

			while (v_iter != NULL)
			{
				v = g_hash_table_lookup(hash, v_iter->data);
				
				if (v != NULL)
				{
					v_index = v->index;
					
					local_score[u_index * x_size + v_index] = score(u, v);
				}

				v_iter = v_iter->next;
			}
			// restart parent iterator from the beginning of the list
			v_iter = g_hash_table_get_keys(hash);
		}

		u_iter = u_iter->next;
	}

	g_list_free(u_iter);
	g_list_free(v_iter);
	u 		= NULL;
	v 		= NULL;
	u_iter 	= NULL;
	v_iter 	= NULL;

	*local_scores = l_scores;
}

double score(NODE *u, NODE *v)
{
	double score = 0;

	return score;
}


