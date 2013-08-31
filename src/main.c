#include "main.h"
#include "library.h"

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
			free(test);
		}
		keys = keys->next;
	}

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
	//time_t current_time;

	int i, j, k, u, v;
	int x_size = 0;
	int m_size = 5; // TODO pass this in as argument, etc.?

	//NODE *x = NULL; // set of all input nodes
	GHashTable* hash = NULL;
	int *m = NULL; // array of permuted indices for the nodes X
	int *g = NULL; // adjacency matrix of the network
	int *c = NULL; // node selection count matrix
	double *l_scores = NULL;

	double max_score = 0;

	// initialize logging
	openlog(PROG_TAG, 0, LOG_USER);
//	syslog(LOG_INFO, "%s LOG INITIALIZED\n", ctime(&current_time));

	// read in data samples from file
	status = read_problem(input_file_name, &hash);

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
	//printf("There are %d keys in the hash\n", x_size);
	//print_nodes(hash);

	g = Malloc(int, x_size*x_size);
	c = Malloc(int, x_size*x_size);

	// initialize c and g matrices <- 0
	// initialize all 1-to-1 scores
	for (i = 0; i < x_size; ++i)
	{
		for (j = 0; j < x_size; ++j)
		{
			g[i * x_size + j] = 0;
			c[i * x_size + j] = 0;
		}
	}

	one_to_one(hash, );
	//
	//
	//
	//BOOKMARK
	//
	//
	//

	// repeat for a set number of iterations
	for (i = 0; i < MAX_ITER; ++i)
	{
		if (m_size > x_size)
		{
			syslog(LOG_INFO, "main(): |M|=%d was larger than |X|=%d", m_size,
					x_size);
			free(g);
			free(c);
			free_nodes(hash);
			return 1;
		}

		m = randperm(m_size, x_size);
		max_score = 0;

		for (j = 0; j < m_size; ++j)
		{
			for (k = j+1; k < m_size; ++k)
			{
				u = m[j];
				v = m[k];
				++c[u * x_size + v];
				++c[v * x_size + u];
			}
		}

		//TODO should the order visited not be linear even though the subset has already been randomized?
		//TODO how to encode the edges before sending to G?
		// check for each operation (add, delete, reverse)
		for (j = 0; j < m_size; ++j)
		{
			for (k = 0; k < m_size; ++k)
			{
				if (j == k) continue; // node cannot be its own parent
				// m = randperm(m_size, m_size); // re-randomize the indices
			//	u = m[j];
			//	v = m[k];
			//	++c[u * x_size + v];
		//		++c[v * x_size + u];

			//	if (score(x[u], x[v]) > max_score) 
				{ // adding edge increases score
				}
			//	else if (score() > max_score)
				{ // deleting edge increases score
				}
			//	else if (score(x[v], x[u]) > max_score)
				{ // reversing edge increases score
				}
			}
		}
		//TODO ensure no cycles?
		//TODO for each selected directed edge in the G_M DAG, add the edge to the adjacency matrix g

	}

	// cleanup before exiting the algorithm
	//print_matrix(c,x_size);
	//print_matrix(g,x_size);
	free(m);
	free(g);
	free(c);
	free_nodes(hash);
	
	syslog(LOG_INFO, "%s", "exiting cleanly\n");
	closelog();
	return 0;
}

int read_problem(const char *filename, GHashTable **x)
{
	FILE *fp = fopen(filename, "r");
	//int isspace(int c);
	//char *endptr;
	char *child_name, *parent_name;

	GHashTable* hash = g_hash_table_new(g_str_hash, g_str_equal);
	NODE *parent = NULL;
	NODE *child = NULL;

	int node_count = 0;

	if (fp == NULL)
	{
		syslog(LOG_INFO, "%s", "read_problem(): couldnt open input file\n");
		return 1;
	}

	max_line_len = 1024;
	line = Malloc(char, max_line_len);

	while (read_line(fp) != NULL)
	{
		child_name = strtok(line, " \t\n");
		parent_name	= strtok(NULL, "\t\n"); //TODO clean this up??
	//	child_name	= strtok(NULL, &endptr);
		if (parent_name == NULL || child_name == NULL)
		{
			syslog(LOG_INFO, "%s", "FAILED::read_problem : empty line\n");
			return 1;
		}
	
		child = g_hash_table_lookup(hash, child_name);
		if (child == NULL)
		{
			child 				= create_node();
			child->index	= node_count++;
			child->name 	= g_strdup(child_name);
 			g_hash_table_insert(hash, g_strdup(child_name), child); //TODO use node->name
		}
		
		parent = g_hash_table_lookup(hash, parent_name);
		if (parent == NULL)
		{
			parent 				= create_node();
			parent->index	= node_count++;
			parent->name 	= g_strdup(parent_name);
 			g_hash_table_insert(hash, g_strdup(parent_name), parent);
		}
	
		parent->edges = g_slist_append(parent->edges, child);
	}

	*x = hash;

	free(line);
	fclose(fp);
	return 0;
}

static char * read_line(FILE *input)
{
	int len;

	if (fgets(line, max_line_len, input) == NULL)
	{
		return NULL;
	}

	while (strrchr(line, '\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line, max_line_len);
		len = (int) strlen(line);
		if (fgets(line + len, max_line_len - len, input) == NULL)
			break;
	}

	return line;
}

NODE * create_node()
{
	NODE *p;
	//p = Malloc(NODE, 1);
	p 				= g_new(NODE, 1);
	p->name 	= NULL; 
	p->index 	= 0;
	p->edges	= NULL;
	return p;
}
