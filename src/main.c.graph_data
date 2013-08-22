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

void free_nodes(NODE *x)
{
	free((*x).data);
	free((*x).name);
	free(x);
	x = NULL;
}

int main(int argc, char *argv[])
{
	int status = 1;
	char input_file_name[1024] =
			"/home/mbarga/Dropbox/school/code/hc/data/toy.txt";
	time_t current_time;

	int i, j, k, u, v;
	int x_size = 0;
	int m_size = 5; // TODO pass this in as argument, etc.?

	NODE *x = NULL; // set of all input nodes
	int *m = NULL; // array of permuted indices for the nodes X
	int *g = NULL; // adjacency matrix of the network
	int *c = NULL; // node selection count matrix

	double max_score = 0;

	// initialize logging
	openlog(PROG_TAG, 0, LOG_USER);
//	syslog(LOG_INFO, "%s LOG INITIALIZED\n", ctime(&current_time));

	// read in data samples from file
	status = read_problem(input_file_name, &x, &x_size);

	if (status != 0)
	{
		syslog(LOG_INFO, "%s", "main(): failed to read input file\n");
		return 1;
	}

	if (x == NULL)
	{
		syslog(LOG_INFO, "%s", "main(): x was returned NULL\n");
		return 1;
	}

/*	
	 printf("--- Printing %d nodes ---\n", x_size);
	 for(i=0; i<x_size; i++) {
		 printf("%s :: ",x[i].name);
		 for (j=0; j<2; ++j) {
			 printf("%f, ",x[i].data[j]);
		 }
		 printf("\n");
	 }
 */

	g = Malloc(int, x_size*x_size);
	c = Malloc(int, x_size*x_size);

	// initialize c and g matrices <- 0
	for (i = 0; i < x_size; ++i)
	{
		for (j = 0; j < x_size; ++j)
		{
			g[i * x_size + j] = 0;
			c[i * x_size + j] = 0;
		}
	}

	/*********************** begin hill climbing algorithm ****************************/
	// repeat for a set number of iterations
	for (i = 0; i < MAX_ITER; ++i)
	{
		if (m_size > x_size)
		{
			syslog(LOG_INFO, "main(): |M|=%d was larger than |X|=%d", m_size,
					x_size);
			free(g);
			free(c);
			free_nodes(x);
			return 1;
		}

		m = randperm(m_size, x_size);
		max_score = 0;
/*
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
*/
		//TODO should the order visited not be linear even though the subset has already been randomized?
		//TODO how to encode the edges before sending to G?
		// check for each operation (add, delete, reverse)
		for (j = 0; j < m_size; ++j)
		{
			for (k = 0; k < m_size; ++k)
			{
				if (j == k) continue; // node cannot be its own parent
				// m = randperm(m_size, m_size); // re-randomize the indices
				u = m[j];
				v = m[k];
				++c[u * x_size + v];
				++c[v * x_size + u];

				if (score(x[u], x[v]) > max_score) 
				{ // adding edge increases score
				}
				else if (score() > max_score)
				{ // deleting edge increases score
				}
				else if (score(x[v], x[u]) > max_score)
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
	free_nodes(x);
	syslog(LOG_INFO, "%s", "exiting cleanly\n");
	closelog();
	return 0;
}

int read_problem(const char *filename, NODE **x, int *n)
{
	FILE *fp = fopen(filename, "r");
	char *p;
	int isspace(int c);
	char *endptr;
	char *val, *label;

	double *data_space = NULL;
	char *name_space = NULL;
	NODE *result = NULL;

	int numFeatures = 0;
	int numNodes = 0;
	int stop = 0;
	int i, j;

	if (fp == NULL)
	{
		syslog(LOG_INFO, "%s", "read_problem(): couldnt open input file\n");
		return 1;
	}

	max_line_len = 1024;
	line = Malloc(char, max_line_len);

	while (read_line(fp) != NULL)
	{
		p = strtok(line, " \t\n"); // toss label
		while (1)
		{
			p = strtok(NULL, "\t"); // read following tokens
			if (p == NULL || *p == '\n')
			{ // check '\n' as ' ' may be after the last feature
				if (stop == 0)
				{
					stop = 1;
					++numFeatures;
				}
				break;
			}
			if (stop == 0)
			{
				++numFeatures;
			}
		}
		++numNodes;
	}

	//TODO count features properly?
	numFeatures--;
	rewind(fp);
	syslog(LOG_INFO, "Found %d nodes, each with %d features\n", numNodes,
			numFeatures);

	result = Malloc(NODE, numNodes);
	data_space = Malloc(double, numNodes*numFeatures);
	name_space = Malloc(char, numNodes*max_name_len);

	for (i = 0; i < numNodes; ++i)
	{
		read_line(fp);
		label = strtok(line, " \t\n");
		if (label == NULL)
		{
			syslog(LOG_INFO, "%s", "FAILED::read_problem : empty line\n");
			return 1;
		}

		// copy node name label to the node itself
		strcpy(&name_space[i * max_name_len], label);

		for (j = 0; j < numFeatures; ++j)
		{
			val = strtok(NULL, "\t");
			//TODO val = strtok(NULL, &endptr);

			if (val == NULL)
				break;

			data_space[i * numFeatures + j] = strtod(val, &endptr);

			if (endptr == val || (*endptr != '\0' && !isspace(*endptr)))
			{
				syslog(LOG_INFO, "%s",
						"read_problem(): endpointer/val mismatch\n");
			}
		}
		result[i].name = &name_space[i * max_name_len];
		result[i].data = &data_space[i * numFeatures];
	}

	*x = result;
	result = NULL;
	*n = numNodes;

	free(line);
	fclose(fp);
	return 0;
}

static char* read_line(FILE *input)
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
