#define _GNU_SOURCE // prevents implicit getline() error on compilation

#include "readfile.h"
#include "util.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h> // size_t type definition
#include <unistd.h> // ssize_t type definition
#include <assert.h>

/*
##############################################################################*/
/**
 * @brief 
 *
 * @param filename
 * @param x
 * @param p
 * @param n
 *
 * @return 
 */
/*
##############################################################################*/
int read_problem(const char *filename, double **x, int *p, int *n)
{
	FILE *fp = fopen(filename, "r");
	int isspace(int c);
	char *line = NULL;
	char *pch = NULL;
	size_t len = 0;
	ssize_t read;

	_CONFIG_ERROR status = E_SUCCESS;

	assert(fp != NULL);

	/* first peek at number of nodes and samples */
	int local_p = 0;
	int local_n = 0;
	while ((read = getline(&line, &len, fp)) != -1)
	{
		pch = strtok (line, " \t\n");
		// this throws away first element in the row (node label)
		pch = strtok (NULL, " \t\n");

		local_n = 0;
		// strtok returns NULL when no tokens left to retrieve
		while (pch != NULL)
		{
			pch = strtok (NULL, " \t\n");
			++local_n;
		}

		++local_p;
	}

	rewind(fp);

	double *dataset = Calloc(double, local_p*local_n);

	/* now write values in file to matrix */
	//for (int i = 0; (read = getline(&line, &len, fp)) != -1; ++i)
	int i = 0;
	while ((read = getline(&line, &len, fp)) != -1)
	{
		// after this line, the node label is in the buffer (e.g. "g1")
		pch = strtok (line, " \t\n");

		// strtok returns NULL when no tokens left to retrieve
		//for (int j = 0; pch != NULL; ++j)
		int j = 0;
		while (pch != NULL)
		{
			pch = strtok (NULL, " \t\n");
			//if ((!isspace(*pch)) && (pch != '\0'))
			if (pch != '\0')
				//TODO revert this
				dataset[i * local_n + j] = (double)atof(pch);
				//dataset[j * local_p + i] = (double)atof(pch);
			++j;
		}

		++i;
	}

	if (line) free(line);

	*p = local_p;
	*n = local_n;
	*x = dataset;

	fclose(fp);

	if (dataset == NULL)
		status = E_NULL_POINTER;

	return status;
}
