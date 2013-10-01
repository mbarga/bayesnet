#include "readfile.h"
#include "library.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

#define PROG_TAG "READFILE"

/* Function
 * -------------------
 *
 *
 */
int
read_problem(const char *filename, double **x, int *n, int *m)
{
    FILE *fp = fopen(filename, "r");
    int isspace(int c);
    char *line = NULL;
    char *pch = NULL;
    size_t len = 0;
    ssize_t read;

    if (fp == NULL)
    {
        errlog("read_problem(): couldnt open input file\n");
        return 1;
    }

    /* first peek at number of nodes and samples */
    int local_n = 0;
    int local_m = 0;
    while ((read = getline(&line, &len, fp)) != -1)
    {
        pch = strtok (line, " \t\n");

        // strtok returns NULL when no tokens left to retrieve
        while (pch != NULL)
        {
          ++local_m;
          pch = strtok (NULL, " \t\n");
        }

        ++local_n;
    }

    local_m = local_m / local_n;

    rewind(fp);

    double *dataset = Malloc(double, local_n*local_m);

    /* now write values in file to matrix */
    for (int i = 0; (read = getline(&line, &len, fp)) != -1; ++i)
    {
        pch = strtok (line, " \t\n");

        // strtok returns NULL when no tokens left to retrieve
        for (int j = 0; pch != NULL; ++j)
        {
            //if ((!isspace(*pch)) && (pch != '\0'))
            if (pch != '\0')
                dataset[i * local_m + j] = (double)atof(pch);

            pch = strtok (NULL, " \t\n");
        }
    }

    if (line) free(line);

    *n = local_n;
    *m = local_m;
    *x = dataset;

    fclose(fp);
    return 0;
}
