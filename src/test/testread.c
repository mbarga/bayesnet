#include "readfile.h"

#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char *argv[])
{
    char input_file_name[1024] =
                "/home/mbarga/Dropbox/a_code/bayesnet/data/toy2.txt";
    double *x = NULL;
    int n;
    int m;

    int status = read_problem(input_file_name, &x, &n, &m);

    if (status == 1)
        return 1;

    printf("n was: %d, m was: %d\n\n",n,m);

    for (int i = 0; i<n; ++i)
    {
        for (int j = 0; j<m; ++j)
            printf("%f, ", x[i * m + j]);

        printf("\n");
    }

    if (x) free(x);
    return 0;
}
