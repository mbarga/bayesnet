#include "probability.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int
tri2dec(int[], int);

/**
 * store contingency tables for all parent set sizes
 * size of table given |p| is: (# categories) ^ |P|
 */

// (distinct tables, states)
// Malloc(double, 3^num_parents)
static const double p_1[5][1][3] =
{   {{0.1, 0.2, 0.7}},
    {{0.2, 0.1, 0.7}},
    {{0.8, 0.1, 0.1}},
    {{0.3, 0.5, 0.2}},
    {{0.1, 0.2, 0.7}}};

//static const double p_2[5][9]   =
//static const double p_3[5][27]  =

/* store conditional probability tables for each size of parent set */
int
get_expression(int *parents, int size)
{
    int seed = time(NULL);
    srand(seed);

    int i = rand() % 5;
    int j = tri2dec(parents, size);

    printf("i,j: %d, %d\n", i,j);

    const double *prob_vector = p_1[i][j];
    printf("sample from prob vector: %f\n", prob_vector[0]);

    // draw an expression level from p(x | parents)
    int e = 1;
    return e;
}

int
tri2dec(int x[], int n)
{
    int sum = 0;
    for (int i=0; i<n; ++i)
    {
        int exp = 3;
        if (i == 0)
            exp = 1;
        else
            for (int j=1; j<i; ++j)
                exp *= 3;

        sum += exp * x[i];
    }

    return sum;
}
