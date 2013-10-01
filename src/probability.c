#include "probability.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int
tri2dec(int[], int);
int
find_ceil(int arr[], int r, int l, int h);
int
my_rand(int arr[], int freq[], int n);

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

//TODO generate larger tables randomly and automatically??

/* store conditional probability tables for each size of parent set */
int
get_expression(int *parents, int size)
{
    srand(time(NULL));

    int i = rand() % 5;
    int j;

    // if no parents
    if (parents == NULL)
        j = 1;
    else
        j = tri2dec(parents, size);

    // draw an expression level from p(x | parents)
    // TODO replace generation based on frequency
    const double *prob_vector = p_1[i][j];
    int freq[3];
    for (int i = 0; i < 3; ++i)
        freq[i] = prob_vector[i] * 100;

    int categories[] =
        { 1, 2, 3 };
    int n = sizeof(categories) / sizeof(categories[0]);
    srand(time(NULL));
    int e = my_rand(categories, freq, n);

    return e;
}

/* Function
 *
 * converts a trinary # represented as an array into a decimal index
 */
int
tri2dec(int x[], int n)
{
    int sum = 0;
    for (int i = 0; i < n; ++i)
    {
        int exp = 3;
        if (i == 0)
            exp = 1;
        else
            for (int j = 1; j < i; ++j)
                exp *= 3;

        sum += exp * x[i];
    }

    return sum;
}

/* Function
 * Source:
 * http://www.geeksforgeeks.org/random-number-generator-in-arbitrary-probability-distribution-fashion/
 *
 * generate random numbers based on frequency
 */
int
find_ceil(int arr[], int r, int l, int h)
{
    int mid;
    while (l < h)
    {
        mid = l + ((h - l) >> 1);  // Same as mid = (l+h)/2
        (r > arr[mid]) ? (l = mid + 1) : (h = mid);
    }
    return (arr[l] >= r) ? l : -1;
}

// The main function that returns a random number from arr[] according to
// distribution array defined by freq[]. n is size of arrays.
int
my_rand(int arr[], int freq[], int n)
{
    // Create and fill prefix array
    int prefix[n], i;
    prefix[0] = freq[0];
    for (i = 1; i < n; ++i)
        prefix[i] = prefix[i - 1] + freq[i];

    // prefix[n-1] is sum of all frequencies. Generate a random number
    // with value from 1 to this sum
    int r = (rand() % prefix[n - 1]) + 1;

    // Find index of ceiling of r in prefix array
    int indexc = find_ceil(prefix, r, 0, n - 1);
    return arr[indexc];
}
