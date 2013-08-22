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

double score(NODE *u, NODE *v)
{
	double score = 0;

	return score;
}
