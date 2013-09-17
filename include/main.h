#ifndef _MAIN_H_
#define _MAIN_H_

#include <glib.h>

// TODO pass in these statics as values to the program?
#define MAX_ITER 40
#define MAX_SUBSETS 10
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

typedef struct node
{
	gchar* 	name;		// node name as string
	gint 		index;	// node probability data array
	int			expression; // C {-1, 0, 1} expression level of gene
	GSList* parents;
} NODE;

void estimate_dag(NODE*, int, int *, int, int**);

#endif
