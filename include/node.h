#ifndef _NODE_H_
#define _NODE_H_

typedef struct node
{
	//gchar* name;        // node name as string
	//gint index;         // node probability data array
	//int expression;     // C {-1, 0, 1} expression level of gene

	//GSList* parents;
	int index;
	int *parents;
	int num_parents;
	double score;
} NODE;

#endif
