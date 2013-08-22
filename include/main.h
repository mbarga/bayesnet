#ifndef _MAIN_H
#define _MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <string.h>
#include <time.h>

#define MAX_ITER 10000
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static char *line = NULL;
static int max_line_len;
static int max_name_len = 128;

typedef struct node
{
	char* 	name;		// node name as string
	double* data;		// node probability data array
} NODE;

int read_problem(const char *filename, NODE **x, int *n);
static char* read_line(FILE *input);

#endif
