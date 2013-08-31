#ifndef _MAIN_H
#define _MAIN_H

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <string.h>
#include <time.h>
#include <glib.h> //TODO will this be supported on all machines?

#define MAX_ITER 10000
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

static char *line = NULL;
static int max_line_len;
//static int max_name_len = 128;

typedef struct node
{
	gchar* 	name;		// node name as string
	gint 		index;	// node probability data array
	GSList* edges;
} NODE;

int read_problem(const char *filename, GHashTable **x);
static char* read_line(FILE *input);
NODE* create_node();

#endif
