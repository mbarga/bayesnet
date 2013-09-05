#include "main.h"
#include "readfile.h"

static char *line = NULL;
static int max_line_len;
//static int max_name_len = 128;

static char * read_line(FILE *input)
{
	int len;

	if (fgets(line, max_line_len, input) == NULL)
	{
		return NULL;
	}

	while (strrchr(line, '\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line, max_line_len);
		len = (int) strlen(line);
		if (fgets(line + len, max_line_len - len, input) == NULL)
			break;
	}

	return line;
}

int read_problem(const char *filename, GHashTable **x)
{
	FILE *fp = fopen(filename, "r");
	//int isspace(int c);
	//char *endptr;
	char *child_name, *parent_name;

	GHashTable* hash = g_hash_table_new(g_str_hash, g_str_equal);
	NODE *parent = NULL;
	NODE *child = NULL;

	int node_count = 0;

	if (fp == NULL)
	{
		syslog(LOG_INFO, "%s", "read_problem(): couldnt open input file\n");
		return 1;
	}

	max_line_len = 1024;
	line = Malloc(char, max_line_len);

	while (read_line(fp) != NULL)
	{
		child_name = strtok(line, " \t\n");
		parent_name	= strtok(NULL, "\t\n"); //TODO clean this up??
	//	child_name	= strtok(NULL, &endptr);
		if (parent_name == NULL || child_name == NULL)
		{
			syslog(LOG_INFO, "%s", "FAILED::read_problem : empty line\n");
			return 1;
		}
	
		child = g_hash_table_lookup(hash, child_name);
		if (child == NULL)
		{
			child 				= create_node();
			child->index	= node_count++;
			child->name 	= g_strdup(child_name);
 			g_hash_table_insert(hash, g_strdup(child_name), child); //TODO use node->name
		}
		
		parent = g_hash_table_lookup(hash, parent_name);
		if (parent == NULL)
		{
			parent 				= create_node();
			parent->index	= node_count++;
			parent->name 	= g_strdup(parent_name);
 			g_hash_table_insert(hash, g_strdup(parent_name), parent);
		}
	
		parent->edges = g_slist_append(parent->edges, child);
	}

	*x = hash;

	free(line);
	child_name = NULL;
	parent_name = NULL;
	parent = NULL;
	child = NULL;
	fclose(fp);
	return 0;
}

NODE * create_node()
{
	NODE *p;
	//p = Malloc(NODE, 1);
	p 				= g_new(NODE, 1);
	p->name 	= NULL; 
	p->index 	= 0;
	p->edges	= NULL;
	return p;
}
