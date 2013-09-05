#ifndef _READFILE_H
#define _READFILE_H

#include <stdio.h>
#include <stdlib.h>
#include <syslog.h>
#include <glib.h>
#include <math.h>

int read_problem(const char *filename, GHashTable **x);
NODE* create_node();

#endif
