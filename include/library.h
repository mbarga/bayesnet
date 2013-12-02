#ifndef _LIBRARY_H_
#define _LIBRARY_H_

#define Calloc(type,n) calloc(n, sizeof(type))
#define matrix(a,dim,a_i,a_j) a[a_i * dim + a_j] // matrix (row,col) -> (a_i,a_j)

typedef enum _config_error
{
  E_SUCCESS = 0,
  E_INVALID_INPUT = 1,
  E_FILE_ERROR = 2,
  E_NULL_POINTER = 3,
  E_NOT_FOUND = 4,
  E_DUPLICATE = 5,
  E_LAST_ERROR
} _CONFIG_ERROR;

/*
struct _errordesc 
{
  int code;
  char *message;
} errordesc[] = {
  { E_SUCCESS, "no error" },
  { E_INVALID_INPUT, "INVALID INPUT" },
  { E_FILE_ERROR, "FILE INPUT ERROR" }
};
*/

void
util_print_dmatrix(double *, int, int);
void
util_print_imatrix(int *, int);
void
util_errlog(char *);

#endif
