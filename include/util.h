#ifndef _UTIL_H_
#define _UTIL_H_

#define Calloc(type,n) calloc(n, sizeof(type))
#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
#define matrix(a,dim,a_i,a_j) a[a_i * dim + a_j] // matrix (row,col) -> (a_i,a_j)

#define INF 999999999 //TODO rewrite this

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

double get_score(void *, int, int *, int);
void util_print_dmatrix(double *, int, int);
void util_print_imatrix(int *, int, const char *);
void util_debuglog(char *);
void util_errlog(char *);
void util_print_score_table(void *);

#endif
