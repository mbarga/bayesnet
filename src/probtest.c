#include "probability.h"

#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char * argv[])
{
    int parents[2] = {0, 0}; // states of the parents TODO in what order??

    int expression = get_expression(parents, 2);

    printf("expression was: %d\n", expression);

    return 0;
}



