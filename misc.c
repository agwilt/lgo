#define _POSIX_C_SOURCE 201710L

#include "misc.h"

#include <stdio.h>
#include <stdlib.h>

void error(int return_code, char *string)
{
	fprintf(stderr, "ERROR: %s\n", string);
	exit(return_code);
}
