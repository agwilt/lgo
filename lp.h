#define _POSIX_C_SOURCE 201710L

#ifndef LP_H
#define LP_H

#include <stddef.h>
#include <stdio.h>

#include "constraint.h"

/*
 * Data Structures
 */

struct LP {
	double *objective;
	struct Constraint *constraints;
	size_t num_constraints;
	size_t _max_num_constraints;
	size_t num_variables;
};

/*
 * Functions
 */

struct LP create_lp_empty(size_t num_variables, size_t max_num_constraints);
struct LP create_lp_from_file(FILE *fp);

void lp_free(struct LP *linear_program);

void lp_print(struct LP const* lin_prog);
void lp_print_nice(struct LP const* lin_prog);

void lp_add_constraint(struct LP *linear_program, struct Constraint constraint);
void lp_remove_constraint(struct LP *lin_prog, size_t index);

#endif // LP_H
