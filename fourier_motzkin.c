#define _POSIX_C_SOURCE 201710L

#include "constraint.h"
#include "misc.h"

#include <stdlib.h>

void eliminate_variable(struct LP *lin_prog);

int main(int argc, char *argv[])
{
	if (argc < 2)
		error(1, "No file name supplied.");
	FILE *fp = fopen(argv[1], "r");
	if (fp == NULL) {
		fprintf(stderr, "ERROR: Cannot open file: %s\n", argv[1]);
		return 1;
	}
	struct LP lin_prog = create_lp_from_file(fp);
	fclose(fp);

	/*
	struct LP lin_prog = create_lp_empty(3, 1);

	lin_prog.objective[0] = -2.0;
	lin_prog.objective[1] = 0.0;
	lin_prog.objective[2] = 8.0;

	struct Constraint c_1 = create_constraint_empty(LESS_EQUAL, 3);
	c_1.linear_combination[0] = 3.5;
	c_1.linear_combination[1] = -2.0;
	c_1.linear_combination[2] = 5.0;
	c_1.value = 3;
	lp_add_constraint(&lin_prog, c_1);

	struct Constraint c_2 = create_constraint_empty(LESS_EQUAL, 3);
	c_2.linear_combination[0] = 0.0;
	c_2.linear_combination[1] = 1.0;
	c_2.linear_combination[2] = -4.0;
	c_2.value = 0;
	lp_add_constraint(&lin_prog, c_2);
	*/

	lp_print(&lin_prog);
	lp_free(&lin_prog);

	return 0;
}

void eliminate_variable(struct LP *lin_prog)
{
	size_t const v = lin_prog->num_variables - 1;

	struct Constraint *old_constraints = lin_prog->constraints;
	size_t const old_num_constraints = lin_prog->num_constraints;

	lin_prog->constraints = malloc(lin_prog->_max_num_constraints * sizeof(struct Constraint));
	lin_prog->num_constraints = 0;

	// Append constraints that stay unchanged, and normalise the rest
	for (size_t c=0; c<old_num_constraints; ++c) {
		if (constraint_normalise_variable(old_constraints + c, v))
			lp_add_constraint(lin_prog, old_constraints[c]);
	}

	// Append sums of relevant constraints
	for (size_t c=0; c<old_num_constraints; ++c) {
		if (old_constraints[c].type == EQUAL)
			continue;
		for (size_t d=0; d<old_num_constraints; ++d)
			if (old_constraints[d].type != EQUAL)
				lp_add_constraint(lin_prog, constraint_sum(old_constraints + c, old_constraints + d));
	}

	// Remove last variable
	lin_prog->num_variables--;
	for (size_t c=0; c<lin_prog->num_constraints; ++c)
		lin_prog->constraints[c].num_variables--;

	free(old_constraints);
}
