#define _POSIX_C_SOURCE 201710L

#include "constraint.h"
#include "misc.h"

#include <stdlib.h>

void eliminate_variable(struct LP *lin_prog);

bool fourier_motzkin_destructive(struct LP *lin_prog);

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
	printf("\n");

	if (fourier_motzkin_destructive(&lin_prog))
		printf("\nFEASIBLE\n\n");
	else
		printf("\nINFEASIBLE\n\n");
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

	/* Normalise constraints */
	for (size_t c=0; c<old_num_constraints; ++c)
		constraint_normalise_variable(old_constraints + c, lin_prog->num_variables, v);
	//lp_prune(lin_prog);

	/* Append (sums of) relevant constraints */
	for (size_t c=0; c<old_num_constraints; ++c) {
		if (old_constraints[c].type == EQUAL)
			lp_add_constraint(lin_prog, old_constraints[c]);
		else if (old_constraints[c].linear_combination[v] == 1) {
			for (size_t d=0; d<old_num_constraints; ++d) {
				if (old_constraints[d].linear_combination[v] == -1) {
					lp_add_constraint(lin_prog, constraint_sum(old_constraints + c,
					                                           old_constraints + d,
					                                           lin_prog->num_variables - 1)
					);
				}
			}
		}
	}

	/* Free old constraints */
	for (size_t c=0; c<old_num_constraints; ++c)
		if (old_constraints[c].type != EQUAL)
			constraint_free(old_constraints + c);

	lin_prog->num_variables--;

	free(old_constraints);
}

bool fourier_motzkin_destructive(struct LP *lin_prog)
{
	while (lin_prog->num_variables > 0) {
		printf("Reducing from %lu to %lu variables. ", lin_prog->num_variables, lin_prog->num_variables-1);
		eliminate_variable(lin_prog);
		printf("Resulting number of constraints: %lu\n", lin_prog->num_constraints);
	}

	/* Look for a constraint of the form "0 <= b", b<0 */
	for (size_t c=0; c<lin_prog->num_constraints; ++c)
		if (lin_prog->constraints[c].value < 0)
			return false;

	/* Nothing found: must be feasible :) */
	return true;
}
