#define _POSIX_C_SOURCE 201710L

#include "constraint.h"
#include "misc.h"

#include <stdlib.h>
#include <math.h>

struct FeasibilityResult {
	bool feasible;
	double *certificate;
	size_t certificate_length;
};

void eliminate_variable(struct LP *lin_prog);

struct FeasibilityResult fourier_motzkin_destructive(struct LP *lin_prog);
double feasible_last_variable(
	struct Constraint *constraints,
	size_t const num_constraints,
	double const* smaller_solution,
	size_t const num_variables
);

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
	lp_print(&lin_prog);
	printf("\n");
	*/

	struct FeasibilityResult result = fourier_motzkin_destructive(&lin_prog);
	if (result.feasible) {
		for (size_t var=0; var<result.certificate_length; ++var) {
			if (var > 0)
				printf(" ");
			printf("%g", result.certificate[var]);
		}
		printf("\n");
	} else {
		printf("INFEASIBLE\n");
	}
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
	lin_prog->num_variables--;

}

struct FeasibilityResult fourier_motzkin_destructive(struct LP *lin_prog)
{
	if (lin_prog->num_variables == 0) {
		struct FeasibilityResult result = {.feasible = true, .certificate = NULL, .certificate_length = 0};
		/* Look for a constraint of the form "0 <= b", b<0 */
		for (size_t c=0; c<lin_prog->num_constraints; ++c)
			if (lin_prog->constraints[c].value < 0)
				result.feasible = false;
		return result;
	} else {
		//fprintf(stderr, "Reducing from %lu to %lu variables. ", lin_prog->num_variables, lin_prog->num_variables-1);
		struct Constraint *old_constraints = lin_prog->constraints;
		size_t const old_num_constraints = lin_prog->num_constraints;
		size_t const old_num_variables = lin_prog->num_variables;
		eliminate_variable(lin_prog);
		//fprintf(stderr, "Resulting number of constraints: %lu\n", lin_prog->num_constraints);
		struct FeasibilityResult result = fourier_motzkin_destructive(lin_prog);
		if (result.feasible) {
			result.certificate_length = old_num_variables;
			result.certificate = realloc(result.certificate, old_num_variables * sizeof(double));
			result.certificate[old_num_variables - 1] = feasible_last_variable(
				old_constraints,
				old_num_constraints,
				result.certificate,
				old_num_variables
			);
		}
		/* Free old constraints */
		for (size_t c=0; c<old_num_constraints; ++c)
			if (old_constraints[c].type != EQUAL)
				constraint_free(old_constraints + c);
		free(old_constraints);
		return result;
	}
}

double feasible_last_variable(
	struct Constraint *constraints,
	size_t const num_constraints,
	double const* smaller_solution,
	size_t const num_variables
)
{
	double lower_bound = -INFINITY;
	double upper_bound = INFINITY;
	for (size_t c=0; c<num_constraints; ++c) {
		double const last_coefficient = constraints[c].linear_combination[num_variables-1];
		if (last_coefficient == 0)
			continue;
		if (constraints[c].type != LESS_EQUAL)
			error(1, "feasible_last_variable does not (yet) support different constraint types than \"<=\"");
		double value = constraints[c].value;
		for (size_t var=0; var<num_variables-1; ++var)
			value -= constraints[c].linear_combination[var] * smaller_solution[var];
		value /= last_coefficient;
		if ((last_coefficient > 0) && (value < upper_bound)) {
			upper_bound = value;
		} else if ((last_coefficient < 0) && (value > lower_bound)) {
			lower_bound = value;
		}
	}
	if (lower_bound > upper_bound)
		error(1, "feasible_last_variable cannot return valid value");
	if (lower_bound != -INFINITY)
		return lower_bound;
	else if (upper_bound != INFINITY)
		return upper_bound;
	else
		return 0;
}
