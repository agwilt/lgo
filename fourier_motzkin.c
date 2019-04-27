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

struct LP eliminate_variable(struct LP *lin_prog);

struct FeasibilityResult fourier_motzkin(struct LP *lin_prog);

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
	lp_print_nice(&lin_prog);
	printf("\n");
	*/

	struct FeasibilityResult result = fourier_motzkin(&lin_prog);
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
	free(result.certificate);

	return 0;
}

struct LP eliminate_variable(struct LP *lin_prog)
{
	struct LP new_lp = create_lp_empty(lin_prog->num_variables - 1, lin_prog->_max_num_constraints);
	size_t const last_var = lin_prog->num_variables - 1;

	//lp_normalise_variable(lin_prog, last_var);
	//lp_prune(lin_prog);

	/* Append (sums of) relevant constraints */
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		if (lin_prog->constraints[c].type != LESS_EQUAL)
			error(1, "eliminate_variable only supports \"<\" constraints for now");
		double const last_coefficient = lin_prog->constraints[c].linear_combination[last_var];
		if (last_coefficient == 0.0)
			lp_add_constraint(&new_lp, constraint_clone(lin_prog->constraints + c, new_lp.num_variables));
		else if (last_coefficient > 0)
			for (size_t d=0; d<lin_prog->num_constraints; ++d)
				if (lin_prog->constraints[d].linear_combination[last_var] < 0)
					lp_add_constraint(
						&new_lp,
						constraint_sum(lin_prog->constraints + c,
					                       lin_prog->constraints + d,
					                       new_lp.num_variables,
							       -lin_prog->constraints[d].linear_combination[last_var],
							       last_coefficient)
					);
	}
	return new_lp;
}

struct FeasibilityResult fourier_motzkin(struct LP *lin_prog)
{
	if (lin_prog->num_variables == 0) {
		struct FeasibilityResult result = {.feasible = true, .certificate = NULL, .certificate_length = 0};
		/* Look for a constraint of the form "0 <= b", b<0 */
		for (size_t c=0; c<lin_prog->num_constraints; ++c)
			if (lin_prog->constraints[c].value < 0)
				result.feasible = false;
		return result;
	} else {
		fprintf(stderr, "Reducing from %lu to %lu variables. ", lin_prog->num_variables, lin_prog->num_variables-1);
		struct LP new_lp = eliminate_variable(lin_prog);
		fprintf(stderr, "Resulting number of constraints: %lu\n", lin_prog->num_constraints);
		struct FeasibilityResult result = fourier_motzkin(&new_lp);
		if (result.feasible) {
			result.certificate_length = lin_prog->num_variables;
			result.certificate = realloc(result.certificate, lin_prog->num_variables * sizeof(double));
			result.certificate[lin_prog->num_variables - 1] = feasible_last_variable(
				lin_prog->constraints,
				lin_prog->num_constraints,
				result.certificate,
				lin_prog->num_variables
			);
		}
		/* Free new LP */
		lp_free(&new_lp);
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
	fprintf(stderr, "\nFinding value for x_%lu\n", num_variables);
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
		fprintf(stderr, "x_%lu in [%g, %g]\n", num_variables, lower_bound, upper_bound);
	}
	if (lower_bound > upper_bound) {
		fprintf(stderr, "ERROR: Feasible_last_variable cannot return valid value in [%g, %g]\n", lower_bound, upper_bound);
		exit(1);
	}
	if ((lower_bound <= 0) && (upper_bound >= 0))
		return 0;
	else if ((lower_bound != -INFINITY) && (upper_bound != 	INFINITY))
		return (upper_bound + lower_bound) / 2.0;
	else if ((lower_bound == -INFINITY) && (upper_bound != INFINITY))
		return round(upper_bound - 1);
	else if ((upper_bound == INFINITY) && (lower_bound != INFINITY))
		return round(lower_bound + 1);
	else {
		error(1, "Something just went horribly wrong.");
		return NAN;
	}
}
