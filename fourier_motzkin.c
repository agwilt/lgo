#define _POSIX_C_SOURCE 201710L

#include "constraint.h"
#include "lp.h"
#include "misc.h"

#include <stdlib.h>
#include <math.h>

struct FeasibilityResult {
	bool feasible;
	double *certificate;
	size_t certificate_length;
};

struct LinCombTuple {
	size_t index1;
	size_t index2;
	double factor1;
	double factor2;
};

struct LP eliminate_variable(struct LP const* lin_prog, struct LinCombTuple **constr_lin_comb);

struct FeasibilityResult fourier_motzkin(struct LP const* lin_prog);

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
	if (! result.feasible)
		printf("empty ");
	/* Print certificate */
	for (size_t i=0; i<result.certificate_length; ++i) {
		if (i > 0)
			printf(" ");
		printf("%g", result.certificate[i]);
	}
	printf("\n");

	if (result.feasible) {
		/* Check solution */
		if (result.certificate_length != lin_prog.num_variables)
			error(1, "Output broken: Certificate has wrong dimension.");
		for (size_t c = 0; c < lin_prog.num_constraints; ++c)
			if (! constraint_fulfilled(lin_prog.constraints + c, result.certificate, lin_prog.num_variables))
				error(1, "Output broken: Certificate invalid.");
	} else {
		/* Check linear combination */
		if (result.certificate_length != lin_prog.num_constraints)
			error(1, "Output broken: Certificate has wrong dimension.");
		for (size_t c = 0; c < result.certificate_length; ++c)
			if (result.certificate[c] < 0)
				error(1, "Output broken: Certificate has negative components.");
		double scalar_prod = 0;
		for (size_t c = 0; c < result.certificate_length; ++c)
			scalar_prod += result.certificate[c] * lin_prog.constraints[c].value;
		if (scalar_prod >= 0)
			error(1, "Output broken: Certificate has non-negative scalar product with b.");
		for (size_t var = 0; var < lin_prog.num_variables; ++var) {
			double entry = 0.0;
			for (size_t c = 0; c < lin_prog.num_constraints; ++c)
				entry += result.certificate[c] * lin_prog.constraints[c].linear_combination[var];
			if (entry != 0)
				error(1, "Output broken: Certificate not orthogonal to im(A).");
		}
	}

	lp_free(&lin_prog);
	free(result.certificate);

	return 0;
}

/* NOTE: constr_lin_comb will be alloc'ed and will have (new_lp.num_constraints) entries */
struct LP eliminate_variable(struct LP const* lin_prog, struct LinCombTuple **constr_lin_comb)
{
	size_t const last_var = lin_prog->num_variables - 1;

	/* Count the number of new constraints; I can't be bothered to implement something like push_back */
	size_t new_num_constraints = 0;
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		double const last_coefficient = lin_prog->constraints[c].linear_combination[last_var];
		if (last_coefficient == 0.0)
			++new_num_constraints;
		else if (last_coefficient > 0)
			for (size_t d=0; d<lin_prog->num_constraints; ++d)
				new_num_constraints += (lin_prog->constraints[d].linear_combination[last_var] < 0);
	}

	struct LP new_lp = create_lp_empty(lin_prog->num_variables - 1, new_num_constraints);
	*constr_lin_comb = realloc(*constr_lin_comb, new_num_constraints * sizeof(struct LinCombTuple));

	size_t new_constraints_added = 0;
	/* Append (sums of) relevant constraints */
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		if (lin_prog->constraints[c].type != LESS_EQUAL)
			error(1, "eliminate_variable only supports \"<\" constraints for now");
		double const last_coefficient = lin_prog->constraints[c].linear_combination[last_var];
		if (last_coefficient == 0.0) {
			lp_add_constraint(&new_lp, constraint_clone(lin_prog->constraints + c, new_lp.num_variables));
			(*constr_lin_comb)[new_constraints_added++] = (struct LinCombTuple) {
				.index1 = c,
				.index2 = c,
				.factor1 = 1.0,
				.factor2 = 0.0
			};
		} else if (last_coefficient > 0) {
			for (size_t d=0; d<lin_prog->num_constraints; ++d) {
				double const other_last_coefficient = lin_prog->constraints[d].linear_combination[last_var];
				if (other_last_coefficient < 0) {
					lp_add_constraint(
						&new_lp,
						constraint_sum(lin_prog->constraints + c,
					                       lin_prog->constraints + d,
					                       new_lp.num_variables,
							       -other_last_coefficient,
							       last_coefficient)
					);
					(*constr_lin_comb)[new_constraints_added++] = (struct LinCombTuple) {
						.index1 = c,
						.index2 = d,
						.factor1 = -other_last_coefficient,
						.factor2 = last_coefficient
					};
				}
			}
		}
	}
	return new_lp;
}

struct FeasibilityResult fourier_motzkin(struct LP const* lin_prog)
{
	if (lin_prog->num_variables == 0) {
		struct FeasibilityResult result = {.feasible = true, .certificate = NULL, .certificate_length = 0};
		/* Look for a constraint of the form "0 <= b", b<0 */
		for (size_t c=0; c<lin_prog->num_constraints; ++c) {
			if (lin_prog->constraints[c].value < 0) {
				result.feasible = false;
				result.certificate_length = lin_prog->num_constraints;
				result.certificate = calloc(result.certificate_length, sizeof(double));
				result.certificate[c] = 1.0;
				break;
			}
		}
		return result;
	} else {
		struct LinCombTuple *lin_combs = NULL;
		struct LP new_lp = eliminate_variable(lin_prog, &lin_combs);
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
		} else {
			double *infeasibility_cert = result.certificate;
			result.certificate_length = lin_prog->num_constraints;
			result.certificate = calloc(lin_prog->num_constraints, sizeof(double));
			for (size_t c=0; c<new_lp.num_constraints; ++c) {
				result.certificate[lin_combs[c].index1] += lin_combs[c].factor1 * infeasibility_cert[c];
				result.certificate[lin_combs[c].index2] += lin_combs[c].factor2 * infeasibility_cert[c];
			}
			free(infeasibility_cert);
		}
		/* Free new LP and linear combination data */
		lp_free(&new_lp);
		free(lin_combs);
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
