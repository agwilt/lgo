#define _POSIX_C_SOURCE 201710L

#include "constraint.h"
#include "misc.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/* Constraint Stuff */

struct Constraint create_constraint_empty(enum ConstraintType const type, size_t const num_variables)
{
	return (struct Constraint) {
		.linear_combination = calloc(num_variables, sizeof(double)),
		.value = 0,
		.type = type,
	};
}

struct Constraint constraint_clone(struct Constraint const* orig, size_t const num_variables)
{
	struct Constraint clone = {
		.linear_combination = malloc(num_variables * sizeof(double)),
		.value = orig->value,
		.type = orig->type,
	};
	memcpy(clone.linear_combination, orig->linear_combination, num_variables*sizeof(double));
	return clone;
}

void constraint_free(struct Constraint *constraint)
{
	free(constraint->linear_combination);
	constraint->linear_combination = NULL;
}

struct Constraint constraint_sum(
	struct Constraint const* lhs,
	struct Constraint const* rhs,
	size_t const num_variables,
	double const factor_lhs,
	double const factor_rhs
)
{
	if (lhs->type != rhs->type)
		error(1, "Tried to add two incompatible constraints");
	if ((factor_lhs < 0) || (factor_rhs < 0))
		error(1, "Invalid factors for constraint_sum");
	struct Constraint sum = {
		.linear_combination = malloc(num_variables * sizeof(double)),
		.value = factor_lhs * lhs->value + factor_rhs * rhs->value,
		.type = lhs->type,
	};
	for (size_t i=0; i<num_variables; ++i)
		sum.linear_combination[i] = factor_lhs * lhs->linear_combination[i] + factor_rhs * rhs->linear_combination[i];
	return sum;
}

bool constraint_fulfilled(struct Constraint const* constraint, double const* assignment, size_t const num_variables)
{
	double lhs = 0;
	for (size_t var = 0; var < num_variables; ++var) {
		lhs += assignment[var] * constraint->linear_combination[var];
	}
	switch (constraint->type) {
		case LESS_EQUAL:
			return lhs <= constraint->value;
		case GREATER_EQUAL:
			return lhs >= constraint->value;
		case EQUAL:
			return lhs == constraint->value;
		default:
			fprintf(stderr, "Error in constraint_fulfilled: Invalid constraint");
			exit(1);
	}
}
