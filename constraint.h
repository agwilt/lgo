#define _POSIX_C_SOURCE 201710L

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <stddef.h>
#include <stdbool.h>

/*
 * Data Structures
 */

enum ConstraintType {
	EQUAL = 0,
	GREATER_EQUAL = 1,
	LESS_EQUAL = -1
};

struct Constraint {
	double *linear_combination;
	double value;
	enum ConstraintType type;
};

/*
 * Functions
 */

struct Constraint create_constraint_empty(
	enum ConstraintType type,
	size_t num_variables
);

struct Constraint constraint_clone(struct Constraint const* orig, size_t num_variables);

void constraint_free(struct Constraint *constraint);

struct Constraint constraint_sum(
	struct Constraint const* lhs,
	struct Constraint const* rhs,
	size_t num_variables,
	double factor_lhs,
	double factor_rhs
);

bool constraint_fulfilled(struct Constraint const* constraint, double const* assignment, size_t num_variables);

#endif // CONSTRAINT_H
