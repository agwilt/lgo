#define _POSIX_C_SOURCE 201710L

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>

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

struct LP {
	double *objective;
	struct Constraint *constraints;
	size_t num_constraints;
	size_t _max_num_constraints;
	size_t num_variables;
};

/* return value: does variable have nonzero coefficient? */
bool constraint_normalise_variable(
	struct Constraint *constraint,
	size_t const num_variables,
	size_t const variable
);

void constraint_multiply(
	struct Constraint *constraint,
	size_t const num_variables,
	double const factor
);

void constraint_free(struct Constraint *constraint);

/* create a new constraint that is the sum of the two given constraints */
struct Constraint constraint_sum(
	struct Constraint *lhs,
	struct Constraint *rhs,
	size_t const num_variables
);

struct Constraint constraint_clone(struct Constraint *orig, size_t const num_variables);

struct Constraint create_constraint_empty(
	enum ConstraintType type,
	size_t const num_variables
);

void lp_print(struct LP *lin_prog);
void lp_print_human_readable(struct LP *lin_prog);

void lp_free(struct LP *linear_program);
void lp_add_constraint(struct LP *linear_program, struct Constraint const constraint);
void lp_remove_constraint(struct LP *lin_prog, size_t index);
struct LP create_lp_empty(size_t const num_variables, size_t const max_num_constraints);
struct LP create_lp_from_file(FILE *fp);
void lp_prune(struct LP *lin_prog);

#endif // CONSTRAINT_H
