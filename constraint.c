#define _POSIX_C_SOURCE 201710L

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "constraint.h"
#include "misc.h"

bool constraint_normalise_variable(struct Constraint *constraint, size_t const variable)
{
	double denominator = fabs(constraint->linear_combination[variable]);
	if (denominator == 0) {
		return false;
	} else {
		for (size_t i=0; i<constraint->num_variables; ++i)
			constraint->linear_combination[i] /= denominator;
		constraint->value /= denominator;
		return true;
	}
}

void constraint_multiply(struct Constraint *constraint, double const factor)
{
	for (size_t i=0; i<constraint->num_variables; ++i)
		constraint->linear_combination[i] *= factor;
	constraint->value *= factor;
	if (factor < 0)
		constraint->type *= -1;
}

void constraint_free(struct Constraint *constraint)
{
	if (constraint->num_variables > 0)
		free(constraint->linear_combination);
	constraint->linear_combination = NULL;
	constraint->num_variables = 0;
}

struct Constraint constraint_sum(struct Constraint *lhs, struct Constraint *rhs)
{
	if ((lhs->type != rhs->type) || (lhs->num_variables != rhs->num_variables))
		error(1, "Tried to add two incompatible constraints");
	struct Constraint sum = {
		.linear_combination = calloc(lhs->num_variables, sizeof(double)),
		.value = lhs->value + rhs->value,
		.type = lhs->type,
		.num_variables = lhs->num_variables
	};
	for (size_t i=0; i<sum.num_variables; ++i)
		sum.linear_combination[i] = lhs->linear_combination[i] + rhs->linear_combination[i];
	return sum;
}

struct Constraint constraint_clone(struct Constraint *orig)
{
	struct Constraint clone = {
		.linear_combination = malloc(orig->num_variables * sizeof(double)),
		.value = orig->value,
		.type = orig->type,
		.num_variables = orig->num_variables
	};
	memcpy(clone.linear_combination, orig->linear_combination, orig->num_variables*sizeof(double));
	return clone;
}

struct Constraint create_constraint_empty(enum ConstraintType type, size_t num_variables)
{
	return (struct Constraint) {
		.linear_combination = calloc(num_variables, sizeof(double)),
		.value = 0,
		.type = type,
		.num_variables = num_variables
	};
}

void lp_print(struct LP *lin_prog)
{
	printf("%lu %lu\n", lin_prog->num_constraints, lin_prog->num_variables);
	// Print objective function
	for (size_t var=0; var<lin_prog->num_variables; ++var)
		printf((var==lin_prog->num_variables-1)?"%g\n":"%g ", lin_prog->objective[var]);
	// Print constraint upper bounds
	for (size_t c=0; c<lin_prog->num_constraints; ++c)
		printf((c==lin_prog->num_constraints-1)?"%g\n":"%g ", lin_prog->constraints[c].value);
	// Print constraint linear combinations
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		if (lin_prog->constraints[c].type != LESS_EQUAL)
			error(1, "LP printing does not support this constraint type.");
		for (size_t var=0; var<lin_prog->num_variables; ++var) {
			printf((var==lin_prog->num_variables-1)?"%g\n":"%g ",
			       lin_prog->constraints[c].linear_combination[var]);
		}
	}
}

void lp_free(struct LP *linear_program)
{
	if (linear_program->num_variables > 0)
		free(linear_program->objective);

	for (size_t c=0; c<linear_program->num_constraints; ++c)
		constraint_free(linear_program->constraints + c);
	if (linear_program->_max_num_constraints > 0)
		free(linear_program->constraints);

	linear_program->objective = NULL;
	linear_program->constraints = NULL;
	linear_program->num_constraints = 0;
	linear_program->_max_num_constraints = 0;
	linear_program->num_variables = 0;
}

void lp_add_constraint(struct LP *linear_program, struct Constraint const constraint)
{
	if (linear_program->_max_num_constraints == 0) {
		linear_program->constraints = malloc(sizeof(struct Constraint));
		linear_program->_max_num_constraints = 1;
	} else while (linear_program->_max_num_constraints <= linear_program->num_constraints) {
		linear_program->_max_num_constraints <<= 1;
		linear_program->constraints = realloc(linear_program->constraints, 
			linear_program->_max_num_constraints * sizeof(struct Constraint));
	}
	linear_program->constraints[linear_program->num_constraints++] = constraint;
}

struct LP create_lp_empty(size_t const num_variables, size_t const max_num_constraints)
{
	return (struct LP) {
		.objective = calloc(num_variables, sizeof(double)),
		.constraints = malloc(max_num_constraints * sizeof(struct Constraint)),
		.num_constraints = 0,
		._max_num_constraints = max_num_constraints,
		.num_variables = num_variables
	};
}

struct LP create_lp_from_file(FILE *fp)
{
	size_t num_variables, num_constraints;

	if (fscanf(fp, "%lu %lu\n", &num_constraints, &num_variables) != 2)
		error(1, "Invalid file format.");

	struct LP lin_prog = create_lp_empty(num_variables, num_constraints);

	// Read objective function
	for (size_t var = 0; var < num_variables; ++var)
		if (fscanf(fp, (var==num_variables-1) ? "%lf\n" : "%lf ", lin_prog.objective+var) != 1)
			error(1, "Cannot read objective function.");

	// Read constraint upper-bounds, create constraints
	lin_prog.num_constraints = num_constraints;
	if (lin_prog.num_constraints > lin_prog._max_num_constraints)
		error(1, "Something just went terribly wrong.");
	for (size_t c = 0; c < num_constraints; ++c) {
		lin_prog.constraints[c].linear_combination = malloc(num_variables * sizeof(double));
		if (fscanf(fp,
		           (c==num_constraints-1) ? "%lf\n" : "%lf ",
		           &(lin_prog.constraints[c].value)) != 1) {
			error(1, "Cannot read constraint bounds.");
		}
		lin_prog.constraints[c].type = LESS_EQUAL;
		lin_prog.constraints[c].num_variables = num_variables;
	}
	printf("Have %lu constraints and %lu variables.\n", lin_prog.num_constraints, lin_prog.num_variables);

	// Read constraint linear combinations
	char *line = NULL;
	size_t len = 0;
	for (size_t c = 0; c < num_constraints; ++c) {
		if (getline(&line, &len, fp) == -1)
			error(1, "Incorrect number of constraints");
		char *p = line;
		for (size_t var = 0; var < num_variables; ++var)
			lin_prog.constraints[c].linear_combination[var] = strtod(p, &p);
	}
	if (line)
		free(line);

	return lin_prog;
}
