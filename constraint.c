#define _POSIX_C_SOURCE 201710L

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "constraint.h"
#include "misc.h"

bool constraint_implies(struct Constraint *lhs, struct Constraint *rhs, size_t const num_variables);

void constraint_normalise_variable(
	struct Constraint *constraint,
	size_t const num_variables,
	size_t const variable
)
{
	if (variable >= num_variables)
		error(1, "Bad call of constraint_normalise_variable");
	double denominator = fabs(constraint->linear_combination[variable]);
	if (denominator != 0) {
		constraint_multiply(constraint, num_variables, 1/denominator);
	}
}

void constraint_multiply(struct Constraint *constraint, size_t const num_variables, double const factor)
{
	for (size_t i=0; i<num_variables; ++i)
		constraint->linear_combination[i] *= factor;
	constraint->value *= factor;
	if (factor < 0)
		constraint->type *= -1;
}

void constraint_free(struct Constraint *constraint)
{
	free(constraint->linear_combination);
	constraint->linear_combination = NULL;
}

struct Constraint constraint_sum(
	struct Constraint *lhs,
	struct Constraint *rhs,
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

struct Constraint constraint_clone(struct Constraint *orig, size_t const num_variables)
{
	struct Constraint clone = {
		.linear_combination = malloc(num_variables * sizeof(double)),
		.value = orig->value,
		.type = orig->type,
	};
	memcpy(clone.linear_combination, orig->linear_combination, num_variables*sizeof(double));
	return clone;
}

struct Constraint create_constraint_empty(enum ConstraintType type, size_t num_variables)
{
	return (struct Constraint) {
		.linear_combination = calloc(num_variables, sizeof(double)),
		.value = 0,
		.type = type,
	};
}

bool constraint_implies(struct Constraint *lhs, struct Constraint *rhs, size_t const num_variables)
{
	for (size_t var=0; var < num_variables; ++var)
		if (lhs->linear_combination[var] != rhs->linear_combination[var])
			return false;
	return (((lhs->type==EQUAL) && (
			(lhs->value == rhs->value) ||
			((rhs->type==LESS_EQUAL) && (lhs->value < rhs->value)) ||
			((rhs->type==GREATER_EQUAL) && (lhs->value > rhs->value))
		)) ||
		((lhs->type==GREATER_EQUAL) && (rhs->type==GREATER_EQUAL) && !(lhs->value < rhs->value)) ||
		((lhs->type==LESS_EQUAL) && (rhs->type==LESS_EQUAL) && !(lhs->value > rhs->value)));
}

void lp_print(struct LP *lin_prog)
{
	printf("%lu %lu\n", lin_prog->num_constraints, lin_prog->num_variables);
	/* Print objective function */
	for (size_t var=0; var<lin_prog->num_variables; ++var)
		printf((var==lin_prog->num_variables-1)?"%g\n":"%g ", lin_prog->objective[var]);
	/* Print constraint upper bounds */
	for (size_t c=0; c<lin_prog->num_constraints; ++c)
		printf((c==lin_prog->num_constraints-1)?"%g\n":"%g ", lin_prog->constraints[c].value);
	/* Print constraint linear combinations */
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		if (lin_prog->constraints[c].type != LESS_EQUAL)
			error(1, "LP printing does not support this constraint type.");
		for (size_t var=0; var<lin_prog->num_variables; ++var) {
			printf((var==lin_prog->num_variables-1)?"%g\n":"%g ",
			       lin_prog->constraints[c].linear_combination[var]);
		}
	}
}

void lp_print_nice(struct LP *lin_prog)
{
#if 0
	printf("Varaibles: ");
	switch (lin_prog->num_variables) {
		case 0:
			printf("none");
			break;
		case 1:
			printf("x_1");
			break;
		case 2:
			printf("x_1, x_2");
			break;
		default:
			printf("x_1,...,x_%lu", lin_prog->num_variables);
	}
	printf("\n");
#endif
	/* Numbers */
	printf("Number of variables: %lu\n", lin_prog->num_variables);
	printf("Number of constraints: %lu\n", lin_prog->num_constraints);
	/* Objective function */
	printf("\nmax     ");
	for (size_t var=0; var<lin_prog->num_variables; ++var) {
		if (lin_prog->objective[var])
			printf("(%g * x_%lu)", lin_prog->objective[var], var+1);
		else
			printf("           ");
		if (var < lin_prog->num_variables - 1)
			printf(" + ");
	}
	printf("\n");
	/* Constraints */
	printf("s.t.    ");
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		struct Constraint *constraint = lin_prog->constraints + c;
		if (c > 0)
			printf("        ");
		for (size_t var=0; var<lin_prog->num_variables; ++var) {
			if (constraint->linear_combination[var])
				printf("(%g * x_%lu)", constraint->linear_combination[var], var+1);
			else
				printf("           ");
			if (var < lin_prog->num_variables - 1)
				printf(" + ");
		}
		if (constraint->type == EQUAL)
			printf(" = ");
		else if (constraint->type == GREATER_EQUAL)
			printf(" >= ");
		else if (constraint->type == LESS_EQUAL)
			printf(" <= ");
		else
			error(1, "Broken constraint while trying to print nicely");
		printf("%g\n", constraint->value);
	}
}

void lp_free(struct LP *linear_program)
{
	free(linear_program->objective);

	for (size_t c=0; c<linear_program->num_constraints; ++c)
		constraint_free(linear_program->constraints + c);
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

void lp_remove_constraint(struct LP *lin_prog, size_t index)
{
	constraint_free(lin_prog->constraints + index);
	if (index < lin_prog->num_constraints - 1) {
		lin_prog->constraints[index] = lin_prog->constraints[lin_prog->num_constraints - 1];
	}
	lin_prog->num_constraints--;
	constraint_free(lin_prog->constraints + lin_prog->num_constraints);
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

	/* Read objective function */
	for (size_t var = 0; var < num_variables; ++var)
		if (fscanf(fp, (var==num_variables-1) ? "%lf\n" : "%lf ", lin_prog.objective+var) != 1)
			error(1, "Cannot read objective function.");

	/* Read constraint upper-bounds, create constraints */
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
	}

	/* Read constraint linear combinations */
	char *line = NULL;
	size_t len = 0;
	for (size_t c = 0; c < num_constraints; ++c) {
		if (getline(&line, &len, fp) == -1)
			error(1, "Incorrect number of constraints");
		char *p = line;
		for (size_t var = 0; var < num_variables; ++var)
			lin_prog.constraints[c].linear_combination[var] = strtod(p, &p);
	}
	free(line);

	return lin_prog;
}

void lp_prune(struct LP *lin_prog)
{
	struct Constraint *constraints = lin_prog->constraints;
	for (size_t c=0; c<lin_prog->num_constraints; ++c) {
		struct Constraint c_constr = constraints[c];

		for (size_t d=c+1; d<lin_prog->num_constraints; ++d) {
			if (constraint_implies(&c_constr, constraints+d, lin_prog->num_variables)) {
				lp_remove_constraint(lin_prog, d--);
			} else if (constraint_implies(&c_constr, constraints+d, lin_prog->num_variables)) {
				lp_remove_constraint(lin_prog, c--);
				break;
			}
		}
	}
}

void lp_normalise_variable(struct LP *lin_prog, size_t const variable)
{
	for (size_t c=0; c<lin_prog->num_constraints; ++c)
		constraint_normalise_variable(lin_prog->constraints + c, lin_prog->num_variables, variable);
}
