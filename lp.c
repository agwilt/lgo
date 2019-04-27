#define _POSIX_C_SOURCE 201710L

#include "lp.h"
#include "misc.h"

#include <stdlib.h>

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
	/* Read number of constraints and variable */
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
		error(1, "It looks like create_lp_empty created an LP with not enough space for constraints");
	for (size_t c = 0; c < num_constraints; ++c) {
		lin_prog.constraints[c].linear_combination = malloc(num_variables * sizeof(double));
		if (fscanf(fp,
		           (c==num_constraints-1) ? "%lf\n" : "%lf ",
		           &(lin_prog.constraints[c].value)) != 1) {
			error(1, "Cannot read constraint upper-bounds.");
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

void lp_print(struct LP const* lin_prog)
{
	/* Print number of constraints and variables */
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

void lp_print_nice(struct LP const* lin_prog)
{
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

void lp_remove_constraint(struct LP *lin_prog, size_t const index)
{
	constraint_free(lin_prog->constraints + index);
	if (index < lin_prog->num_constraints - 1) {
		lin_prog->constraints[index] = lin_prog->constraints[lin_prog->num_constraints - 1];
	}
	lin_prog->num_constraints--;
	constraint_free(lin_prog->constraints + lin_prog->num_constraints);
}
