#!/usr/bin/env python3

from sys import argv
import random

if len(argv) < 3:
    print("Usage:", argv[0], "num_constraints num_variables")
    exit(1)

num_constraints = int(argv[1])
num_variables = int(argv[2])
print(num_constraints, num_variables)
# Objective
print(" ".join([str(random.random()) for i in range(num_variables)]))
# Constraint upper-bounds
print(" ".join([str(random.randint(-20, 20)) for i in range(num_constraints)]))
# Constraint linear combinations
for c in range(num_constraints):
    print(" ".join([str(random.randint(-5,5)) for i in range(num_variables)]))

