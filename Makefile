CC = gcc
CFLAGS = -Wall -Wextra -g -pedantic -std=c17 -lm

DEPS = Makefile
OBJ = obj

$(OBJ)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: fourier_motzkin

fourier_motzkin: fourier_motzkin.c $(OBJ)/misc.o $(OBJ)/constraint.o $(OBJ)/lp.o $(DEPS)
	$(CC) $(CFLAGS) fourier_motzkin.c $(OBJ)/misc.o $(OBJ)/constraint.o $(OBJ)/lp.o -o fourier_motzkin

.PHONY: clean bin_clean obj_clean
clean: bin_clean obj_clean
bin_clean:
	rm -fv fourier_motzkin
obj_clean:
	rm -rfv "$(OBJ)/*"
