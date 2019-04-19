CC = gcc
CFLAGS = -Wall -Wextra -Ofast -pedantic -std=c17

DEPS = Makefile
OBJ = obj

$(OBJ)/%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: fourier_motzkin

fourier_motzkin: fourier_motzkin.c $(OBJ)/misc.o $(OBJ)/constraint.o $(DEPS)
	$(CC) $(CFLAGS) fourier_motzkin.c $(OBJ)/misc.o $(OBJ)/constraint.o -o fourier_motzkin

.PHONY: clean bin_clean obj_clean
clean: bin_clean obj_clean
bin_clean:
	rm -fv fourier_motzkin
obj_clean:
	rm -rfv "$(OBJ)"
