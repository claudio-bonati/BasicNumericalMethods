CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries= example_integration \
	example_linear_equations \
	example_ode_initial \
	example_ode_tpbvp \
	example_ode_eigen \
	example_root_find \
        example_spectrum \
	test_double \
	test_float 

all: $(binaries)

example_integration: ../src/example_integration.c integration.o
	$(CC) $(CFLAGS) -c ../src/example_integration.c 
	$(CC) $(CFLAGS) integration.o example_integration.o -o $@ -lm 

example_linear_equations: ../src/example_linear_equations.c linear_algebra.o
	$(CC) $(CFLAGS) -c ../src/example_linear_equations.c 
	$(CC) $(CFLAGS) linear_algebra.o example_linear_equations.o -o $@ -lm 

example_ode_eigen: ../src/example_ode_eigen.c linear_algebra.o ode_init.o ode_eigen.o root_find.o 
	$(CC) $(CFLAGS) -c ../src/example_ode_eigen.c 
	$(CC) $(CFLAGS) linear_algebra.o ode_init.o ode_eigen.o example_ode_eigen.o root_find.o -o $@ -lm 

example_ode_initial: ../src/example_ode_initial.c ode_init.o
	$(CC) $(CFLAGS) -c ../src/example_ode_initial.c 
	$(CC) $(CFLAGS) ode_init.o example_ode_initial.o -o $@ -lm 

example_ode_tpbvp: ../src/example_ode_tpbvp.c linear_algebra.o ode_init.o ode_tpbvp.o root_find.o 
	$(CC) $(CFLAGS) -c ../src/example_ode_tpbvp.c 
	$(CC) $(CFLAGS) linear_algebra.o ode_init.o ode_tpbvp.o example_ode_tpbvp.o root_find.o -o $@ -lm 

example_root_find: ../src/example_root_find.c root_find.o
	$(CC) $(CFLAGS) -c ../src/example_root_find.c 
	$(CC) $(CFLAGS) root_find.o example_root_find.o -o $@ -lm 

example_spectrum: ../src/example_spectrum.c linear_algebra.o
	$(CC) $(CFLAGS) -c ../src/example_spectrum.c 
	$(CC) $(CFLAGS) linear_algebra.o example_spectrum.o -o $@ -lm 

integration.o: ../lib/integration.c ../include/integration.h
	$(CC) $(CFLAGS) -c ../lib/integration.c 

linear_algebra.o: ../lib/linear_algebra.c ../include/linear_algebra.h
	$(CC) $(CFLAGS) -c ../lib/linear_algebra.c 

ode_eigen.o: ../lib/ode_eigen.c ../include/ode_eigen.h
	$(CC) $(CFLAGS) -c ../lib/ode_eigen.c

ode_init.o: ../lib/ode_init.c ../include/ode_init.h
	$(CC) $(CFLAGS) -c ../lib/ode_init.c

ode_tpbvp.o: ../lib/ode_tpbvp.c ../include/ode_tpbvp.h
	$(CC) $(CFLAGS) -c ../lib/ode_tpbvp.c

root_find.o: ../lib/root_find.c ../include/root_find.h
	$(CC) $(CFLAGS) -c ../lib/root_find.c 

test_double: ../src/test_double.c
	$(CC) $(CFLAGS) ../src/test_double.c -o $@ -lm

test_float: ../src/test_float.c
	$(CC) $(CFLAGS) ../src/test_float.c -o $@ -lm

.PHONY: clean cleanobj
clean:
	rm -f $(binaries) *.o *.dat

cleanobj:
	rm -f *.o
