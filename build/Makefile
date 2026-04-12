CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries= example_integration \
	example_linear_algebra \
	example_root_find \
	test_double \
	test_float 

all: $(binaries)

example_integration: ../src/example_integration.c integration.o
	$(CC) $(CFLAGS) -c ../src/example_integration.c 
	$(CC) $(CFLAGS) integration.o example_integration.o -o $@ -lm 

example_linear_algebra: ../src/example_linear_algebra.c linear_algebra.o
	$(CC) $(CFLAGS) -c ../src/example_linear_algebra.c 
	$(CC) $(CFLAGS) linear_algebra.o example_linear_algebra.o -o $@ -lm 

example_root_find: ../src/example_root_find.c root_find.o
	$(CC) $(CFLAGS) -c ../src/example_root_find.c 
	$(CC) $(CFLAGS) root_find.o example_root_find.o -o $@ -lm 

integration.o: ../lib/integration.c ../include/integration.h
	$(CC) $(CFLAGS) -c ../lib/integration.c 

linear_algebra.o: ../lib/linear_algebra.c ../include/linear_algebra.h
	$(CC) $(CFLAGS) -c ../lib/linear_algebra.c 

root_find.o: ../lib/root_find.c ../include/root_find.h
	$(CC) $(CFLAGS) -c ../lib/root_find.c 

test_double: ../src/test_double.c
	$(CC) $(CFLAGS) ../src/test_double.c -o $@ -lm

test_float: ../src/test_float.c
	$(CC) $(CFLAGS) ../src/test_float.c -o $@ -lm

.PHONY: clean cleanobj
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
