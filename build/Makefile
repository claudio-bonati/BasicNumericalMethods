CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries=test_float test_double example_root_find example_integration

all: $(binaries)

test_float: ../src/test_float.c
	$(CC) $(CFLAGS) ../src/test_float.c -o $@ -lm

test_double: ../src/test_double.c
	$(CC) $(CFLAGS) ../src/test_double.c -o $@ -lm

root_find.o: ../lib/root_find.c ../include/root_find.h
	$(CC) $(CFLAGS) -c ../lib/root_find.c 

example_root_find: ../src/example_root_find.c root_find.o
	$(CC) $(CFLAGS) -c ../src/example_root_find.c 
	$(CC) $(CFLAGS) root_find.o example_root_find.o -o $@ -lm 

integration.o: ../lib/integration.c ../include/integration.h
	$(CC) $(CFLAGS) -c ../lib/integration.c 

example_integration: ../src/example_integration.c integration.o
	$(CC) $(CFLAGS) -c ../src/example_integration.c 
	$(CC) $(CFLAGS) integration.o example_integration.o -o $@ -lm 


#pi_no_rng: ../src/pi_no_rng.c
#	$(CC) $(CFLAGS) ../src/pi_no_rng.c -o $@ -lm
#
#random.o: ../lib/random.c ../include/pcg32min.h ../include/random.h
#	$(CC) $(CFLAGS) -c ../lib/random.c
#
#pcg32min.o: ../lib/pcg32min.c ../include/pcg32min.h 
#	$(CC) $(CFLAGS) -c ../lib/pcg32min.c
#
#pi_ss: ../src/pi_ss.c random.o pcg32min.o
#	$(CC) $(CFLAGS) -c ../src/pi_ss.c  
#	$(CC) $(CFLAGS) pcg32min.o random.o pi_ss.o -o $@ -lm
#

.PHONY: clean
clean:
	rm -f $(binaries) *.o

cleanobj:
	rm -f *.o
