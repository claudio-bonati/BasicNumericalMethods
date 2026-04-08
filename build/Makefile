CC=gcc
CFLAGS=-O3 -Wall -Wextra -Werror -pedantic -Wconversion -std=c99 

binaries=test_float test_double

all: $(binaries)

test_float: ../src/test_float.c
	$(CC) $(CFLAGS) ../src/test_float.c -o $@ -lm

test_double: ../src/test_double.c
	$(CC) $(CFLAGS) ../src/test_double.c -o $@ -lm

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
