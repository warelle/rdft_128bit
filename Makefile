# Makefile

PROGRAM = main
SUB_OBJ = lu.o lib.o gen.o rdft.o matlib.o givens.o solve.o
OBJ     = main.o $(SUB_OBJ)
TESTER  = test.h
CC      = g++
CFLAGS  = -Wall -O2 -std=gnu++11
LIB     = -lquadmath -lm

LIBOBJ  = eigen.o
LIBCC   = g++
LIBFLG  = -Wall -O2 -std=gnu++11
INCL    = -I./eigen

.PHONY: all
all: main

$(PROGRAM): $(OBJ) $(TESTER) $(LIBOBJ)
	$(CC) $(CFLAGS) -o main $(OBJ) $(LIBOBJ) $(LIB)

$(LIBOBJ): eigen.cpp eigen.h $(TESTER)
	$(LIBCC) $(LIBFLG) $(INCL) -c $<

.c.o: $(TESTER)
	$(CC) $(CFLAGS) -c $<
.cpp.o: $(TESTER)
	$(CC) $(CFLAGS) -c $<

zero: test.cxx $(SUB_OBJ) $(LIBOBJ)
	$(CC) $(CFLAGS) -o test test.cxx $(SUB_OBJ) $(LIBOBJ) $(LIB)

.PHONY: clean
clean:
	rm $(PROGRAM) $(OBJ) $(LIBOBJ)
