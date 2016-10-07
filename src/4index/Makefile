# Makefile for the '4index' module

CC=mpicc
CCFLAGS=-Wall -O3 -fopenmp -Werror=return-type -Werror=implicit-function-declaration
CLEAN=rm -rf *.o *.a
OBJ=

#all:
#	$(CC) $(CCFLAGS) scf.c -o scf -llapack -lscalapack -lblas


#analysis.o: analysis.c
#	$(CC) $(CCFLAGS) -c analysis.c -o analysis.o -lcblas

lib4index.a: $(OBJ)
	ar rcs lib4index.a $(OBJ)
	
clean:
	$(CLEAN)

