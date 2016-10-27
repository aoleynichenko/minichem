# Makefile for the 'scf' module

CC=mpicc
CCFLAGS=-Wall -O3 -fopenmp -Werror=return-type -Werror=implicit-function-declaration
CLEAN=rm -rf *.o *.a
OBJ=scf.o 1e.o 2e.o boys.o ints_print.o scf_input.o molden.o eht.o analysis.o diis.o rhf.o uhf.o

#all:
#	$(CC) $(CCFLAGS) scf.c -o scf -llapack -lscalapack -lblas

1e.o: 1e.c
	$(CC) $(CCFLAGS) -c 1e.c -o 1e.o

2e.o: 2e.c
	$(CC) $(CCFLAGS) -c 2e.c -o 2e.o
	
boys.o: boys.c
	$(CC) $(CCFLAGS) -c boys.c -o boys.o

ints_print.o: ints_print.c
	$(CC) $(CCFLAGS) -c ints_print.c -o ints_print.o

scf.o: scf.c
	$(CC) $(CCFLAGS) -c scf.c -o scf.o  -llapack -lscalapack -lblas -lcblas

scf_input.o: scf_input.c
	$(CC) $(CCFLAGS) -c scf_input.c -o scf_input.o
	
molden.o: molden.c
	$(CC) $(CCFLAGS) -c molden.c -o molden.o

eht.o: eht.c
	$(CC) $(CCFLAGS) -c eht.c -o eht.o
	
uhf.o: uhf.c
	$(CC) $(CCFLAGS) -c uhf.c -o uhf.o
	
rhf.o: rhf.c
	$(CC) $(CCFLAGS) -c rhf.c -o rhf.o
	
diis.o: diis.c
	$(CC) $(CCFLAGS) -c diis.c -o diis.o

analysis.o: analysis.c
	$(CC) $(CCFLAGS) -c analysis.c -o analysis.o -lcblas

libscf.a: $(OBJ)
	ar rcs libscf.a $(OBJ)
	
clean:
	$(CLEAN)

