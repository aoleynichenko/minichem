# Makefile for 'lib' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o *.a
OBJ=io/libio.o chem/libchem.o
HEADERS=

auxlib.o: make-io make-chem
	ld -r $(OBJ) -o auxlib.o

make-io:
	$(MAKE) -C io libio.o CXX=$(CXX)

make-chem:
	$(MAKE) -C chem libchem.o CXX=$(CXX)

clean:
	$(CLEAN)
	$(MAKE) -C io clean
	$(MAKE) -C chem clean
