# Makefile for 'lib' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o *.a
OBJ=io/libio.o chem/libchem.o basis/libbasis.o
HEADERS=

auxlib.o: make-io make-chem make-basis
	ld -r $(OBJ) -o auxlib.o

make-io:
	$(MAKE) -C io libio.o CXX=$(CXX)

make-chem:
	$(MAKE) -C chem libchem.o CXX=$(CXX)

make-basis:
	$(MAKE) -C basis libbasis.o CXX=$(CXX)

clean:
	$(CLEAN)
	$(MAKE) -C io clean
	$(MAKE) -C chem clean
	$(MAKE) -C basis clean
