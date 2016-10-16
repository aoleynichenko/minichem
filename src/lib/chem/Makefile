# Makefile for 'lib/chem' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=Molecule.o
HEADERS=

libchem.o: $(OBJ)
	ld -r $(OBJ) -o libchem.o

Molecule.o: Molecule.cc
	$(CXX) $(CXXFLAGS) -c -o Molecule.o Molecule.cc

clean:
	$(CLEAN)
