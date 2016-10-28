# Makefile for 'lib/chem' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall -std=c++11
CLEAN=rm -rf *.o
OBJ=Molecule.o Elements.o
HEADERS=

libchem.o: $(OBJ)
	ld -r $(OBJ) -o libchem.o

Molecule.o: Molecule.cc
	$(CXX) $(CXXFLAGS) -c -o Molecule.o Molecule.cc

Elements.o: Elements.cc
	$(CXX) $(CXXFLAGS) -c -o Elements.o Elements.cc

clean:
	$(CLEAN)
