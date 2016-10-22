# Makefile for 'lib/basis' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall -std=c++11
CLEAN=rm -rf *.o
OBJ=BasisSet.o
HEADERS=

libbasis.o: $(OBJ)
	ld -r $(OBJ) -o libbasis.o

BasisSet.o: BasisSet.cc
	$(CXX) $(CXXFLAGS) -c -o BasisSet.o BasisSet.cc

clean:
	$(CLEAN)
