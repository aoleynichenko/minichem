# Makefile for 'scf' module

CC=gcc
#CXX=g++
FC=gfortran
CCFLAGS=-Wall
#CXXFLAGS=-Wall
CLEAN=rm -rf *.o *.a
OBJ=rhf.o
HEADERS=

libscf.o: $(OBJ)
	ld -r $(OBJ) -o libscf.o

rhf.o: rhf.cc
	$(CXX) $(CXXFLAGS) -c -o rhf.o rhf.cc

clean:
	$(CLEAN)
