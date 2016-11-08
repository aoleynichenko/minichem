# Makefile for 'scf' module

CC=gcc
#CXX=g++
FC=gfortran
CCFLAGS=-Wall
#CXXFLAGS+=-Wall
CLEAN=rm -rf *.o *.a
OBJ=rhf.o uhf.o Diis.o
HEADERS=

.cc.o:
	$(CXX) $(CXXFLAGS) -c -o $@ $<

rhf.o: rhf.cc
Diis.o: Diis.cc

libscf.o: $(OBJ)
	ld -r $(OBJ) -o libscf.o

#rhf.o: rhf.cc
#	$(CXX) $(CXXFLAGS) -c -o rhf.o rhf.cc

#Diis.o: Diis.cc
#	$(CXX) $(CXXFLAGS) -c -o Diis.o Diis.cc

clean:
	$(CLEAN)
