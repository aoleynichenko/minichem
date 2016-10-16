# Makefile for 'lib/io' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=OutputStream.o
HEADERS=

libio.o: $(OBJ)
	ld -r $(OBJ) -o libio.o

OutputStream.o: OutputStream.cc
	$(CXX) $(CXXFLAGS) -c -o OutputStream.o OutputStream.cc

clean:
	$(CLEAN)
