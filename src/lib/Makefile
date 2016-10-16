# Makefile for 'lib' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o *.a
OBJ=io/libio.o
HEADERS=

auxlib.o: make-io
	ld -r $(OBJ) -o auxlib.o

make-io:
	$(MAKE) -C io libio.o CXX=$(CXX)

clean:
	$(CLEAN)
	$(MAKE) -C io clean
