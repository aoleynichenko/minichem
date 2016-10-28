# Makefile for 'lib' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o *.a
OBJ=io/libio.o chem/libchem.o basis/libbasis.o wf/libwf.o \
	qscript/libqscript.o
HEADERS=

auxlib.o: make-io make-chem make-basis make-wf make-qscript
	ld -r $(OBJ) -o auxlib.o

make-io:
	$(MAKE) -C io libio.o CXX=$(CXX)

make-chem:
	$(MAKE) -C chem libchem.o CXX=$(CXX)

make-basis:
	$(MAKE) -C basis libbasis.o CXX=$(CXX)

make-wf:
	$(MAKE) -C wf libwf.o CXX=$(CXX)

make-qscript:
	$(MAKE) -C qscript libqscript.o CXX=$(CXX)

clean:
	$(CLEAN)
	$(MAKE) -C io clean
	$(MAKE) -C chem clean
	$(MAKE) -C basis clean
	$(MAKE) -C wf clean
	$(MAKE) -C qscript clean


