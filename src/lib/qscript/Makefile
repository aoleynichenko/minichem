# Makefile for 'lib/qscript' module
# QuantumSqript essentials

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=QS_Object.o
HEADERS=

libqscript.o: $(OBJ)
	ld -r $(OBJ) -o libqscript.o

QS_Object.o: QS_Object.cc
	$(CXX) $(CXXFLAGS) -c -o QS_Object.o QS_Object.cc

clean:
	$(CLEAN)
