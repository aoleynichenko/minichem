# Minichem's top-level makefile
# C++ version

CC=mpicc
CXX=g++
FC=mpif90
CCFLAGS=-Wall
CXXFLAGS=-Wall
SRC=main.c
OBJS=main.o OutputStream.o
TARGET=minichem
CLEAN=rm -rf *.o

.PHONY: all clean install uninstall

all: $(TARGET)
clean:
	$(CLEAN)
install:
	mkdir ../bin ../lib
	cp $(TARGET) ../bin
uninstall:
	rm -r ../bin ../lib

$(TARGET): main.o OutputStream.o
	$(CXX) $(CXXFLAGS) -lm -fopenmp -o $(TARGET) $(OBJS) -llapack -lblas -lcblas -llapacke
		
main.o: main.cc
	$(CXX) $(CXXFLAGS) -c -o main.o main.cc

OutputStream.o: OutputStream.cc
	$(CXX) $(CXXFLAGS) -c -o OutputStream.o OutputStream.cc


