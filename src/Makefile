# Minichem's top-level makefile
# C++ version

CC=mpicc
CXX=g++
FC=mpif90
CCFLAGS=-Wall
CXXFLAGS=-Wall
SRC=main.c
OBJS=lib/auxlib.o main.o
TARGET=minichem
CLEAN=rm -rf *.o

.PHONY: all clean install uninstall

all: $(TARGET)
clean:
	$(CLEAN)
	$(MAKE) -C lib clean
install:
	mkdir ../bin ../lib
	cp $(TARGET) ../bin
uninstall:
	rm -r ../bin ../lib

$(TARGET): make-lib main.o
	$(CXX) $(CXXFLAGS) -lm -fopenmp -o $(TARGET) $(OBJS) -llapack -lblas -lcblas -llapacke
		
main.o: main.cc
	$(CXX) $(CXXFLAGS) -c -o main.o main.cc

make-lib:
	$(MAKE) -C lib auxlib.o CXX=$(CXX)

