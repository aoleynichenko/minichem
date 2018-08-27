# Minichem top-level makefile

CC=mpicc
FC=mpif90
CCFLAGS=-Wall
SRC=main.c
OBJS=util/util.o input/input.o scf/libscf.a main.o
TARGET=minichem.x
CLEAN=rm -rf *.o

.PHONY: all clean install uninstall

all: $(TARGET)
clean:
	$(CLEAN)
	$(MAKE) -C util clean
	$(MAKE) -C input clean
	$(MAKE) -C scf clean
install:
	mkdir ../bin ../lib
	cp $(TARGET) ../bin
uninstall:
	rm -r ../bin ../lib

$(TARGET): make-util make-input make-scf main.o
	$(CC) $(CCFLAGS) -lm -fopenmp -o $(TARGET) $(OBJS) -llapack -lblas -lcblas -llapacke
		
main.o: main.c
	$(CC) $(CCFLAGS) -c -I util -I input -o main.o main.c

make-util:
	$(MAKE) -C util util.o CC=$(CC)

make-input:
	$(MAKE) -C input input.o CC=$(CC)

make-scf:
	$(MAKE) -C scf libscf.a CC=$(CC)





