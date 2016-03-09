# Makefile for 'util' module

CC=gcc
FC=gfortran
CCFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=errors.o decor.o buildno.o memory.o
HEADERS=util.h

errors.o: errors.c
	$(CC) $(CCFLAGS) -c errors.c -o errors.o

memory.o: memory.c
	$(CC) $(CCFLAGS) -c memory.c -o memory.o

decor.o: decor.c
	$(CC) $(CCFLAGS) -c decor.c -o decor.o

buildno.o: buildno.c
	$(CC) $(CCFLAGS) -c buildno.c -o buildno.o

util.o: $(OBJ)
	ld -r $(OBJ) -o util.o

clean:
	$(CLEAN)
