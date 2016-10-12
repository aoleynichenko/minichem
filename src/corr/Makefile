# Makefile for the 'corr' module

CC=gcc
CCFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=ao2mo.o

ao2mo.o: ao2mo.c
	$(CC) $(CCFLAGS) -c ao2mo.c -o ao2mo.o
	
corr.o: $(OBJ)
	ld -r $(OBJ) -o corr.o

clean:
	$(CLEAN)
