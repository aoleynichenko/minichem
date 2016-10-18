# Makefile for 'lib/io' module

CC=gcc
CXX=g++
FC=gfortran
CCFLAGS=-Wall
CXXFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=OutputStream.o Lexer.o Token.o Log.o
HEADERS=

libio.o: $(OBJ)
	ld -r $(OBJ) -o libio.o

Lexer.o: Lexer.cc
	$(CXX) $(CXXFLAGS) -c -o Lexer.o Lexer.cc

Log.o: Log.cc
	$(CXX) $(CXXFLAGS) -c -o Log.o Log.cc

OutputStream.o: OutputStream.cc
	$(CXX) $(CXXFLAGS) -c -o OutputStream.o OutputStream.cc

Token.o: Token.cc
	$(CXX) $(CXXFLAGS) -c -o Token.o Token.cc

clean:
	$(CLEAN)
