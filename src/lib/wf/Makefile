# Makefile for 'lib/wf' module

#CC=gcc
#CXX=g++
#FC=gfortran
#CCFLAGS=-Wall
#CXXFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=Wavefunction.o
HEADERS=

libwf.o: $(OBJ)
	ld -r $(OBJ) -o libwf.o

Wavefunction.o: Wavefunction.cc
	$(CXX) $(CXXFLAGS) -c -o Wavefunction.o Wavefunction.cc

RhfWavefunction.o: RhfWavefunction.cc
	$(CXX) $(CXXFLAGS) -c -o RhfWavefunction.o RhfWavefunction.cc

clean:
	$(CLEAN)
