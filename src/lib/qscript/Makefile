# Makefile for 'lib/qscript' module
# QuantumSqript essentials

#CC=gcc
#CXX=g++
#FC=gfortran
#CCFLAGS=-Wall
#CXXFLAGS=-Wall
CLEAN=rm -rf *.o
OBJ=QS_Object.o QS_Scope.o QS_RuntimeError.o
HEADERS=

libqscript.o: $(OBJ)
	@echo building QuantumScript essentials
	ld -r $(OBJ) -o libqscript.o

QS_Object.o: QS_Object.cc
	$(CXX) $(CXXFLAGS) -c -o QS_Object.o QS_Object.cc

QS_RuntimeError.o: QS_RuntimeError.cc
	$(CXX) $(CXXFLAGS) -c -o QS_RuntimeError.o QS_RuntimeError.cc

QS_Scope.o: QS_Scope.cc
	$(CXX) $(CXXFLAGS) -c -o QS_Scope.o QS_Scope.cc

clean:
	$(CLEAN)
