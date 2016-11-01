LIBINT_HOME=/usr/local/libint/2.2.0-beta1
LIBEIGEN_INCLUDE=/usr/include/eigen3

CXX=g++
CXXFLAGS=-std=c++11 -I $(LIBINT_HOME)/include -I $(LIBINT_HOME)/include/libint2 \
		 -I $(LIBEIGEN_INCLUDE)
LINKFLAGS=-L $(LIBINT_HOME)/lib -lint2
SRC=hartree-fock.c
OBJS=hartree-fock.o
TARGET=hartree-fock
CLEAN=rm -rf *.o
INCLUDES=

.PHONY: all clean

all: $(TARGET)
clean:
	$(CLEAN)
	rm -f hartree-fock

(TARGET): hartree-fock.o
	g++ $(CXXFLAGS) -lm -o $(TARGET) $(OBJS) $(LINKFLAGS)

hartree-fock.o: hartree-fock.cc
	$(CXX) $(CXXFLAGS) -c -o hartree-fock.o hartree-fock.cc

