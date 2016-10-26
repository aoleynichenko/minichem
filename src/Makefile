# minichem's top-level makefile
# C++ version of minichem code
#
# Alexander Oleynichenko, 2016
# For questions, alexvoleynichenko@gmail.com
#
# Libraries which are needed to be installed on your machine before
# build minichem:
#  - Eigen  >= 3.2.10
#  - Libint >= 2.2.0
# To properly build minichem, you should set LIBINT_HOME and LIBEIGEN_HEADERS
# path variables specific for your machine (examples are given).
# Attention: minichem has been tested only for x86_64 platform!


# To evaluate integrals over gaussian atomic basis function
# minichem uses libint2 library by Edward Valeev and coworkers.
# LIBINT_HOME contains header files and static libraries lib/libint2.a
# and lib/libint2.la. You can easily detect LIBINT_HOME directory on
# your computer, it consists of folders: include, lib, share

LIBINT_HOME=/usr/local/libint/2.2.0-beta1

# Directory with Eigen's headers and source code.
# On my laptop, there are following folders in $(LIBEIGEN_INCLUDE):
#    Eigen  signature_of_eigen3_matrix_library  unsupported

LIBEIGEN_INCLUDE=/usr/include/eigen3

CC=mpicc
CXX=g++
FC=mpif90
CCFLAGS=-Wall
CXXFLAGS=-Wall -std=c++11 -I $(LIBINT_HOME)/include -I $(LIBINT_HOME)/include/libint2 \
		 -I $(LIBEIGEN_INCLUDE)
LINKFLAGS=-L $(LIBINT_HOME)/lib -lint2
SRC=main.c
OBJS=lib/auxlib.o Kernel.o main.o
TARGET=minichem
CLEAN=rm -rf *.o
INCLUDES=

# common variables for all subdirectories
export CXX
export CXXFLAGS

.PHONY: all clean install uninstall

all: $(TARGET)
clean:
	$(CLEAN)
	rm -f minichem
	$(MAKE) -C lib clean
install:
	@mkdir ../bin ../lib
	@cp $(TARGET) ../bin
	@mkdir ../lib/basissets
	@cp basissets/* ../lib/basissets/
	@echo Congratulations! Now you have minichem on your computer!
	@echo Set MINICHEM_HOME variable to enable minichem\'s basis set library
	@echo Proper MINICHEM_HOME = `pwd`/..
uninstall:
	@rm -r ../bin ../lib

$(TARGET): make-lib Kernel.o main.o
	echo $(CXXFLAGS)
	$(CXX) $(CXXFLAGS) -lm -fopenmp -o $(TARGET) $(OBJS) $(LINKFLAGS)

main.o: main.cc
	$(CXX) $(CXXFLAGS) -c -o main.o main.cc

Kernel.o: Kernel.cc
	$(CXX) $(CXXFLAGS) -c -o Kernel.o Kernel.cc

make-lib:
	$(MAKE) -C lib auxlib.o -e CXX=$(CXX)
