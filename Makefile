# Makefile for TB-QuantumTransport
# Copyright (C) 2016 Jose Hugo Garcia Aguilar <adamecius@gmail.com>
#

#--------------------------------------------------------------
# This makefile compiles the project TB-QuantumTransp which is 
# projected to use a combination of C, Cx11, intel mkl, nvidia cuda,
# mpi and openmp libraries, therefore the compilation and linking
# is very herogeneus. Due to that, I have decided to separte into
# three different compilations, which later may become three different 
# makefiles. 
#--------------------------------------------------------------


#---------------SPECIAL DIRECTORIES----------------------------
ICCDIR =$(INTEL_HOME)
MKLDIR =$(ICCDIR)/mkl
MPIDIR= $(I_MPI_ROOT)/intel64
INTEL_SEQ =$(MKLDIR)
INTEL_PAR =$(MKLDIR)


#There are different flavors
#for the compilation
FLAVOR=mpi+openmp

#SPECIAL COMPILATION DIRECTIVES
MODULES=  kpm kpm/parallel lattice fourier random

#look for include files in
#each of the modules
CFLAGS+= -I$(MKLDIR)/include $(patsubst %,-I%,$(MODULES)) -Iinclude  
LDFLAGS+= -L$(MKLDIR)/lib/intel64 $(patsubst %,-L%,$(MODULES)) 
CPPFLAGS=
LIBS=
#include the description for
#each module
include $(patsubst %,%/module.mk,$(MODULES))
#---------------COMPILERS AND FLAGS----------------------------

# define the C++ compiler to use
#if the MPI COMPILER is defined it takes precedence over the sequential
CC=icpc # define the C compile-time flags
ifneq ($(MPICC),"")
CC=$(MPICC)
endif
 
CFLAGS+=-DBENCHMARK -O3 -ipo -fpic   

#---------------LOCAL SOURCES AND TARGETS-------------------------
LOC_SRCS:=$(wildcard src/*.cpp) #$(addprefix src/, lattice_index.cpp  TB-QuantumTransp.cpp  regular_hamiltonian.cpp kpm.cpp  lattice.cpp onsite_disorder.cpp irregular_hamiltonian.cpp)  # should be .cpp
LOC_OBJS:=$(addprefix obj/, $(notdir $(LOC_SRCS:.cpp=.o)) ) 


SRCS+=$(LOC_SRCS)
OBJS+=$(LOC_OBJS)
	
# define the executable file 
MAIN = conductivity_moments.mpi+omp

#---------------THE REST OF THE MAKE SHOULD BE GENERIC------------------

.PHONY: depend clean

all:    $(MAIN)
		@echo $(LIBS)
#		@echo  the program $(MAIN) was compiled 

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $(MAIN)

obj/%.o : src/%.cpp
	$(CC) $(CFLAGS) -c $<  -o $@

clean:
	$(RM) obj/*.o*  *~ $(MAIN) $(patsubst %,%/*.o*,$(MODULES))

depend: $(SRCS)
		makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
