# Makefile for TB-QuantumTransport
# Copyright (C) 2016 Jose Hugo Garcia Aguilar <adamecius@gmail.com>
#

#--------------------------------------------------------------
# This makefile compiles the example of the KPM module.
# which uses a combination of C, Cx11, intel mkl, nvidia cuda,
# mpi and openmp libraries, therefore the compilation and linking
# is very heterogeneus. Due to that, I have decided to separte into
# three different compilations, which later may become three different 
# makefiles. 
#--------------------------------------------------------------

CC=icpc
CFLAGS = -I../include 
LDFLAGS= 
SRCS   = chebsyhev_set.cpp kpm_linalg.cpp kpm.cpp example-benchmark.cpp
OBJS   = $(SRCS:.cpp=.o) 
MAIN   = example-benchmark

CFLAGS+= -qopenmp #-I/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/include
#LDFLAGS+=-L/opt/intel/compilers_and_libraries_2016.3.210/linux/mkl/lib/intel64/
#LIBS+=-lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl

#---------------THE REST OF THE MAKE SHOULD BE GENERIC------------------

all:    $(MAIN)
		@echo  the program $(OBJS) was compiled 

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o $(MAIN)

%.o : %.cpp
	$(CC) $(CFLAGS) -c $<  -o $@

clean:
	$(RM) *.o*  *~ $(MAIN)


# DO NOT DELETE THIS LINE -- make depend needs it
