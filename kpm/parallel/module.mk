
ifeq ($(FLAVOR),mpi+openmp)

MPICC=mpiicpc

CFLAGS+=-I$(MPIDIR)/include64  -qopenmp -mt_mpi
LDFLAGS+=-L$(MPIDIR)/lib64 

endif

