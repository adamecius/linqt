#hopping.cpp hopping.hpp  

FOU_LOC_SRCS := fourier/fourier_transform.cpp 
FOU_LOC_OBJS := $(LAT_LOC_SRCS:.cpp=.o) 
SRCS+=$(LAT_LOC_SRCS)
OBJS+=$(LAT_LOC_OBJS)


FOU_CC=icc 
FOU_OPT=-O3 -fpic  -ipo 
FOU_LDFLAGS= -L/opt/intel/compilers_and_libraries/linux/mkl/lib/intel64
FOU_CFLAGS= -I../include -I/opt/intel/compilers_and_libraries/linux/mkl/include 
FOU_LIBS=-limf -lirc -qopenmp -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -ldl -liomp5
#MKL_LIB="-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl "


all:	fourier_transform.o

%.o : %.cpp
	$(FOU_CC) $(FOU_OPT) $(FOU_CFLAGS) $(FOU_LDFLAGS) $(LIBS) -c $< -o $@ 
