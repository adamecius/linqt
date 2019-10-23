#hopping.cpp hopping.hpp  

FOU_LOC_SRCS := fourier/fourier_transform.cpp 
FOU_LOC_OBJS := $(FOU_LOC_SRCS:.cpp=.o) 
SRCS+=$(FOU_LOC_SRCS)
OBJS+=$(FOU_LOC_OBJS)
LIBS+= -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm -ldl
fourier/%.o : fourier/%.cpp
	$(CC) $(CFLAGS) -qopenmp -mt_mpi  -c $< -o $@ 
