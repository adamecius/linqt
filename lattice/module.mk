#hopping.cpp hopping.hpp  

LAT_LOC_SRCS := lattice/onsite_disorder.cpp \
		lattice/lattice.cpp \
		lattice/irregular_hamiltonian.cpp  \
	
LAT_LOC_OBJS := $(LAT_LOC_SRCS:.cpp=.o) 
SRCS+=$(LAT_LOC_SRCS)
OBJS+=$(LAT_LOC_OBJS)

lattice/%.o : lattice/%.cpp
	$(CC) $(CFLAGS) -c $<  -o $@

