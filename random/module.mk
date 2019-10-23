#hopping.cpp hopping.hpp  

RAND_LOC_SRCS := random/random.cpp 
	
RAND_LOC_OBJS := $(RAND_LOC_SRCS:.cpp=.o) 
SRCS+=$(RAND_LOC_SRCS)
OBJS+=$(RAND_LOC_OBJS)

random/%.o : random/%.cpp
	$(CC) $(CFLAGS) -c $<  -o $@

