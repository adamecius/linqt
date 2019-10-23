KPM_LOC_SRCS := kpm/kpm.cpp kpm/chebyshev_set.cpp kpm/kpm_utilities.cpp kpm/kpm_linalg.cpp
KPM_LOC_OBJS := $(KPM_LOC_SRCS:.cpp=.o) 
SRCS+=$(KPM_LOC_SRCS)
OBJS+=$(KPM_LOC_OBJS)

kpm/%.o : kpm/%.cpp
	$(CC) $(CFLAGS) -c $<  -o $@

