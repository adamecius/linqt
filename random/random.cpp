



#include "random.hpp"

unsigned int NumCal::random::random_seed()
{
	unsigned int pid=getpid(); // get it as per your OS
	unsigned int seed= time(NULL)*pid;
	return seed;
}
