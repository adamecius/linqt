#ifndef TIMING_HPP
#define TIMING_HPP

#include <ctime>
#include <cstdio>

#ifndef BENCHMARK
static inline void
cycletime (int ncycle)
{
}
;
//dummy function
#endif

#ifdef BENCHMARK

static inline void cycletime(int _ncycle)
{

  const int ncall_total= _ncycle; //expected total number of calls
  static int ncalls;// current number of calls
  static time_t start,end;//time variables 											//Static variables that will save the time
  static bool first_call=true;//bolean variables which tells you if this is a begin or a end call													//
  static double total_time;
  static double print_time;	//ammount of time for each print
  if(_ncycle==-1)//When casted with -1, reset the variables
    {
      first_call=true;
      ncalls=0;
      total_time=0;
      print_time=0;
    }
  else
    {
      if(first_call)
	{
	  start=clock();
	  first_call=false;
	}
      else
	{
	  end=clock();
	  first_call=true;	//this will catch the time at the end of cycle
	  const double cycle_time=(double)((double)(end-start)/(double)CLOCKS_PER_SEC)/KPM_OMP_NUM_THREADS;
	  print_time=print_time+cycle_time;
	  total_time =total_time+cycle_time;//this will calculate the diference  in seconds
	  ncalls=ncalls+1;
	  //this will store the  number of cycles
	  //std::cout<<"time enlapsep="<<total<<" ncalls="<<ncalls<<" "<<ncycle<<std::endl;
	  const double tcycle=total_time/(double)ncalls;//time per cycle
	  const double rest =tcycle*((double)(ncall_total-ncalls));//time until endç
	  if(print_time >1 )
	    {
	      if(rest <60 )
		{							//seconds
		  std::cout	<<"The cycle mean time is  "
		      <<1000*tcycle<<" ms. The whole set wil end in "
		      <<rest<<" sec, "
		      <<(ncalls*100.0f/ncall_total)
		      <<" completed"<<std::endl;
		}
	      else if(rest <3600)
		{							//minutes
		  std::cout	<<"The cycle mean time is  "
		      <<1000*tcycle<<" ms. The whole set wil end in "
		      <<rest/60<<" min, "
		      <<(ncalls*100.0f/ncall_total)
		      <<" completed"<<std::endl;
		}
	      else if( rest <86400)
		{								//hours
		  std::cout	<<"The cycle mean time is  "
		      <<1000*tcycle<<" ms. The whole set wil end in "
		      <<rest/3600<<" hours, "
		      <<(ncalls*100.0f/ncall_total)
		      <<" completed"<<std::endl;
		}
	      else
		{							//days
		  std::cout	<<"The cycle mean time is  "
		      <<1000*tcycle<<" ms. The whole set wil end in "
		      <<rest/86400<<" days, "
		      <<(ncalls*100.0f/ncall_total)
		      <<" completed"<<std::endl;
		}
	      print_time=0;
	    }

	}
    }

};
#endif

#endif
