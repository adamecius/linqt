#ifndef TIME_STATION_HPP
#define TIME_STATION_HPP


#include<chrono>
#include<iostream>


typedef double  r_value_t;

class time_station{
        public:
        time_station(){  start_ = std::chrono::steady_clock::now();};

        time_station(int time, std::string msg){
          print_msg(time, msg);
        };//for whenever all has been added previously

  
        inline void stop(){
          end_ = std::chrono::steady_clock::now();
          time_microSec_ =  std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count() ;
        }


       inline void stop(std::string msg ){    
         end_ = std::chrono::steady_clock::now();
         time_microSec_ =  std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count() ;

         print_msg(time_microSec_, msg);
       };

  
       inline void stop_add(int* time, std::string msg ){

         end_ = std::chrono::steady_clock::now();
         time_microSec_ =  std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count() ;

         *time += time_microSec_;

         print_msg(time_microSec_, msg);
       };
  
       inline void stop_add(int* time ){
         end_ = std::chrono::steady_clock::now();
         time_microSec_ =  std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count() ;

         *time += time_microSec_;
       };
  
  
       inline void print_msg(int time, const std::string& msg ){
    
         int sec   =  time/1000000;
         int min   = sec / 60;
         int hour  = min / 60;
         int reSec = sec % 60;
         int reMin = min % 60;
  
         std::cout<<msg;

         if(hour>0)
           std::cout<<hour<<" hrs, "<<reMin<<" mins;"<<std::endl;
         else if(min>0)
           std::cout<<min<<" mins, "<<reSec<<" secs;"<<std::endl;
         else if(sec>0)
           std::cout<<sec<<" secs;"<<" ("<< r_value_t(time)/r_value_t(1000)<<"ms) "<<std::endl;
         else
          std::cout<< r_value_t(time)/r_value_t(1000)<<" ms "<<std::endl;    
       }

  
       private:
       std::chrono::steady_clock::time_point start_, end_;
       int time_microSec_;
  
};

#endif //TIME_STATION_HPP
