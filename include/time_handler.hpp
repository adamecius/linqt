#ifndef TIME_HANDLER
#define TIME_HANDLER

#include <string> 
#include <libconfig.h++>

namespace timeConst{ const int INF = -1; };

class time_unit
{
	
	private:
	template <typename T> 
	std::string to_string( T x ) 
	{
		std::stringstream ss;
		ss << x;
		return ss.str();
	}
	public:
	time_unit(): real_value(0),  string_value(""){}
	
	time_unit(const double x) 			// conversion from constructor:
	{
		setTime( x );
	}			
	void operator= (const double  x)	// conversion from  (assignment):
	{
		setTime( x );
	} 	

	bool Infinite( )
	{
		return string_value.compare( "INF" ) == 0;
	}

	void setTime( const double x )
	{
		if ( string_value.compare( "INF" ) != 0)
		{  
			real_value = x;	
			if ( x < 0 )
			{
				real_value = timeConst::INF;
				string_value="INF" ;
			}
			else
				string_value=to_string(x);
		}
	}


    operator double() const { return real_value; }
    operator std::string() const { return string_value; }


	void ReadTime(	libconfig::Config& cfg, const std::string timeId)
	{
		try
		{
			cfg.lookupValue(timeId,string_value);
		}
		catch(const libconfig::ParseException &)
		{
			std::cerr<<" ERRRO WHILE READING "<<timeId<<std::endl;
		}
		if ( string_value.compare( "INF")==0 )
			real_value   = timeConst::INF;
		else
		{
			try
			{
				cfg.lookupValue(timeId,real_value);
				string_value = to_string(real_value);
			}
			catch(const libconfig::ParseException &)
			{
				std::cerr<<" ERRRO WHILE READING "<<timeId<<std::endl;
			}
		}
		if( real_value< 0 )
		{
			real_value= timeConst::INF;
			string_value = "INF";
		}
	}
	
	private:
	double real_value;
	std::string string_value;
};


#endif
