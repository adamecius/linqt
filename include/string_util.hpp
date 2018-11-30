#ifndef PARSER_STRING_UTIL_HPP
#define PARSER_STRING_UTIL_HPP

#include <vector>
#include <string>
#include <sstream>

namespace qt
{

	template <typename T>
	std::string to_string(const T x)
	{
		std::string Result;		// string which will contain the result
		std::ostringstream convert;	// stream used for the conversion
		convert << x;			// insert the textual representation of 'Number' in the characters in the stream
		Result = convert.str();		// set 'Result' to the contents of the stream
		return Result;
	}

	template <typename T>
	bool str2val(const std::string str, T& val)
	{
		std::istringstream ss(str);
		
		while (!ss.eof())
		{
			if (ss >> val)
			return true;
		}
		return false; // There is no integer!
	}

	template <typename T>
	bool str2vector(const std::string line, std::vector<T>& vec)
	{

		//tokenizer
		std::string token("");
		std::vector<std::string> token_array;

		for(int i=0;i<line.size();i++)
		if( isdigit(char(line[i])) || line[i]=='.' )
			token+=line[i];
		else if( token.size()!=0)
		{
			token_array.push_back(token);
			token="";
		}	
		vec=std::vector<T>(token_array.size() );
		
		bool isSuccess=true;
		for( int tk=0;tk<token_array.size();tk++)
		{
			std::istringstream ss(token_array[tk]);
			while (!ss.eof())
				isSuccess*=(bool)(ss >> vec[tk]);
		}
		 
	return isSuccess;
	}


}







#endif
