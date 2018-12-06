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

		int i=0;
		while( i <line.size())
		{
			for(i;i < line.size() ; i++) 
			if( std::isdigit(char(line[i])) || line[i]=='.' || line[i]=='+' || line[i]=='-' )
				token+=line[i];
			else 
				break;
			token_array.push_back(token);
			token="";		
			i++;
		}
		
		//Array for the tokents
		const int tkDim = token_array.size();
		if( tkDim == 0 ) return  false;

		vec=std::vector<T>(tkDim);
		for( int tk=0;tk<tkDim;tk++)
		{
			std::istringstream ss(token_array[tk]);
			ss >> vec[tk];
		}
	return true;
	}
}







#endif
