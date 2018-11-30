#ifndef PARSER_UTIL_HPP
#define PARSER_UTIL_HPP

#include <string>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <algorithm>
#include "string_util.hpp" //qt:: conversion functions

namespace parser
{

const bool fromInit=true;


const bool CurrLine=false;

const bool Optional=true;

const bool NotOptional=false;



size_t 
FindHeader(	std::ifstream& 
			infile, const std::string text, 
			bool fromInitLine=true);



std::string 
RemoveComments( std::string inLine);	

bool 
fileExists (const std::string& name);


std::string 
TrimLine( std::string inLine);


void 
SeparateFieldValue( std::string line, std::string& field, std::string& value );


std::ifstream& 
GetBoolOption(std::ifstream& infile, const std::string text, bool& value);


void 
Tokenize( std::string line, std::string sep,  std::vector<std::string>& toklist );

};
#endif
