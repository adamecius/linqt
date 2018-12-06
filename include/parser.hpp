#ifndef PARSER_HPP
#define PARSER_HPP

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include "string_util.hpp" //qt:: conversion functions
#include "parser_util.hpp"

namespace parser
{

//This function tries to get the requested value, if it fails it return false
template <typename T>
bool 
try_GetValue(std::ifstream& infile, const std::string tag, T& val);

template <typename T>
bool 
GetValue(std::ifstream& infile, const std::string text, T& val, bool  optional=false  , bool  fromInit=true);

template <typename T>
bool 
GetVector(std::ifstream& infile, const std::string text, std::vector<T>& val,int dim, bool  optional=false  );


//This function tries to get the values from a block, if it fails it return false
template <typename T> 
bool 
GetBlock(std::ifstream& infile, const std::string text, std::vector< std::vector<T> >& block,const int dim0,const int  dim1, bool  optional=false  );


template <typename T>
std::ifstream& 
GetBoolOption(std::ifstream& infile, const T text, bool& value);

}


//template implementations
#include "../src/parser/parser_get_value.tpp"
#include "../src/parser/parser_get_block.tpp"
#include "../src/parser/parser_get_vector.tpp"
#include "../src/parser/parser_get_boolean.tpp"

#endif
