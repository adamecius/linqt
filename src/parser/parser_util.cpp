#include "parser_util.hpp" //qt:: conversion functions



size_t 
parser::FindHeader(	std::ifstream& infile, 
					const std::string text, 
					bool fromInitLine)					
{
	size_t line_pos=std::string::npos;
	
	if( fromInitLine )
	{
		infile.clear();
		infile.seekg (0, std::ios::beg);//Look from the start
	}
	
	std::string  line;
	if ( infile ) // Check if file is valid
	try
	{
		while	(  getline(infile , line ) )
		{
			//Remove everything after a comment if present
			size_t  com_pos;
			if( (com_pos=line.find('#', 0) )==std::string::npos ) com_pos=line.size();
			line=line.substr (0,com_pos);

			if(line.find(text, 0) != std::string::npos)
			{
					line_pos =(int) infile.tellg() - line.size() -1;
					return line_pos;
			}
		}
	}
	catch(std::ifstream::failure e)
	{
		if( infile.bad() == 1)
			std::cerr<<"FindHeader throw status iostream::bad. ABORTING"<<std::endl;
		infile.clear();
		return std::string::npos;
	}
	else
		std::cerr<<"FindHeader tried to read an invalid file. ABORTING"<<std::endl;

	return std::string::npos;
};




bool
parser::fileExists (const std::string& name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
};


std::string
parser::RemoveComments( std::string inLine)
{
	//Remove coments
	size_t pos = inLine.find("#");
	std::string outLine = inLine.substr (0,pos-1);
	return outLine;
};


std::string
parser::TrimLine( std::string inLine)
{
	size_t pos;
	std::string outLine;
	std::string::iterator end_pos;
	//Remove coments
	pos = inLine.find("#");
	outLine = inLine.substr (0,pos-1);
	//Remove spaces and tabulators
	end_pos = remove(outLine.begin(), outLine.end(),' ');
	outLine.erase(end_pos, outLine.end());
	end_pos = remove(outLine.begin(), outLine.end(),'\t');
	outLine.erase(end_pos, outLine.end());
	return outLine;
};




void
parser::SeparateFieldValue( std::string line, std::string& field, std::string& value )
{
		//Find where is the equal sign
		size_t 	equ_pos;
		if( (equ_pos=line.find('=', 0) )==std::string::npos ) equ_pos=line.size();
		//Assign what is before the equal to field, and after to value
		field=line.substr (0,equ_pos);
		value=line.substr (equ_pos+1,-1);
};


void
parser::Tokenize( std::string line, std::string sep,  std::vector<std::string>& toklist )
{
		replace(line.begin(), line.end(), '\t', ' ');

		//Find where is the equal sign
		std::string token;
		size_t 	sep_pos;
		while( true )
		{
			sep_pos=line.find(sep, 0);
			token=line.substr (0,sep_pos);
			line =line.substr (sep_pos+1,-1);
			TrimLine(token);
			if(  token.size() != 0)
				toklist.push_back(token);

			if( sep_pos==std::string::npos )
				break;
		}
};
