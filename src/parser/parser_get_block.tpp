

//This function tries to get the values from a block, if it fails it return false
template <typename T>
bool 
parser::try_GetBlock(std::ifstream& infile, const std::string tag,  std::vector<T>& val)
{
	std::string begBlock="BEGIN_"+tag;
	std::string endBlock="END_"+tag;
	std::vector<std::string> toklist;
	bool found_beg=false, found_end=false;
	infile.clear();
	std::string	line; //variable to store each line of the file
	if ( infile ) // Check if file is valid
	try
	{
		//Check if the beginning and end of the block are present
		while	(  getline(infile , line ) )
		{	
			//Remove comments and spaces
			line = TrimLine( line );
			found_beg +=  line.find(begBlock, 0) 	!= std::string::npos ;
			found_end +=  line.find(endBlock, 0) 	!= std::string::npos ;
			if(  found_beg && found_end )
				break;
		}
		infile.clear(); 
		infile.seekg(0,std::ios::beg);
		found_beg=false; found_end=false;
		while	(  getline(infile , line ) )
		{	
			//Remove comments and spaces
			std::string header = TrimLine( line );
			//Get the tokens if is not beginning or end lnie
			found_end +=  header.find(endBlock, 0) 	!= std::string::npos ;
			if( found_beg && !found_end )
				Tokenize( line, " ",  toklist );			
			found_beg +=  header.find(begBlock, 0) 	!= std::string::npos ;
			//Break before reaching eof
			if(  found_beg && found_end )
				break;
		}
		val = std::vector<T>( toklist.size()) ;
		for(int i=0; i<toklist.size(); i++)
			qt::str2val(toklist[i],val[i]);
		return found_beg && found_end;
	}
	

	catch(std::ifstream::failure e)
	{
		std::cerr<<"GetBlock function could find the beginning or end of Block: "<<tag<<std::endl;
		infile.clear();
		return false;
	}
	return false;
};


template <typename T>
bool 
parser::GetBlock(std::ifstream& infile, const std::string text, std::vector<T>& val, bool  optional=false  )
{
	infile.clear(); 
	infile.seekg(0,std::ios::beg);

	bool wasfound=try_GetBlock(infile,text,val);
	if( !wasfound )
	{
		wasfound=true;
		if( !optional)
			std::cerr<<"Field: "<<text<<" not found in .CFG file"<<" aborting "<<std::endl;
		
	}

	infile.clear(); 
	infile.seekg(0,std::ios::beg);
	
	return wasfound;
};
