
//This function tries to get the requested value, if it fails it return false
template <typename T>
bool parser::try_GetValue(std::ifstream& infile, const std::string tag, T& val)
{
	infile.clear();
	std::string text=tag;
	std::transform(text.begin(), text.end(), text.begin(), ::tolower);
	std::string	line; //variable to store each line of the file
	if ( infile ) // Check if file is valid
	try
	{
		while	(  getline(infile , line ) )
		{	
			/*
						std::string field, value;
			//Remove comments and spaces
			line = TrimLine( line );
			//Assign what is before the equal to field, and after to value			
			SeparateFieldValue(line,field,value);
			std::cout<<line<<" "<<" "<<field<<" "<<value<<std::endl;
			std::exit(-1);
			*/
			size_t 
			equ_pos,txt_pos, com_pos; //if not comment is found return the final line
			//tolower everything before equal
			if( (equ_pos=line.find('=', 0) )==std::string::npos ) equ_pos=line.size();
			std::transform(line.begin(), line.begin()+equ_pos, line.begin(), ::tolower);
			//Remove everything after a comment
			if( (com_pos=line.find('#', 0) )==std::string::npos ) com_pos=line.size();
			line=line.substr (0,com_pos); 
			//return the value
			if	(
					( txt_pos=line.find(text, 0) )	!= std::string::npos && 
					txt_pos<equ_pos
				)
			{
				std::string str 		= TrimLine(line.substr(equ_pos+1,line.size()));  
				return qt::str2val(str,val);
			}
		}

	}
	catch(std::ifstream::failure e)
	{
		if( infile.bad() == 1)
				std::cerr<<"GetValue throw status iostream::bad. ABORTING"<<std::endl;
		infile.clear();
		return false;
	}
	else
		std::cerr<<"GetValue tried to read an invalid file. ABORTING"<<std::endl;

	return false;
};


template <typename T>
bool parser::GetValue(std::ifstream& infile, const std::string text, T& val, bool  optional=false  , bool  fromInit=true)
{
	long curr_pos;
	
	if( fromInit)
		curr_pos=0;
	else
		curr_pos= infile.tellg(); //if optional=true and fromInit=False, the end position in the file is the initial one

	infile.clear(); 
	infile.seekg(curr_pos,std::ios::beg);

	bool wasfound=try_GetValue(infile,text,val);
	if( !wasfound )
	{
		wasfound=true;
		if( !optional)
			std::cerr<<"Field: "<<text<<" not found in .CFG file"<<" aborting "<<std::endl;

	}

	infile.clear(); 
	infile.seekg(curr_pos,std::ios::beg);
	
	return wasfound;
};

