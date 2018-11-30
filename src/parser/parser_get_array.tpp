



//This function tries to get the requested value, if it fails it return false
template <typename T>
bool 
parser::try_GetArray(std::ifstream& infile, const std::string tag, T* val,int dim, std::string header="")
{
	//transform the tag into its lowercap form
	std::string lower_tag=tag;
	std::transform(tag.begin(), tag.end(), lower_tag.begin(), ::tolower);

	//Look at the beginning of the file
	infile.clear();
	infile.seekg(0, std::ios::beg);
	if(header!="")
		FindHeader(infile,header);
	

	std::string	line; //variable to store each line of the file
	if ( infile ) // Check if file is valid
	while	(  getline(infile , line ) )
	{
		size_t 
		equ_pos,txt_pos, com_pos; //if not comment is found return the final line
	
		//tolower everything before equal
		if( (equ_pos=line.find('=', 0) )==std::string::npos ) equ_pos=line.size();
		std::transform(line.begin(), line.begin()+equ_pos, line.begin(), ::tolower);
		//Remove everything after a comment
		if( (com_pos=line.find('#', 0) )==std::string::npos ) com_pos=line.size();
		line=line.substr (0,com_pos); 
		//Check for the tag(everything in lower)
		if	(
				( txt_pos=line.find(lower_tag, 0) )	!= std::string::npos && 
				txt_pos<equ_pos
			)
		{
			//if found tag, keep everything after=
			std::string str 		= TrimLine(line.substr(equ_pos+1,line.size()));  
			infile.clear();
			infile.seekg(0, std::ios::beg);
			std::vector<T> vec_val;
			bool isSuccess=str2vector(str,vec_val);
			isSuccess*=(bool)( vec_val.size()== dim );
			if( isSuccess )
				for(int i=0;i< dim; i++) val[i]=vec_val[i];
			else
				std::cerr	<<"One of the parameters could not be read."<<std::endl
							<<"Or the dimension of the vector:\t"<<vec_val.size()<<std::endl
							<<"is not the same as of the array:\t"<<dim<<std::endl;

			return isSuccess;
			
		}
	}
	
	infile.clear();
	infile.seekg(0, std::ios::beg);
	return false;

};


template <typename T>
bool 
parser::GetArray(std::ifstream& infile, const std::string text, T* val,int dim, std::string header, bool  optional=false  )
{
	bool wasfound=try_GetArray(infile,text,val,dim,header);
	if( !wasfound )
	{
		wasfound=true;
		if( !optional)
		{
			if(IS_WORLD_ROOT())
				std::cerr<<"Field: "<<text<<" not found in .CFG file"<<" aborting "<<std::endl;
			SAFE_MPI_FINALIZE();
		}
	}

	infile.clear(); 
	infile.seekg(std::ios::beg);
	
	return wasfound;


};


