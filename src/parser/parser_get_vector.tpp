



//This function tries to get the requested value, if it fails it return false
template <typename T>
bool 
parser::GetVector(std::ifstream& infile, const std::string tag,std::vector<T>& val,int dim, bool  optional=false  )
{
	//transform the tag into its lowercap form
	std::string lower_tag=tag;
	std::transform(tag.begin(), tag.end(), lower_tag.begin(), ::tolower);

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
			//reset the value before procedding
			infile.clear();
			infile.seekg(0, std::ios::beg);

			//if found tag, keep everything after=
			std::string str 		= TrimLine(line.substr(equ_pos+1,line.size()));  
			std::vector<T> vec_val;	//vector to store the array

			//Transform the vector and return false if fail
			//or the dimension is incorrect
			if( !qt::str2vector(str,vec_val) ||  vec_val.size()!= dim ) 
				return false;
			//Pass the vall
			for(int i=0;i<dim;i++)
				val[i] = vec_val[i];
			return true;	
		}
	}
	
	infile.clear();
	infile.seekg(0, std::ios::beg);
	return false;

};


