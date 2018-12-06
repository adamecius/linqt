template <typename T>
bool 
parser::GetBlock(std::ifstream& infile, const std::string tag, std::vector< std::vector<T> >& block,const int dim0,const int  dim1,  bool  optional=false  )
{
	infile.clear();
	infile.seekg(0, std::ios::beg);
	std::string begBlock="BEGIN_"+tag;
	std::string endBlock="END_"+tag;
	std::string	line; //variable to store each line of the file

	bool found_beg=false,found_end =false;
	try
	{
		if ( infile ) // Check if file is valid
		while(!found_beg && getline(infile,line) )
		{
			found_beg = line.find(begBlock, 0)!= std::string::npos;
			while( found_beg && !found_end && getline(infile,line) )
			{
				found_end = line.find(endBlock, 0)!= std::string::npos;

				if ( !found_end )
				{
					std::vector<T> vec_val;	
					std::string str	= TrimLine(line);  
					if( !qt::str2vector(str,vec_val) ||  vec_val.size()!= dim1 ) 
						return false;
					//Pass the value to the block
					block.push_back(vec_val);				
				}
			}
		}
	}
	catch(std::ifstream::failure e)
	{//Reset to the begining of the file
		infile.clear();
		infile.seekg(0, std::ios::beg);
	}

	if( block.size() != dim0 )
		return false;
							//reset the value before procedding
			return true;	

	return true;
};
