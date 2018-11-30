#ifndef SHARED_VALIDATION_HPP
#define SHARED_VALIDATION_HPP

// C++ libraries
#include <string>
#include <vector>


//LINUX-ONLY LIBRARIES
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#define ABORT_IF_FAILED(x) do { int stat = (x); if (FAILED(stat)) return stat; } while(0)

namespace validation
{
	bool CheckFolders(const size_t TotDirs, const std::string* dir_name);
	bool CreateFolder(const std::string dir_name);
	bool CheckFiles(const size_t TotDirs, const std::string* dir_name);


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>//
	bool CheckFolders(const size_t TotDirs, const std::string* dir_name)
	{	
		bool success=true;
		for(size_t i=0;i<TotDirs;i++)
		if(NULL == opendir(dir_name[i].c_str()) )
		{
			std::cout<<"Directory: "<<dir_name[i]<<" not found."<<std::endl;
			std::cout<<"Creating the directory"<<std::endl;
			success=success*CreateFolder(dir_name[i]);
			if( !success )
				std::cerr<<"Error creating directory: "<<dir_name[i]<<std::endl;
		}
		return success;
	};

	inline bool CreateFolder(std::string dir_name)
	{
//		const int dir_err = mkdir(dir_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
//		if (-1 == dir_err)
//			return false;
		return true;
	};


	inline bool CheckFiles(const size_t TotFiles, const std::string* file_name)
	{	
		bool success=true;
		struct stat buffer;   
		for(size_t i=0;i<TotFiles;i++)
		{
			success=success* (stat (file_name[i].c_str(), &buffer) == 0); 
			if( !success )
				std::cerr<<"The File "<<file_name[i]<<" was not found."<<std::endl;
		}
		return success;
	};

}


#endif
