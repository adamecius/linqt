#ifndef MPI_OFILE
#define MPI_OFILE

#include <vector>
#include <string>
#include <fstream>


void concatenate_moments(const int Mom,std::vector<std::string> momName_array, 
						 std::string outputname)
{
		
		const kpm::complex NULLVAL=02051988.; //defines NULL VALUE
		std::string header;	//Variable to store the header of the moment file
		std::vector<kpm::complex> mu2D(Mom*Mom,NULLVAL); //moment array

		const int numFiles=momName_array.size();
		for(int id = 0 ; id< numFiles; id ++)
		{
			std::ifstream mom2Dtmp( momName_array[id].c_str()  );
			if ( id == 0 )
				std::getline (mom2Dtmp,header);
			while( !mom2Dtmp.eof() )
			{
				int m,n;
				kpm::real ReVal, ImVal;
				mom2Dtmp>>m>>n>>ReVal>>ImVal;
				mu2D[m*Mom+n] = kpm::complex(ReVal,ImVal);
			}
			mom2Dtmp.close();			
		}

		std::ofstream mom2D( outputname.c_str()  );
		mom2D<<header<<std::endl;
		for(int m = 0 ; m < Mom ; m++ )
		for(int n = 0 ; n < Mom ; n++ )
		if ( mu2D[m*Mom+n] == NULLVAL) 
			std::cerr<<"The moment mu("<<m<<","<<n<<") was not calculated. Value="<< mu2D[m*Mom+n]<<std::endl;
		else
			mom2D<<m<<" "<<n<<" "<<mu2D[m*Mom+n].real()<<" "<<mu2D[m*Mom+n].imag()<<std::endl;
		mom2D.close();			
};

#endif




