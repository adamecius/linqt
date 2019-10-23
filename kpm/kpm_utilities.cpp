
#include "kpm_utilities.hpp"

my::integer
kpm_util::GetMomentsPerNode( const my::integer numMoms)
{
	//Check if the number of processor is not too large
	assert (KPM_MPI_GETPROC() <= std::floor (0.5 * numMoms ) );

	return std::ceil( (my::real) numMoms / (my::real) KPM_MPI_GETPROC() );
}

void  kpm_util::GetMPIStatus(const my::integer numMoms)
{
	if( KPM_MPI_GETRANK() ==0 )
	std::cout<<"\nThe KPM-MPI enviroment is running "<<KPM_MPI_GETPROC()<<" processes."<<std::endl
			 <<"each node will calculate "<<kpm_util::GetMomentsPerNode(numMoms)<< " moments"<<std::endl
			 <<"for a grand total of "<<kpm_util::GetMomentsPerNode(numMoms)	*KPM_MPI_GETPROC()<< " moments"<<std::endl;
};

void  kpm_util::GetOMPStatus()
{
	if( KPM_MPI_GETRANK() ==0 )
	std::cout<<"\nThe KPM-OMP enviroment is using "<<KPM_OMP_NUM_THREADS<<" threads per process."<<std::endl;
};


std::string
kpm_util::GetNodeLabel()
{
  std::string sNodeID;          // string which will contain the result
  std::ostringstream convert;   // stream used for the conversion
  convert << KPM_MPI_GETRANK(); 	// insert the textual representation of 'Number' in the characters in the stream
  return "_proc"+convert.str (); 		// set 'Result' to the contents of the stream
}

