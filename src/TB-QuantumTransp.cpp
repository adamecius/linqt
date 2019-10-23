
#include <string>
#include <iostream>
#include <fstream>


#include "types_definitions.hpp"
#include "lattice.hpp"
#include "lattice_io.hpp"

#include "kpm.hpp"
#include "kpm_parallel.hpp"
#include "onsite_disorder.hpp"

#include <sys/time.h>
#include "random.hpp"
#include "string_utilities.hpp"

int main (int argc, char** argv)
{
	KPM_PARALLEL_INIT();

	if(argc != 2 )
	{
		NumCal::cerr<<" Please submit a config file, and only a config file"<<NumCal::endl;
		return 0;
	}

	//Read the name of the configuration file
	const std::string config_filename(argv[1]);

	//Initialize the Lattice class using the configuration file
	NumCal::Lattice my_lattice( config_filename ); //Checked

	//Print the information read from the configuration file
	NumCal::lattice_io::PrintInfFromCfg(config_filename,my_lattice); //checked

	//Read the hopping list from the hopping file
	my_lattice.ReadHoppingList(); //checked

	NumCal::integer
	estimated_hoppings= my_lattice.TotalOfOrbitals()+10;
	my_lattice.ReserveIrrHamSpace( estimated_hoppings );
	NumCal::cout<<"Number of stimated entries for sparse matrix : " <<estimated_hoppings<< NumCal::endl;
	NumCal::cout<<"Reserved memory for sparse matrix : " <<estimated_hoppings*sizeof(NumCal::scalar)/pow(2,20)<<" MB."<< NumCal::endl;


	std::cout<<"LOADING THE PUDDLES IN she_gapWS2N2048M3072PuddleUp2.8xp0.001eVRp3.dis"<<std::endl;
	std::ifstream puddle_file( ( my_lattice.Label()+".dis" ).c_str());


	std::cout<<" "<<estimated_hoppings<<" "<<my_lattice.TotalOfOrbitals()<<" "<<std::endl;

	int counter=0;
	while ( !puddle_file.eof()  || my_lattice.TotalOfOrbitals() == counter)
	{

//	for(int k0=0;k0<estimated_hoppings/2;k0++){

		NumCal::integer k0;
		NumCal::real Ei;
//		NumCal::real Ei= 2.*( (my::real)rand()/(my::real)RAND_MAX -0.5 );

		puddle_file>>k0>>Ei;
//		std::cout<<k0<<" "<<estimated_hoppings<<" "<<my_lattice.TotalOfOrbitals()<<" "<<Ei<<std::endl;

		for(my::integer is=0; is < my_lattice.SpinNumber()    ; is++ )
				my_lattice.IrregHam().SetHopping( IrregularHamiltonian::Hopping (k0,k0,Ei) );
	counter++;
	}
	puddle_file.close();


	my_lattice.SetTotalHamiltonian(); //This has to come before rescale

	NumCal::cout<<NumCal::endl<<"-------------BEGINNING OF KPM SECTION-------------"<<NumCal::endl<<NumCal::endl;
	NumCal::Kpm	kpm(my_lattice.TotalOfOrbitals(),config_filename);
	NumCal::cout<<"Truncation order: "<<kpm.TruncOrder()<<NumCal::endl;
	NumCal::cout<<"Energy Bounds: [ "<<kpm.MinBound()<<" , "<<kpm.MaxBound()<<" , "<<kpm.CutOff()<<" ]"<<NumCal::endl;
	NumCal::cout<<"Rescaling the regular and irregular part of the hamiltonian"<<NumCal::endl;
	my_lattice.RescaleHamiltonian(kpm.MinBound() , kpm.MaxBound() ,kpm.CutOff() );

	///The first Calculation will be a density
	///of states
	std::string momDOSOutput;
	momDOSOutput = "tmp/"+my_lattice.Label()+".momDOS.tmp";
//	kpm.DOSMoments(my_lattice,momDOSOutput);

	std::vector<std::string> momcondOutputFile (2);
	momcondOutputFile[0] = "tmp/"+my_lattice.Label()+".momxxcond.tmp";
	momcondOutputFile[1] = "tmp/"+my_lattice.Label()+".momxycond.tmp";



	kpm.SpinConductivityMoments(my_lattice,momcondOutputFile[0], momcondOutputFile[1]);




	KPM_PARALLEL_FINALIZE();


	return 0;
};
