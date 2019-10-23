

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_01.hpp>
#include <ctime>
#include <sys/time.h>
#include <unistd.h>


int
IndexesToIndex(	const int i0, const int dim0,
				const int i1, const int dim1,
				const int i2, const int dim2,
				const int io, const int norb,
				const int is, const int nspin)
				{
					return ( ( ( (i0+dim0)%dim0*dim1 + (i1+dim1)%dim1 )*dim2 + (i2+dim2)%dim2 )*norb + io )*nspin + is;
				}

int main (int argc, char** argv)
{

	if(argc != 4+1 )
	{
		std::cerr<<" The function :"<<argv[0]<<" should be called with the following parameters: "<< std::endl
				 <<" config_file Puddle_Height[eV] Puddle_Concentration[0-1] Puddle_Range[lat_unit]."<<std::endl;
		return 0;
	}

	//Read the name of the configuration file
	const std::string config_filename(argv[1]);
	const double puddHeight=atof(argv[2]);
	const double puddConcen=atof(argv[3]);
	const double puddRange =atof(argv[4]);


	std::ifstream config_file(config_filename.c_str());
	///Look for the line lattice_info
	std::string line;
	bool found=false;
	for (int line_num = 0; std::getline(config_file, line); ++line_num)
		if(line=="lattice_info")
		{
			found=true;
			break;
		}
	if(!found)
	{
		config_file.close();
		std::cerr<<"The config file does not posses a disorder_info section. The simulation cannot proceed"<<std::endl;
		return 0;
	}

		std::string label;
		int OrbitalNum=0;
		int SpinNum=0;
		int coord=0;
		std::vector<int> SitesInDir(3);
		std::vector < std::vector<double> > lat(3);
		for(int i=0;i<3;i++)
			lat[i]=std::vector<double>(3,0);

		config_file	>>label
					>>lat[0][0]>>lat[0][1]>>lat[0][2]
					>>lat[1][0]>>lat[1][1]>>lat[1][2]
	                >>lat[2][0]>>lat[2][1]>>lat[2][2]
	                >>SitesInDir[0]>>SitesInDir[1]>>SitesInDir[2]
				    >>OrbitalNum>>SpinNum>>coord;

		const double Delta[] = { (lat[0][0]+lat[1][0])/3. , (lat[0][1]+lat[1][1])/3. , 0 };			// THIRD LATTICE VECTOR


		std::string output_name=label+"PuddleHeight"+
								argv[2]+"eVConc"+argv[3]+
								"Range"+argv[4]+"nm.dis";

		std::ofstream output_file(output_name.c_str());


		const double epsilon=std::numeric_limits<double>::epsilon();
		if(puddHeight <= epsilon || puddConcen  <= epsilon)
		{
			output_file.close();
			return 0;
		}


		boost::random::mt19937 rng;         // produces randomness out of thin air
		boost::random::uniform_01<double> dice;
		boost::random::uniform_real_distribution<double> RandomEnergy( -1.0 , 1.0);


		///Number of orbitals in the simulation
		const int TotalOfOrbitals=SitesInDir[0]*SitesInDir[0]*SitesInDir[2]*OrbitalNum*SpinNum;

		///Stimated number of impurities
		const int  estImpNum= (1.3)*puddConcen * TotalOfOrbitals/SpinNum;//The two is to overstimate the value
		std::vector< std::vector< double > > impPosition;
		impPosition.reserve( estImpNum );
		std::cout<<"The estimated number of puddles is set as: "<<estImpNum<<std::endl;


		std::vector< double >
		this_imp_pos(3);										//a singl impurity position

		int impurity_count= 0;


		// first we go through all the possible positions
		for(int i2=0; i2 < SitesInDir[2]  ; i2++ )
			for(int i1=0; i1 < SitesInDir[1]  ; i1++ )
				for(int i0=0; i0 < SitesInDir[0] ; i0++ )
					for(int io=0; io < OrbitalNum ; io++ )
						if( dice(rng) < puddConcen )
						{
							this_imp_pos[0]=lat[0][0]*i0 + lat[1][0]*i1 + lat[2][0]*i2 + io*Delta[0];
							this_imp_pos[1]=lat[0][1]*i0 + lat[1][1]*i1 + lat[2][1]*i2 + io*Delta[1];
							this_imp_pos[1]=lat[0][2]*i0 + lat[1][2]*i1 + lat[2][2]*i2 + io*Delta[2];
							impPosition.push_back(this_imp_pos);
						}

		int counter=0;
		std::vector<double> diagElements( TotalOfOrbitals ,0 );
		//We incorporate the onsite energy of the puddles.
		for(int i2=0; i2 < SitesInDir[2]  ; i2++ )
			for(int i1=0; i1 < SitesInDir[1]  ; i1++ )
				for(int i0=0; i0 < SitesInDir[0] ; i0++ )
					for(int io=0; io < OrbitalNum ; io++ )
					{
						for(int it0=-1; it0 <= 1; it0++ )	// tile index 1
							for(int it1=-1; it1 <= 1; it1++ )	// tile index 2. This is used to take into account the periodic boundary condition
							{
								const int I0=i0+SitesInDir[0]*it0 ;
								const int I1=i1+SitesInDir[1]*it1 ;
								const double
								r[2]={
										lat[0][0]*I0+lat[1][0]*I1+io*Delta[0],
										lat[0][1]*I0+lat[1][1]*I1+io*Delta[1]
								};
								for(int  imp=0;imp< impPosition.size(); imp++)
								{
									const double
									dist = 	( r[0]-impPosition[imp][0] )*( r[0]-impPosition[imp][0] ) +
											( r[1]-impPosition[imp][1] )*( r[1]-impPosition[imp][1] );

									const double
									Ei=RandomEnergy(rng)*puddHeight*exp(-0.5*dist/puddRange/puddRange);
									//Created for periodic puddles
									if( std::abs(Ei) > epsilon )
										for(int is=0; is <SpinNum   ; is++ )
										{
											//compute the site where the energy is going to be changed
											const int
											k0= IndexesToIndex(	I0, SitesInDir[0],
																I1, SitesInDir[1],
																i2, SitesInDir[2],
																io, OrbitalNum,
																is, SpinNum);
											diagElements[k0]=diagElements[k0]+Ei;
										}
								}
							}
						++counter;
						const int checkpoint= int((double)TotalOfOrbitals/SpinNum/100);
						if(counter%checkpoint == 0 )
							std::cout<<"Setting the puddles, there are :"<<counter<<"/"<<TotalOfOrbitals/SpinNum<<" left."<<std::endl;
						}
		// std::cout<<"Setting the hamiltonian in the sparse matrix"<<std::endl;

		double mean=0;
		double stdDev=0;

		for(int k0=0;k0<diagElements.size();k0++)
			if(std::abs(diagElements[k0])>epsilon)
			{
				output_file<<k0<<" "<<diagElements[k0]<<std::endl;
				mean=mean+diagElements[k0];
				stdDev= stdDev +diagElements[k0]*diagElements[k0];
			}

		mean=mean/diagElements.size();
		stdDev=stdDev/diagElements.size();
		stdDev=sqrt(stdDev- mean*mean);

/*		for(int k0=0;k0<diagElements.size();k0++)
			if(std::abs(diagElements[k0])>epsilon)
				stdDev=stdDev+pow(mean-diagElements[k0],2.);
		stdDev=sqrt(stdDev/diagElements.size());
 */
		std::cout<<"The puddle have mean:"<<mean<<" StdDev: "<<stdDev<<std::endl;


	output_file.close();
 return 0;}
