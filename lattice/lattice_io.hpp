/*
 * lattice_io.hpp
 *
 *  Created on: 05/09/2016
 *      Author: jgarcia
 */

#ifndef LATTICE_LATTICE_IO_HPP_
#define LATTICE_LATTICE_IO_HPP_

#include <iostream>
#include <string>

#include "mpi_util.hpp"
#include "lattice.hpp"

namespace NumCal
{
	namespace lattice_io
	{

		void PrintInfFromCfg(std::string _config_filename, NumCal::Lattice _lattice)
		{
			NumCal::cout<<"-------------READING FROM INPUT FILES-------------"<<NumCal::endl<<NumCal::endl;

			NumCal::cout<<"The following fields were read from the file: "<<_config_filename<<NumCal::endl<<NumCal::endl;

			NumCal::cout<<"\tSimulation Label: \""<<_lattice.Label()<<"\""<<NumCal::endl;

			for(short v=0;v<3 ;v++)
			{
				NumCal::cout<<"\tLattice Vector "<<v<<" : ( ";
				for(short i=0;i<2 ;i++)
					NumCal::cout<<_lattice.LatticeVector(v)[i]<<",\t";
				NumCal::cout<<_lattice.LatticeVector(v)[2]<<" )"<<NumCal::endl;

				NumCal::cout<<"\tNumber of UnitCells in this direction: "<<_lattice. CellsInDir(v)<<NumCal::endl;
			}

			NumCal::cout<<"\tNumber of orbitals per UnitCell: "<<_lattice.OrbitalNumber()<<NumCal::endl;
			NumCal::cout<<"\tThe spin of the system: "<<_lattice.OrbitalNumber()<<NumCal::endl<<NumCal::endl;

			NumCal::cout<<"-------------CALCULATED FROM THIS INPUTS-------------"<<NumCal::endl<<NumCal::endl;

			NumCal::cout<<"\tThe volume of the sample: "<<_lattice.Volume()<<NumCal::endl;

			for(short v=0;v<3 ;v++)
			{
				NumCal::cout<<"\tReciprocal Lattice Vector "<<v<<" : ( ";
				for(short i=0;i<2 ;i++)
					NumCal::cout<<_lattice.ReciprocalLatticeVector(v)[i]<<",\t";
				NumCal::cout<<_lattice.ReciprocalLatticeVector(v)[2]<<" )"<<NumCal::endl;
			}
			NumCal::cout<<NumCal::endl;



		};
	};
};
#endif /* LATTICE_LATTICE_IO_HPP_ */
