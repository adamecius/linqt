

#ifndef TBOperator_HPP
#define TBOperator_HPP

#include "types_definitions.hpp"

struct mat_tuple
{
	kpm::integer col;
	kpm::complex val;

	//<<<<<<<<<<<<<<<CONSTRUCTOR>>>>>>>>>>>>>>>>>>>//
	mat_tuple():
	 col(0), val(0) {}
}

class TBOperator()
{

	//<<<<<<<<<<<<<<<CONSTRUCTORS>>>>>>>>>>>>>>>>>>>//
	TBOperator():
	dim_(0), max_coord(0)
	Ham_triplet( std::vector< std::vector< mat_tuple> > >() ) {}
	
	TBOperator(const int _dim, const int _max_coord):
	dim_(_dim), max_coord(_max_coord)
	Ham_triplet( std::vector< std::vector< mat_tuple> > >(_dim) ) { }
	

	//<<<<<<<<<<<<<<<PUBLIC METHODS>>>>>>>>>>>>>>>>>>>//
	void ReadElemsFromFile(const std::string filename)
	{
		//maximum possible number of entries
		const int max_num_entries= Dimensions()*max_coord;
		

		std::ifstream HamFile  (filename);	//note, I am not checking if the file is open
		std::ifstream CoordFile(filename);	//note, I am not checking if the file is open

		
		for( int m=0; m < Dimension() ; ++m )
		{
			int coord_row,coord_num;
			CoordFile>>coord_row>>coord_num;

			Ham_triplet[row]= std::vector< mat_tuple>(coord_num);

			for( int n=0; coord_num ; ++n )
			{
				kpm::integer col,row;
				kpm::real Reval, Imval;
				inputFile>>row>>col>>Reval>>Imval;
				
				if( row != coord_row)
				{
					std::ceer<<"Problem matching the Hamiltonian file";
					std::ceer<<"with the coordination file"<<std::endl;
					std::exit(-1);
				}
				
				Ham_triplet[row][n].col=col;
				Ham_triplet[row][n].val=val;				
			}
		}
		HamFile.close();
		CoordFile.close();	
	}

	
	//<<<<<<<<<<<<<<<GETTERS AND SETTERS>>>>>>>>>>>>>>>>>>>//
	inline 
	kpm::integer Dimension() const
	{
		return dim_;
	}  
	
	public:
	std::vector< std::vector< mat_tuple> > > Ham_triplet;
	std::vector< mat_tuple > Ham_triplet;
	kpm::integer dim_;
	
	private
	kpm::integer max_coord;
}





#endif
