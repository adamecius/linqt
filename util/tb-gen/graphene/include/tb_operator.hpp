

#ifndef TBOperator_HPP
#define TBOperator_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <limits>       // std::numeric_limits
#include "lattice_geometry.h"

struct mat_tuple
{
	int col;
	Complex val;

	//<<<<<<<<<<<<<<<CONSTRUCTOR>>>>>>>>>>>>>>>>>>>//
	mat_tuple():
	 col(0), val(0) {}

	mat_tuple(const int _col, Complex _val  ):
	 col(_col), val(_val) {}

};


struct by_col { 
    bool operator()(mat_tuple const &a, mat_tuple const &b) const { 
        return a.col < b.col ;
    }
};


class TBOperator
{
	public:
	//<<<<<<<<<<<<<<<CONSTRUCTORS>>>>>>>>>>>>>>>>>>>//
	TBOperator():
	dim_(0), max_coord(0),
	Op_triplet( std::vector< std::vector< mat_tuple> >() ) {}
	
	TBOperator(const int _dim, const int _max_coord):
	dim_(_dim), max_coord(_max_coord),
	Op_triplet( std::vector< std::vector< mat_tuple> >(_dim) )
	{
		for(int i=0; i< _dim ; i++) 
			Op_triplet[i].reserve( _max_coord );
	}

	//<<<<<<<<<<<<<<<PUBLIC METHODS>>>>>>>>>>>>>>>>>>>//
	void AddEntry(const int row,const int col,const Complex val )
	{
		if( row >=  Dimension() )
		{
			std::cerr<<" The row in AddEntry method is larger than Dimension "
					 <<" of the system"<<std::endl;
			std::exit(-1);
		}


		bool elemFound=false;
		if ( Op_triplet[row].size() == 0) 
		{
			Op_triplet[row].push_back(mat_tuple(col,val));
			elemFound =true;
		}
		else 
		for(int i=0;i < Op_triplet[row].size() ;i++)
		if( Op_triplet[row][i].col == col )
		{
			Op_triplet[row][i].val+= val ;
			elemFound =true;
			goto outTheLoop;
		}	
		outTheLoop:

		if (! elemFound )
			Op_triplet[row].push_back(mat_tuple(col,val));

		size_t currentNumElem= Op_triplet[row].size();
		if( currentNumElem > max_coord )
		{
			std::cerr<<	"The number of conections in row "
					 <<row<<" is "
					 <<currentNumElem
					 <<" and exceed or is equal to conectivity number"
					 <<max_coord<<" "<<std::endl;
			std::exit(-1);
		}
		
	}
	

	inline 
	int WriteIntoFile(std::string label ) const
	{
		Real NumZERO=std::numeric_limits<Real>::epsilon();

		std::string hamFilename = label + ".OP" ;
		std::ofstream hamFile( hamFilename.c_str(),  std::ofstream::binary );

		const size_t dim=Dimension();
		
		size_t NNZ=0;
		for(size_t row=0;row < dim ; row++)
		{
			const size_t elemPerRow = Op_triplet[row].size() ;
			for(size_t n=0;n < elemPerRow ; n++)
			//if( std::norm( Op_triplet[row][n].val ) > NumZERO )
				NNZ+=1;
		}

		hamFile	<<Dimension()<<" "<<NNZ<<" ";
		for(size_t row=0;row < dim ; row++)
		{
			const size_t elemPerRow = Op_triplet[row].size() ;
			
			//Get the number of non-zero elements in the row
			std::vector<mat_tuple> unsortedCol;

			size_t nnz=0;
			for(size_t n=0;n < elemPerRow ; n++)
			//if( std::norm( Op_triplet[row][n].val ) > NumZERO )
			{	
				unsortedCol.push_back(Op_triplet[row][n]);
				nnz+=1;
			}
			hamFile<<nnz<<" ";

			std::sort(unsortedCol.begin(), unsortedCol.end(), by_col());			

			for(int n=0;n < unsortedCol.size() ; n++)
				hamFile	<<unsortedCol[n].col<<" "
						<<unsortedCol[n].val.real()<<" "
						<<unsortedCol[n].val.imag()<<" ";
		}
		hamFile.close();

	}  


	inline 
	int saveTxtCoo(std::string label ) const
	{
		Real NumZERO=std::numeric_limits<Real>::epsilon();

		std::string hamFilename = label + ".COO" ;
		std::ofstream hamFile( hamFilename.c_str()  );

		
		for(int row=0;row < Dimension() ; row++)
		{
			std::vector<mat_tuple> unsortedCol;
			for(int n=0;n < Op_triplet[row].size() ; n++)
//			if( std::norm( Op_triplet[row][n].val ) > NumZERO )
			{
				mat_tuple elem;
				elem.col=Op_triplet[row][n].col;
				elem.val=Op_triplet[row][n].val;
				unsortedCol.push_back(elem);
			}
			std::sort(unsortedCol.begin(), unsortedCol.end(), by_col());


			for(int n=0;n < unsortedCol.size() ; n++)
				hamFile	<<row<<" "
						<<unsortedCol[n].col<<" "
						<<unsortedCol[n].val.real()<<" "
						<<unsortedCol[n].val.imag()<<std::endl;
		}
		hamFile.close();
	}

  

	
	//<<<<<<<<<<<<<<<GETTERS AND SETTERS>>>>>>>>>>>>>>>>>>>//
	inline 
	int Dimension() const
	{
		return dim_;
	}  
	
	
		std::vector< mat_tuple>&
	operator [] (const int i)
	{
		return Op_triplet[i];
	}

	public:
	std::vector< std::vector< mat_tuple> > Op_triplet;
	int dim_;
	
	private:
	int max_coord;
};








#endif


