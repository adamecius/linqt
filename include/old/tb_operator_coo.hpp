

#ifndef TBOperator_HPP
#define TBOperator_HPP

#include <mpi.h>
#include <omp.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include "types_definitions.hpp"
#include "parser.hpp"

#include "kpm_memory.hpp" //aligned_malloc,aligned_free



struct mat_tuple
{
	kpm::integer col;
	kpm::complex val;

	//<<<<<<<<<<<<<<<CONSTRUCTOR>>>>>>>>>>>>>>>>>>>//
	mat_tuple():
	 col(0), val(0) {}
};

class TBOperator
{
	public:
	//<<<<<<<<<<<<<<<CONSTRUCTORS>>>>>>>>>>>>>>>>>>>//
	TBOperator():
	dim_(0), max_coord(0),
	Op_triplet( std::vector< std::vector< mat_tuple> >() ) {};
	
	TBOperator(const kpm::integer _dim, const kpm::integer _max_coord):
	dim_(_dim), max_coord(_max_coord),
	Op_triplet( std::vector< std::vector< mat_tuple> >(_dim) ) { }
	

	//<<<<<<<<<<<<<<<PUBLIC METHODS>>>>>>>>>>>>>>>>>>>//
	void ReadElemsFromFile(const std::string opLabel)
	{
		//Cheack if the file existas and open the file
		std::string
		localFilename= opLabel;	
	
		if ( ! fileExists(localFilename) )
		{
			std::cerr<<"File: "<<localFilename<<" not found, return -1;";
			std::exit(-1);
		}
		std::ifstream OpFile  (localFilename.c_str(), std::ofstream::binary);
	
		//Read the dimension of the matrix, and the total number 
		size_t  OpDim,NNZ=0,ElemsPerRow=0 ;
		OpFile>>OpDim; SetDim(OpDim);
		

		aligned_malloc( (void**)&pntrb ,(OpDim+1)*sizeof(*pntrb),ALIGN); 


		Op_triplet= std::vector< std::vector< mat_tuple> >(OpDim);
		for( size_t row = 0; row < OpDim; ++row )
		{
			pntrb[row]=NNZ;
			OpFile>>ElemsPerRow;  NNZ+=ElemsPerRow;

			Op_triplet[row]=std::vector< mat_tuple>(ElemsPerRow);
			for( size_t elem=0; elem < ElemsPerRow ; ++elem )
			{
				size_t col; kpm::real Reval, Imval;
				OpFile>>col>>Reval>>Imval;
				mat_tuple t; t.col=col; t.val = kpm::complex(Reval,Imval);
				Op_triplet[row][elem]= t;
			}		
		}
		pntrb[OpDim]=NNZ;


		aligned_malloc( (void**)&val ,NNZ*sizeof(*val),ALIGN); size_t idx=0;
		aligned_malloc( (void**)&col ,NNZ*sizeof(*col),ALIGN);
		for( size_t row = 0; row < OpDim; ++row )
		{
			ElemsPerRow= Op_triplet[row].size();
			for( size_t elem=0; elem < ElemsPerRow ; ++elem )
			{
				col[idx]=Op_triplet[row][elem].col;
				val[idx]=Op_triplet[row][elem].val;
				idx+=1;		
			}
		}
		
		OpFile.close();
	};

	//<<<<<<<<<<<<<<<PUBLIC METHODS>>>>>>>>>>>>>>>>>>>//
	void ReadOpAsCSR(const std::string opLabel)
	{
		//Cheack if the file existas and open the file
		std::string
		localFilename= opLabel;	
	
		if ( ! fileExists(localFilename) )
		{
			std::cerr<<"File: "<<localFilename<<" not found, return -1;";
			std::exit(-1);
		}
		std::ifstream OpFile  (localFilename.c_str(), std::ofstream::binary);
	
		//Read the dimension of the matrix, and the total number 
		size_t  OpDim,NNZ=0,ElemsPerRow=0,idx=0;
		OpFile>>OpDim>>NNZ; SetDim(OpDim);
		

		aligned_malloc( (void**)&pntrb ,(OpDim+1)*sizeof(*pntrb),ALIGN); 
		aligned_malloc( (void**)&val ,NNZ*sizeof(*val),ALIGN); 
		aligned_malloc( (void**)&col ,NNZ*sizeof(*col),ALIGN); NNZ=0;
		

		for( size_t row = 0; row < OpDim; ++row )
		{
			pntrb[row]=NNZ;
			OpFile>>ElemsPerRow;  NNZ+=ElemsPerRow;
			for( size_t elem=0; elem < ElemsPerRow ; ++elem )
			{
				size_t _col; kpm::real _Reval, _Imval;
				OpFile>>_col>>_Reval>>_Imval;
				col[idx]=_col;
				val[idx]=kpm::complex(_Reval,_Imval);;
				idx+=1;		
			}
		}
		pntrb[OpDim]=NNZ;

		OpFile.close();
	};
	
	
	
	void Multiply(	const kpm::real _a, 
					const kpm::complex* x, 
					const kpm::real _b,
					kpm::complex* y ) const 
	{
		const int dim=Dimension();
		const kpm::real a=_a,b=_b;
		size_t row,idx;
		kpm::complex ytmp;
		int chunk = 1000;
		
		#pragma omp parallel for private(row,idx,ytmp) firstprivate(a,b)   schedule(static,chunk) 
		for( row=0; row < dim ; row++ )
		{
			ytmp = b*y[row] ;
			for( idx=pntrb[row]; idx < pntrb[row+1] ; idx++ )
			{
				const kpm::complex val_ =  val[idx];
				const kpm::integer col_ =  col[idx];
				ytmp+= x[col_]*a*val_;
			}
			y[row]=ytmp;
		}
		
		 
	}			
	


	void Rescale(	const kpm::real _scal )   
	{
		const kpm::real scal=_scal;
		const int dim=Dimension();
		size_t row,idx;

		#pragma omp parallel for firstprivate(scal) private(row,idx)
		for( row=0; row < dim ; row++ )
			for( idx=pntrb[row]; idx < pntrb[row+1] ; idx++ )
				val[idx]*=scal;
			
	}			
	
	//<<<<<<<<<<<<<<<GETTERS AND SETTERS>>>>>>>>>>>>>>>>>>>//
	size_t SetDim(size_t _dim) 
	{
		dim_=_dim;
	}  


	size_t Dimension() const
	{
		return dim_;
	}  
	size_t Dim() const
	{
		return dim_;
	}  
	
	//
	kpm::complex get (const kpm::integer row, const kpm::integer col)
	{

		for( kpm::integer n=0; n<Op_triplet[row].size() ; ++n )
		{

			if( Op_triplet[row][n].col==col ) 
				return Op_triplet[row][n].val ;				
		}
		
		std::cout<<"The element was not found, return 0"<<std::endl;
		return kpm::complex(0.) ;
					
	}
	
	public:
	std::vector< std::vector< mat_tuple> > Op_triplet;
	kpm::integer dim_;
	kpm::complex *val;
	kpm::integer *pntrb, *col;	


	kpm::integer *pntrb, *col;	


	private:
	kpm::integer max_coord;

};





#endif
