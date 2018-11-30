
#ifndef DENSE_TBOP_H
#define DENSE_TBOP_H

#include <vector>
#include <string>
#include <fstream>
#include "kpm_dense_matrix.hpp"
#include "types_definitions.hpp"

class TBOperator
{
	
	public:
	TBOperator(void){};	

	TBOperator(std::string inputName)
	{
		readOperator(inputName);
	};
	
	void readOperator( std::string inputName)	
	{
		std::ifstream inputFile;
		inputFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		try  {	inputFile.open(inputName.c_str());	}
		catch (std::ifstream::failure e)
		{
			inputFile.close();
			std::cerr 	<<"Exception:"<<e.what()
						<<"\n while opening k-op  info sfile: "
						<< inputName<<"."<<std::endl;
			std::exit(-1);
		}
		inputFile>>OpDim>>numEntries;
		Op = kpm::complex_matrix(OpDim,OpDim);
		
		val = std::vector< kpm::complex > (numEntries);
		row = std::vector< int > (numEntries);
		col = std::vector< int > (numEntries);
		pos = std::vector< std::vector<double> > (numEntries);
		for(int n=0; n<numEntries ; n++)
			pos[n] = std::vector<double>(3);

		for(int n=0; n<numEntries ; n++)
		{
			kpm::real Reval,Imval;
			inputFile>>row[n]>>col[n]>>Reval>>Imval>>pos[n][0]>>pos[n][1]>>pos[n][2];
			val[n] = kpm::complex (Reval, Imval);
		}
		
		kpm::real k[3]={0,0,0};
		Addk_dependence( k );

	}

	void readLatticeInfo( std::string inputName)	
	{

		lat = std::vector< std::vector< double> >(3);
		rec = std::vector< std::vector< double> >(3);
		for(int i=0;i<3;i++)
		{
			lat[i] = std::vector< double>(3);
			rec[i] = std::vector< double>(3);
		}

		std::ifstream inputFile;
		inputFile.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		try  {	inputFile.open(inputName.c_str());	}
		catch (std::ifstream::failure e)
		{
			inputFile.close();
			std::cerr 	<<"Exception:"<<e.what()
						<<"\n while opening lattice info sfile: "
						<< inputName<<"."<<std::endl;
			std::exit(-1);
		}
		
		inputFile>>latconst;
		
		for(int m=0; m<3 ; m++)
		for(int n=0; n<3 ; n++)
		{
			inputFile>>lat[m][n];
			lat[m][n]=lat[m][n]*latconst;
		}
		inputFile.close();
		
		
		latVol =-lat[0][2]*lat[1][1]*lat[2][0] + 
				 lat[0][1]*lat[1][2]*lat[2][0] + 
				 lat[0][2]*lat[1][0]*lat[2][1] - 
				 lat[0][0]*lat[1][2]*lat[2][1] - 
				 lat[0][1]*lat[1][0]*lat[2][2] +
				 lat[0][0]*lat[1][1]*lat[2][2];


		rec[0][0]=-lat[1][2]*lat[2][1] + lat[1][1]*lat[2][2]; 
		rec[0][1]= lat[1][2]*lat[2][0] - lat[1][0]*lat[2][2]; 
		rec[0][2]=-lat[1][1]*lat[2][0] + lat[1][0]*lat[2][1];

		rec[1][0]= lat[0][2]*lat[2][1] - lat[0][1]*lat[2][2]; 
		rec[1][1]=-lat[0][2]*lat[2][0] + lat[0][0]*lat[2][2]; 
		rec[1][2]= lat[0][1]*lat[2][0] - lat[0][0]*lat[2][1];

		rec[2][0]=-lat[0][2]*lat[1][1] + lat[0][1]*lat[1][2]; 
		rec[2][1]= lat[0][2]*lat[1][0] - lat[0][0]*lat[1][2]; 
		rec[2][2]=-lat[0][1]*lat[1][0] + lat[0][0]*lat[1][1];

		for(int m=0; m<3 ; m++)
		for(int n=0; n<3 ; n++)
			rec[m][n]=rec[m][n]*2.0*M_PI/latVol;
	}

	kpm::complex_matrix Mat() 
	{
		return Op; 
	}

	kpm::complex operator()( const int i, const int j) const
	{
		return Op(i,j); 
	}

	double Lat( const int i, const int j) const
	{	
		return lat[i][j]; 
	}

	double LatConst() const { return latconst; }
	double Rec( const int i, const int j) const
	{
		return rec[i][j]; 
	}
	
	double LatVol() const {return latVol;}

	inline
	void Addk_dependence( const kpm::real* k )
	{
		kpm::complex x=0;
		for( int i=0; i < OpDim; i++)
		for( int j=0; j < OpDim; j++)
			Op(i,j) = 0.0;

		for( int n=0; n < numEntries ; n++ )
		{
			const kpm::real dotv= pos[n][0]*k[0] + pos[n][1]*k[1] + pos[n][2]*k[2];
			kpm::complex expv= kpm::complex( cos(dotv) , -sin(dotv) );
			kpm::complex v= val[n];
			kpm::complex O= kpm::complex( v.real()*expv.real() - v.imag()*expv.imag() , v.real()*expv.imag() + v.imag()*expv.real() );
			Op(row[n],col[n])+=O;
		}
	}

	inline
	void Multiply(	const kpm::real a,const kpm::complex* X, 
					const kpm::real b, kpm::complex* Y) const 
	{
		for( int i=0; i < OpDim; i++)
		{
			kpm::complex t= 0.0;
			for( int j=0; j < OpDim; j++)
			{
				const kpm::complex O=Op(i,j);
				const kpm::complex x=X[j];
				t.real(t.real() + x.real()*O.real() - x.imag()*O.imag()) ;
				t.imag(t.imag() + x.real()*O.imag() + x.imag()*O.real()) ;
			}
			kpm::complex y=Y[i];
			y.real( b*y.real() + a*t.real());
			y.imag( b*y.imag() + a*t.imag());
			Y[i]= y;
		}	
	}

	void Rescale(	const kpm::real a)
	{
		for( int i=0; i < val.size(); i++)
			val[i]=a*val[i];
		kpm::real k[3]={0,0,0};
		Addk_dependence( k );

	}


	int Dimension() const { return OpDim;}

	private:
	kpm::complex_matrix Op;
	int OpDim,  numEntries;
	std::vector< kpm::complex > val;
	std::vector< int > row, col;
	std::vector< std::vector<double> > pos;
	std::vector< std::vector<double> > lat;
	std::vector< std::vector<double> > rec;
	double latVol, latconst;
};

#endif
