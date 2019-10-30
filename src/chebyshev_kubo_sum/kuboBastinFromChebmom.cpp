// C & C++ libraries
#include <iostream>		/* for std::cout mostly */
#include <string>		/* for std::string class */
#include <fstream>		/* for std::ofstream and std::ifstream functions classes*/
#include <vector>		/* for std::vector mostly class*/
#include <complex>		/* for std::vector mostly class*/

// MKL LIBRARIES
#define MKL_Complex16 std::complex<double>
#include "mathimf.h"
#include "mkl.h"
#include <numeric>
#include <algorithm>

typedef std::complex<double> complex; 

struct GLFunctor
{
	GLFunctor(double _theta): theta(_theta) {};
    void operator()(double& m) const { m= cos(theta*m)/(M_PI*sin(theta));  }
	const double theta;
};

struct GRFunctor
{
	GRFunctor(double _theta):theta(_theta), I(complex(0,1)){};
    void operator()(complex& m) const
    { 
		const double 
		stheta= sin(theta),
		ctheta= cos(theta);
		m= I*exp(-I*theta*m)*( m*stheta - I*ctheta )/stheta/stheta/stheta;  
	}
	const complex I; 
	const double theta;
};

struct chebMom2D
{
	chebMom2D():numMom1(0),numMom0(0){};

	chebMom2D(const int m0,const int m1):numMom1(m1),numMom0(m0),mu( std::vector<complex>(numMom1*numMom0, 0.0) ){};
	 
	complex& operator()(const int m0,const int m1)
	{
		return mu[ m0*numMom1 + m1 ];
	}
	int numMom1,numMom0;
	std::vector< complex<double> > mu;
};

struct KuboFunctor
{
	//CONSTRUCTOR OF THE FUNCTOR
	KuboFunctor(const int numMom):
	numMom(_numMom),
	GL( std::vector< double >(numMom) ), 
	GR( std::vector< complex<double> >(numMom) )
	{};

	//THE FUNCTOR
    complex operator()(const double theta) 
    {
		const int M = GL.size();
		GLFunctor GLfun(theta); GRFunctor GRfun(theta);

		//Compute the operator at the left
		for(int m=0; m <numMom; m++)
			GL[m] = m;
		for_each(GL.begin(), GL.end(), GLfun);

		//Compute the operator at the left
		for(int m=0; m <numMom; m++)
			GR[m] = m;
		for_each(GR.begin(), GR.end(), GRfun);

		//Construct the Gamma_mn matrix
		const complex<double> sum;
		for(int m =0 ; m < numMom ; m++ )
		{
			const complex<double> partial_sum;
			void cblas_zdotu_sub (numMom,&GR[0],1,&(chebmom->mu[0]),1, partial_sum);
			sum+=partial_sum*GL[m];
		}
		return -sum/sin(theta);
	}

	std::vector< complex<double> > GR;
	std::vector< double> GL;
	complex<double> *chebMoms; 
	const int numMom;
	
};

int main(int argc, char *argv[])
{	
	if (argc != 7)
	{
//		printHelpMessage();
		return 0;
	}
//	else
//		printWelcomeMessage();

	const std::string
	momfilename = argv[1];
	double
	temperature = atof(argv[3]),
	Emin 		= atof(argv[4]),
	Emax 		= atof(argv[5]),
	dE   		= atof(argv[6]);

	std::ifstream momfile( momfilename.c_str() );

	double HalfWidth,BandCenter; int dim; 
	int numMoms0,numMoms1;
	momfile>>dim>>HalfWidth>>BandCenter; HalfWidth=HalfWidth/2.0;//in the file what you have is the bandwidth
	momfile>>numMoms0>>numMoms1;
	const int maxNumMom = ( (numMoms0 > numMom1) ? numMom0 : numMom1 ) ;
 	chebMom mu(maxNumMom,maxNumMom); double rmu,imu;
	for( int m0 = 0 ; m0 < numMoms0 ; m0++)
	for( int m1 = 0 ; m1 < numMoms1 ; m1++)
	{	
		momfile>>rmu>>imu;
		mu(m0,m1) = std::complex( rmu,imu);
	}
	const int maxNumMom = ( (numMoms0 > numMom1) ? numMom0 : numMom1 ) ;
	momfile.close();


	double lambda = maxNumMom *broadening/HalfWidth/1000.;
	
	//By performing the transformation x = Cos(theta)
	//all nodes of the matGamma matrix are equally spaced
	//on the interval (0,pi) and the separation is given by 1/M.
	//Hence, in this domain of integration, one would need at least
	// 10*M energy point to capture the oscilattion properly
	const int num_angles = 10*maxNumMom;
	vector< double >  angles(num_angles,0);
	
	for( int i=0; i < num_angles; i++)
	{
		const double 
		xbound = 0.9,
		theta_min = acos( xbound),
		theta_max = acos(-xbound),
		theta = theta_min + i*(theta_max-theta_min)/(num_angles-1) ;
		angles[i] = theta ;
	}
	
	
	std::string
	outputName  =momfilename+".OUT";
	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::ofstream outputfile( outputName.c_str() );
	KuboFunctor kuboFun( maxNumMom ); kuboFun.chebMoms = &mu;
	for( vector< double >::iterator it =  angles.begin();
									it!=  angles.end();
								    it++)
	{
		outputfile<<cos(*it)*HalfWidth + BandCenter <<" "<<kuboFun(*it)/HalfWidth/HalfWidth<<std::endl;
	}


return 0;
}		
		
		
cd