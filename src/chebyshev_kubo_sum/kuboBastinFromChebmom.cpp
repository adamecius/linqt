
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

const std::complex<double> I(0,1); 


void printHelpMessage();

void printWelcomeMessage();



struct GLFunctor
{
	GLFunctor(double _theta): theta(_theta) {};
    void operator()(double& m) const { m= cos(theta*m)/(M_PI*sin(theta));  }
	const double theta;
};

struct GRFunctor
{
	GRFunctor(double _theta): theta(_theta){};
    void operator()(std::complex<double>& m) const
    { 
		const double 
		stheta= sin(theta),
		ctheta= cos(theta);
		m= exp(-I*theta*m)*( m*stheta - I*ctheta )/stheta/stheta/stheta;  
	}
	const double theta;
};

struct chebMom2D
{
	chebMom2D():numMom1(0),numMom0(0){};

	chebMom2D(const int m0,const int m1):numMom1(m1),numMom0(m0),mu( std::vector<std::complex<double> >(numMom1*numMom0, 0.0) ){};
	 
	std::complex<double>& operator()(const int m0,const int m1)
	{
		return mu[ m0*numMom1 + m1 ];
	}
	int numMom1,numMom0;
	std::vector< std::complex<double> > mu;
};

struct KuboFunctor
{
	//CONSTRUCTOR OF THE FUNCTOR
	KuboFunctor(const int _numMom):
	numMom(_numMom),
	GL( std::vector< double >(_numMom) ), 
	GR( std::vector< std::complex<double> >(_numMom) )
	{};

	//THE FUNCTOR
    double operator()(const double theta) 
    {
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
		std::complex<double> sum=0.0;
		for(int m =0 ; m < numMom ; m++ )
		{
			std::complex<double> partial_sum;
			cblas_zdotu_sub (numMom,&GR[0],1,chebMoms+m*numMom,1, &partial_sum);
			sum+=partial_sum*GL[m];
		}
		return 2.0*(I*sum).real()*sin(theta);//the sin(theta) is due to changing the integration from dx -> sin(theta)dtheta
	}

	std::vector< std::complex<double> > GR;
	std::vector< double> GL;
	std::complex<double> *chebMoms;
	const int numMom;

};

int main(int argc, char *argv[])
{	
	if (argc != 3)
	{
		printHelpMessage();
		return 0;
	}
	else
		printWelcomeMessage();

	std::size_t pos = std::string(argv[1]).find(".chebmom2D"); 
	if( pos== std::string::npos){ std::cerr<<"The first argument does not seem to be a valid .chebmom2D file"<<std::endl; return -1;}
	const std::string
	momfilename = argv[1],
	simulation_label = momfilename.substr(0,pos);     // get from "live" to the end
	const double
	broadening = atof(argv[2])/1000.;

	std::ifstream momfile( momfilename.c_str() );
	double HalfWidth,BandCenter; int dim;
	int numMoms0,numMoms1;
	momfile>>dim>>HalfWidth>>BandCenter; HalfWidth=HalfWidth/2.0;//in the file what you have is the bandwidth
	momfile>>numMoms0>>numMoms1;
	const int maxNumMom = ( (numMoms0 > numMoms1) ? numMoms0 : numMoms1 ) ;
 	chebMom2D mu(maxNumMom,maxNumMom); double rmu,imu;
	for( int m0 = 0 ; m0 < numMoms0 ; m0++)
	for( int m1 = 0 ; m1 < numMoms1 ; m1++)
	{
		const double
		phi_J = M_PI/(maxNumMom+1),
		g_D_m0=( (maxNumMom-m0+1)*cos( phi_J*m0)+ sin(phi_J*m0)*cos(phi_J)/sin(phi_J) )/(maxNumMom+1),
		g_D_m1=( (maxNumMom-m1+1)*cos( phi_J*m1)+ sin(phi_J*m1)*cos(phi_J)/sin(phi_J) )/(maxNumMom+1);
		momfile>>rmu>>imu;
		mu(m0,m1) = std::complex( rmu,imu)*g_D_m0*g_D_m1;
	}
	momfile.close();


	double lambda = maxNumMom *broadening/HalfWidth/1000.;

	//By performing the transformation x = Cos(theta)
	//all nodes of the matGamma matrix are equally spaced
	//on the interval (0,pi) and the separation is given by 1/M.
	//Hence, in this domain of integration, one would need at least
	// 10*M energy point to capture the oscilattion properly
	const int num_angles = 10*maxNumMom;
	std::vector< double >  angles(num_angles,0);



	const double
	xbound = 0.9,
    	theta_min = acos( xbound),
    	theta_max = acos(-xbound);
	for( int i=0; i < num_angles; i++)
		angles[i] = theta_min + i*(theta_max-theta_min)/(num_angles-1) ;

	std::string
	outputName  ="KuboBastin_sum_"+simulation_label+".dat";

	std::cout<<"Saving the data in "<<outputName<<std::endl;
	std::ofstream outputfile( outputName.c_str() );
	KuboFunctor kuboFun( maxNumMom ); kuboFun.chebMoms = &mu(0,0);
	double acc=0;
	for( std::vector< double >::iterator it =  angles.begin();
	   				it!=  angles.end();
								    it++)
	{
		acc += dim*kuboFun(*it)/HalfWidth/HalfWidth*(theta_max-theta_min)/(num_angles-1);
		outputfile<<cos(*it)*HalfWidth + BandCenter <<" "<<acc<<std::endl;
//		outputfile<<cos(*it)*HalfWidth + BandCenter <<" "<<dim*kuboFun(*it)/HalfWidth/HalfWidth<<std::endl;
	}
	outputfile.close();


	std::cout<<"The program finished succesfully."<<std::endl;
return 0;
}
	

void printHelpMessage()
{
	std::cout << "The program should be called with the following options: moments_filename broadening(meV)" << std::endl
			  << std::endl;
	std::cout << "moments_filename will be used to look for .chebmom2D file" << std::endl;
	std::cout << "broadening in (meV) will define the broadening of the delta functions" << std::endl;
};

void printWelcomeMessage()
{
	std::cout << "WELCOME: This program will compute the chebyshev sum of the kubo-bastin formula for non equlibrium properties" << std::endl;
};
