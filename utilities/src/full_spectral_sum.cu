
#include "full_spectral_sum.hpp"


void 
utility::kernel::JacksonFunction(const my::real m , const my::real M, my::scalar& mu)
{
	const my::real theta= M_PI/(M+1.);
	mu=mu*my::real( 
		  ( 
			(M-m+1.)*
			cos(theta*m) + 
			sin(theta*m)*cos(theta)/sin(theta) 
		  )/(M+1.) );
}; 

void
utility::kernel::LorentzFunction(const my::real,const my::real lambda, 
					const my::real m , 
					const my::real M, my::scalar& mu)
	{
	  mu=mu*my::real( sinh( lambda*(1.-m/M)) /sinh(lambda) );
	}; // end plus




void 
utility::sum::SpectralConductivity(const std::string moment_filename , my::integer M0,  my::integer M1,
						const my::integer NE,
						const std::string output_filename)
{
	if(M1>M0)
	{
		std::cout<<"Warning M1 cannot be larger than M0"<<std::endl;
	}
	std::ifstream moment_file( moment_filename.c_str() );
	std::ofstream output_file( output_filename.c_str() );

	std::string  label_out;
	//Read M from file;
	my::integer M;
	moment_file>>label_out>>M;
	if(M0==0)
		M0=M;
	if(M1==0)
		M1=M;
		
	//Read the bounds and cutoff from file 
	my::real MaxBound, MinBound,  CutOff,norb,spin,area;
	moment_file>>MinBound>>MaxBound>>CutOff>>norb>>spin>>area;
	std::cout<<"The moments will be computed using : M=("<<M0<<","<<M1<<") "<<MinBound<<" "<<MaxBound<<" "<<CutOff<<" "<<norb<<" "<<spin<<" "<<area<<std::endl;
	const my::real Emin=MinBound;
	const my::real Emax=MaxBound;

	//Read the moments from file
	my::real mure, muim;
	my::scalar mutmp;
	my::integer m0,m1;
	thrust::host_vector< thrust_complex > h_mu(M0*M1);

	for(my::integer m=0;m<M;m++)
	for(my::integer n=0;n<M;n++)
	{		
		moment_file>>m0>>m1>>mure>>muim;
		if( m0<M0 && m1< M1)
		{
			mutmp= my::scalar(mure,muim);
			mutmp=mutmp*my::real(4.0);
			if(m0==0)
				mutmp=mutmp*my::real(0.5);
			if(m1==0)
				mutmp=mutmp*my::real(0.5);
			kernel::JacksonFunction(m0,M0,mutmp);
			kernel::JacksonFunction(m1,M1,mutmp);
			h_mu[m0*M0+m1]= thrust_complex( mutmp );
		}
	}
	
	moment_file.close();

	thrust::device_vector< thrust_complex > d_mu( h_mu );
	
	//Creates operatirs  that will be used by the binary_op
	thrust::counting_iterator<my::integer> index(0);						//Definimos el inidice del vector como iterador
	thrust_complex zero= 0.0;	
	thrust::plus<thrust_complex>	binary_sum;
	//The custom binary operator
	SpectralConductivity_binary cond_bin(M0,M1);							//Definimos la operacion binaria del kernel (mu,Index)

//	std::cout<<"Multiply the result for unknown 10 factor"<<std::endl;
	const my::real scal0 = 4*pow(2.*CutOff/(MaxBound-MinBound),1.)/area/spin;	
	const my::real dE =(Emax-Emin)/NE;
	for(my::real E=Emin;E<=Emax;E=E+dE)
	{
		const my::real En   = ( E-(MaxBound+MinBound)*0.5 )*2.*CutOff/(MaxBound-MinBound) ;
		const my::real scal = scal0/pow( 1 - En*En ,2 );

		cond_bin.En=En;
		const my::scalar 
		SIGMA=thrust::inner_product(d_mu.begin(),d_mu.end(), index, zero,binary_sum,cond_bin)*scal;	
		output_file<< E<<" "<<SIGMA.real()<<" "<<SIGMA.imag()<<std::endl;
	} 
	output_file.close();
};
