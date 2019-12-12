#ifndef CHEBYSHEV_KUBO
#define CHEBYSHEV_KUBO
/*
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
		delta_Functor GLfun(theta); greenR_Functor GRfun(theta);

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
			cdot(numMom,&GR[0],chebMoms+m*numMom,&partial_sum);
			sum+=partial_sum*GL[m];
		}
		return -2.0*(I*sum).real()*sin(theta);//the sin(theta) is due to changing the integration from dx -> sin(theta)dtheta
	}

	std::vector< std::complex<double> > GR;
	std::vector< double> GL;
	std::complex<double> *chebMoms;
	const int numMom;
};
*/


#endif
