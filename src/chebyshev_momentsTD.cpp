#include "chebyshev_moments.hpp"


int chebyshev::MomentsTD::Evolve( vector_t& Phi)
{
	const auto dim = this->SystemSize();

	if( this->Chebyshev0().size()!= dim )
		this->Chebyshev0() = Moments::vector_t(dim,Moments::value_t(0)); 

	if( this->Chebyshev1().size()!= dim )
		this->Chebyshev1() = Moments::vector_t(dim,Moments::value_t(0)); 
	//From now on this-> will be discarded in Chebyshev0() and Chebyshev1()

	const auto I = Moments::value_t(0, 1);
	const auto x = this->TimeDiff()*this->ChebyshevFreq();
	const double tol = 1e-10;
	double momcutOff = 2.0*tol;
	
		//Initial block
	linalg::copy(Phi , Chebyshev0() );
	this->Hamiltonian().Multiply(Chebyshev0(), Chebyshev1());

	int n = 0;
	auto nIp =Moments::value_t(1, 0);
	double  Jn = 0.5*chebyshev::besselJ(n,x);
	linalg::scal(0,Phi); //Set to zero
	while( momcutOff > tol)
	{
		//---Save J0 to Phi. n = 0,1,2,3,..,
		linalg::axpy( nIp*2*Jn , Chebyshev0(), Phi);

		// Evolve n to n+1
		this->Hamiltonian().Multiply(2.0,  Chebyshev1(), -1.0,  Chebyshev0() );
		Chebyshev0().swap(Chebyshev1());

		nIp*=-I ;
		n++;
		Jn  = chebyshev::besselJ(n,x);
		momcutOff = std::fabs(Jn) ;
	}
	const auto exp0= exp(-I*this->ChebyshevFreq_0()*this->TimeDiff() );
	linalg::scal(exp0,Phi); //Set to zero

  return 0;
};


chebyshev::MomentsTD::MomentsTD( std::string momfilename)
{
  //Check if the input_momfile have the right extension 
  std::size_t ext_pos = string( momfilename ).find(".chebmomTD"); 
  if( ext_pos == string::npos )
    { std::cerr<<"The first argument does not seem to be a valid .chebmomTD file"<<std::endl; assert(false);}

  //if it does, use it to get the extension
  this->SystemLabel( momfilename.substr(0,ext_pos) ); 

  //and then try to open the file
  std::ifstream momfile( momfilename.c_str() );
  assert( momfile.is_open() );

  //if succesful, read the header
  int ibuff; double dbuff;
  momfile>>ibuff; this->SystemSize(ibuff);
  momfile>>dbuff; this->BandWidth(dbuff);
  momfile>>dbuff; this->BandCenter(dbuff);
  momfile>>ibuff; this->MaxTimeStep(ibuff);
  momfile>>dbuff; this->TimeDiff(dbuff);

  //create the moment array and read the data
	
  momfile>>this->_numMoms>>this->_maxTimeStep;

  this->MomentVector( Moments::vector_t( this->HighestMomentNumber() * this->MaxTimeStep(), 0.0) );
  double rmu, imu;
  for( int m = 0; m < this->HighestMomentNumber() ; m++)
  for( int n = 0; n < this->MaxTimeStep() ; n++)
  { 
	  momfile>>rmu>>imu;
	  this->operator()(m,n) = Moments::value_t(rmu, imu);
  }
  momfile.close();
};

void chebyshev::MomentsTD::saveIn(std::string filename)
{
  typedef std::numeric_limits<double> dbl;
  ofstream outputfile(filename.c_str());
  outputfile.precision(dbl::digits10);
  outputfile << this->SystemSize() << " " << this->BandWidth() << " " << this->BandCenter() << " " 
	     << this->MaxTimeStep() << " " << this->TimeDiff() << " " << std::endl;
  //Print the number of moments for all directions in a line
  outputfile << this->HighestMomentNumber() << " " << this->MaxTimeStep() << " " << std::endl;

  for ( auto mom : this->MomentVector() )
    outputfile << mom.real() << " " << mom.imag() << std::endl;
  outputfile.close();
};

void chebyshev::MomentsTD::Print()
{
  std::cout<<"\n\nCHEBYSHEV TD MOMENTS INFO"<<std::endl;
  std::cout<<"\tSYSTEM:\t\t\t"<<this->SystemLabel()<<std::endl;
  if( this-> SystemSize() > 0 )
    std::cout<<"\tSIZE:\t\t\t"<<this-> SystemSize()<<std::endl;

  std::cout<<"\tMOMENTS SIZE:\t\t"<<"("
	   <<this->HighestMomentNumber()<< " x " <<this->MaxTimeStep()<<")"<<std::endl;
  std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
  std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
  std::cout<<"\tENERGY SPECTRUM:\t("
	   <<-this->HalfWidth()+this->BandCenter()<<" , "
	   << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;
  std::cout<<"\tTIME STEP:\t\t"<<this->MaxTimeStep()<<std::endl;
  std::cout<<"\tTIME DIFF:\t"<<this->TimeDiff()<<std::endl;
  
};

void chebyshev::MomentsTD::MomentNumber(const size_t numMoms )
{ 
	assert( numMoms <= this->HighestMomentNumber() );
	const auto maxtime = this->MaxTimeStep();
	
	//Copy all previous moments smaller or equal than the new size into new vector;
	chebyshev::MomentsTD new_mom( numMoms, maxtime);
	for( size_t m = 0 ; m < numMoms  ; m++)
	for( size_t n = 0 ; n < maxtime ; n++)
      new_mom(m, n) = this->operator()(m, n); 

	this->_numMoms = new_mom._numMoms;
	this->MomentVector( new_mom.MomentVector() );
};

void chebyshev::MomentsTD::ApplyJacksonKernel( const double broad )
{
  assert( broad > 0);
  const double eta = 2.0*broad/1000/this->BandWidth();
  int maxMom =  ceil(M_PI/eta);
  
  if(  maxMom > this->HighestMomentNumber() )
	maxMom = this->HighestMomentNumber();
  std::cout << "Kernel reduced the number of moments to " << maxMom <<" for a broadening of "<<M_PI/maxMom << std::endl;
  this->MomentNumber(maxMom);

  const double phi_J = M_PI/(double)(maxMom+1.0);
  double g_D_m;

  for( size_t m = 0 ; m < maxMom ; m++)
  {
	  g_D_m = ( (maxMom - m + 1) * cos(phi_J * m) + sin(phi_J * m) * cos(phi_J) / sin(phi_J) ) * phi_J/M_PI;
	  for( size_t n = 0 ; n < this->MaxTimeStep() ; n++) this->operator()(m, n) *= g_D_m;
  }
}
