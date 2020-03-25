#include "chebyshev_moments.hpp"

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
  momfile>>dbuff; this->TimeStep(dbuff);
  momfile>>dbuff; this->TimeCoeff(dbuff);

  //create the moment array and read the data
	
  momfile>>this->numMoms>>this->numTimes;

  this->MomentVector( Moments::vector_t(numMoms * numTimes, 0.0) );
  double rmu, imu;
  for( int m = 0; m < numMoms; m++)
    for( int n = 0; n < numTimes; n++)
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
	     << this->TimeStep() << " " << this->TimeCoeff() << " " << std::endl;
  //Print the number of moments for all directions in a line
  outputfile << numMoms << " " << numTimes << " " << std::endl;

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
	   <<this->HighestMomentNumber()<< " x " <<this->HighestTimeNumber()<<")"<<std::endl;
  std::cout<<"\tSCALE FACTOR:\t\t"<<this->ScaleFactor()<<std::endl;
  std::cout<<"\tSHIFT FACTOR:\t\t"<<this->ShiftFactor()<<std::endl;
  std::cout<<"\tENERGY SPECTRUM:\t("
	   <<-this->HalfWidth()+this->BandCenter()<<" , "
	   << this->HalfWidth()+this->BandCenter()<<")"<<std::endl<<std::endl;
  std::cout<<"\tTIME STEP:\t\t"<<this->TimeStep()<<std::endl;
  std::cout<<"\tTIME COEFFICIENT:\t"<<this->TimeCoeff()<<std::endl;
  
};

void chebyshev::MomentsTD::MomentNumber(const size_t mom)
{ 
  assert(mom <= numMoms);
  chebyshev::MomentsTD new_mom( mom, this->HighestTimeNumber() );
  for( size_t m = 0 ; m < mom ; m++)
    for( size_t n = 0 ; n < numTimes ; n++)
      new_mom(m, n) = this->operator()(m, n); 
  this->numMoms = new_mom.MomentNumber();
  this->MomentVector( new_mom.MomentVector() );
};

void chebyshev::MomentsTD::ApplyJacksonKernel( const double broad )
{
  assert( broad > 0);
  const double eta = 2.0*broad/1000/this->BandWidth();
		
  int maxMom =  ceil(M_PI/eta);
  
  if(  maxMom > numMoms ) maxMom = numMoms;
  std::cout << "Kernel reduced the number of moments to " << maxMom << std::endl;
  this->MomentNumber(maxMom);

  const double phi_J = M_PI/(double)(numMoms+1.0);
  double g_D_m;

  for( size_t m = 0 ; m < numMoms ; m++)
    {
      g_D_m = ( (numMoms - m + 1) * cos(phi_J * m) + sin(phi_J * m) * cos(phi_J) / sin(phi_J) ) * phi_J/M_PI;
      for( size_t n = 0 ; n < numTimes ; n++) this->operator()(m, n) *= g_D_m;
    }
}
