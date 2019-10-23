

#include "kpm.hpp"
#include "timing.hpp"



bool IsNumber(const std::string& s) {
	std::string::const_iterator it = s.begin();
	while (it != s.end() && std::isdigit(*it))
		++it;
	return !s.empty() && it == s.end();
}
;


NumCal::integer
GetMomentsPerNode( const NumCal::integer numMoms)
{
	return std::ceil( (NumCal::real) numMoms / (NumCal::real) KPM_MPI_GETPROC() );
}

namespace NumCal
{
NumCal::real
Kpm::EnormToEn(const NumCal::real Enorm) const
{
	return EnergyFactor()*Enorm+0.5*(MaxBound()+MinBound());
};

NumCal::real
Kpm::EnToEnorm(const NumCal::real En) const
{
	return (En-0.5*(MaxBound()+MinBound()) )/EnergyFactor();
};


/**************************SETTERS*************************************/
void
Kpm::SetTruncOrder(const NumCal::size_t _trunc_order)
{
	if( _trunc_order <= 0 )
	{
		std::cerr<<" The trucation order :"<<_trunc_order;
		std::cerr<<" should be higher than zero"<<std::endl;
		std::exit(-1);
	}

	trunc_order_=_trunc_order;
};

void
Kpm::SetBounds(	const NumCal::real _min_bound,
				const NumCal::real _max_bound,
				const NumCal::real _cutoff)
{
	if( _max_bound <= _min_bound)
	{
		std::cerr<<"The bounds:"<<_max_bound<<" "<<_min_bound;
		std::cerr<<" are incorrect"<<std::endl;
		std::exit(-1);
	}
	if( _cutoff >= 1 )
	{
		std::cerr<<"The cutoff:"<<_cutoff;
		std::cerr<<" is incorrect"<<std::endl;
		std::exit(-1);
	}
	min_bound_=_min_bound;
	max_bound_=_max_bound;
	cutoff_=_cutoff;
	scale_factor_= (MaxBound()-MinBound())/( 2.*CutOff() );
};

void
Kpm::DOSMoments(Lattice _lattice, std::string _dosMomFile)
{
	std::cout<<"Warning: The density if states is only calculated in node 0"<<std::endl;
	if( KPM_MPI_GETRANK() == 0)
	{
			//Get the KPM-OMP STATUS
			kpm_util::GetOMPStatus(  );
			_dosMomFile+= kpm_util::GetNodeLabel();
			std::ofstream dos_output(_dosMomFile.c_str() );
			//The define the Set0 as m
			const NumCal::integer mSet=0;
			//Initialize the Vector Jm0 with a random vector
			for(NumCal::integer i=0;i< VecDim() ; i++)
				chebSet.ChebVec(mSet,0)[i]= distribution(generator);
			//Normalize the Vector
			linalg.normalize( chebSet.ChebVec(mSet,0) );
			//Initialize the vector Jm1
			linalg.copy(chebSet.ChebVec(mSet,0), chebSet.ChebVec(mSet,1) );
			//Iterate to calculate the density of states
			for(NumCal::integer m=0;m<TruncOrder();m++)
			{
				NumCal::scalar chebMom;
				chebSet.ChebyshevIteration(_lattice , mSet, m);
				chebMom=linalg.dot(chebSet.ChebVec(mSet,0), chebSet.ChebVec(mSet,1) );
				dos_output<<m<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			}
	dos_output.close();
	}
}


void
Kpm::SpinConductivityMoments(Lattice _lattice, std::string _condMomXXFile,std::string _condMomXYFile)
{
	NumCal::random::generator random_generator;

	NumCal::scalar *Chebt;
	///Get the id of the process
	const NumCal::integer id = KPM_MPI_GETRANK() ;
	///Get the total number of processes
	const NumCal::integer nproc = KPM_MPI_GETPROC();

	///Distribute the moments between the processes
	const NumCal::integer dMom = kpm_util::GetMomentsPerNode( TruncOrder() );
	const NumCal::integer myMom= id * dMom;

	SetTruncOrder (dMom * nproc);
	//Get the KPM-MPI STATUS
	kpm_util::GetMPIStatus( TruncOrder() );
	//Get the KPM-OMP STATUS
	kpm_util::GetOMPStatus(  );

	_condMomXXFile+= kpm_util::GetNodeLabel();
	_condMomXYFile+= kpm_util::GetNodeLabel();

	std::ofstream outputxx(_condMomXXFile.c_str() );
	std::ofstream outputxy(_condMomXYFile.c_str() );

	if(id==0)
	{
		outputxx<<_lattice.Label()+"conductivityxx"<<" "
				<<TruncOrder()<<" "<<MinBound()<<" "<<MaxBound()<<" "
				<<CutOff()<<" "<<_lattice.OrbitalNumber()<<" "
				<<_lattice.SpinNumber()<<" "<<AREA<<std::endl;

		outputxy<<_lattice.Label()+"conductivityxy"<<" "
				<<TruncOrder()<<" "<<MinBound()<<" "<<MaxBound()<<" "
				<<CutOff()<<" "<<_lattice.OrbitalNumber()<<" "
				<<_lattice.SpinNumber()<<" "<<AREA<<std::endl;
	}


	//The define the Set0 as m and the Set 1 as n
	const NumCal::integer mSet=0,nSet=1;

	//Initialize the Vector Jm0 with a random vector
	for(NumCal::integer i=0;i< VecDim() ; i++)
		chebSet.ChebVec(mSet,0)[chebSet.MemSep()*i]= distribution(generator);

	//Normalize the Vector
	linalg.normalize( chebSet.ChebVec(mSet,0) );
	//The initial moment is then copy into the initial
	//moment of the n-set vector to iterate
	linalg.copy(chebSet.ChebVec(mSet,0),
				chebSet.ChebVec(nSet,0)
				);
	//Multiply the initial vector by the velocity operator
	///Needs to be fixed for chebyshev sets
	_lattice.ApplyVelocity(chebSet.MemSep(), 0, chebSet.ChebVec(mSet,0), chebSet.ChebVec(mSet,1) );

	///Iterate the vector until  the moment corresponding to the process
	for(NumCal::integer m=0;m< myMom ;m++)
		chebSet.ChebyshevIteration( _lattice , mSet, m);
	//for optimization purposes we will use the Jm0 vector
	// as a buffer for this calculation
	Chebt=chebSet.ChebVec(mSet,0);

	//When reaching the correct initial moment for the
	//m set, proceed to initialize the n ster
	if(id == 0)
		cycletime(-1);
	for(NumCal::integer m=myMom;m<myMom+dMom;m++)
	{
		NumCal::scalar chebMom;
		chebSet.ChebyshevIteration( _lattice , mSet, m);
		//The initial moment is copy into the element Jn1
		linalg.copy(chebSet.ChebVec(nSet,0), chebSet.ChebVec(nSet,1) );
		for(NumCal::integer n=0;n< TruncOrder();n++)
		{
			if(id == 0)
				cycletime(dMom*TruncOrder());

			chebSet.ChebyshevIteration( _lattice , nSet, n);
			_lattice.ApplyVelocity( chebSet.MemSep(), 0,chebSet.ChebVec(nSet,1), Chebt );

			chebMom= linalg.dot( chebSet.ChebVec(mSet,1), Chebt );
			outputxx<<m<<" "<<n<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			//for optimization purposes we will use the Jm0 vector
			// as a buffer for this calculation
			_lattice.ApplySpinZ_Velocity( chebSet.MemSep(),1, chebSet.ChebVec(nSet,1), Chebt );
			chebMom= linalg.dot( chebSet.ChebVec(mSet,1), Chebt );
			outputxy<<m<<" "<<n<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			if(id == 0)
				cycletime(dMom*TruncOrder());

		}
	}
	outputxx.close();
	outputxy.close();
};


void
Kpm::ConductivityMoments(Lattice _lattice, std::string _condMomXXFile,std::string _condMomXYFile)
{
	NumCal::random::generator random_generator;

	NumCal::scalar *Chebt;
	///Get the id of the process
	const NumCal::integer id = KPM_MPI_GETRANK() ;
	///Get the total number of processes
	const NumCal::integer nproc = KPM_MPI_GETPROC();

	///Distribute the moments between the processes
	const NumCal::integer dMom = kpm_util::GetMomentsPerNode( TruncOrder() );
	const NumCal::integer myMom= id * dMom;

	SetTruncOrder (dMom * nproc);
	//Get the KPM-MPI STATUS
	kpm_util::GetMPIStatus( TruncOrder() );
	//Get the KPM-OMP STATUS
	kpm_util::GetOMPStatus(  );

	_condMomXXFile+= kpm_util::GetNodeLabel();
	_condMomXYFile+= kpm_util::GetNodeLabel();

	std::ofstream outputxx(_condMomXXFile.c_str() );
	std::ofstream outputxy(_condMomXYFile.c_str() );

	if(id==0)
	{
		outputxx<<_lattice.Label()+"conductivityxx"<<" "
				<<TruncOrder()<<" "<<MinBound()<<" "<<MaxBound()<<" "
				<<CutOff()<<" "<<_lattice.OrbitalNumber()<<" "
				<<_lattice.SpinNumber()<<" "<<AREA<<std::endl;

		outputxy<<_lattice.Label()+"conductivityxx"<<" "
				<<TruncOrder()<<" "<<MinBound()<<" "<<MaxBound()<<" "
				<<CutOff()<<" "<<_lattice.OrbitalNumber()<<" "
				<<_lattice.SpinNumber()<<" "<<AREA<<std::endl;
	}


	//The define the Set0 as m and the Set 1 as n
	const NumCal::integer mSet=0,nSet=1;

	//Initialize the Vector Jm0 with a random vector
	for(NumCal::integer i=0;i< VecDim() ; i++)
		chebSet.ChebVec(mSet,0)[chebSet.MemSep()*i]= distribution(generator);

	//Normalize the Vector
	linalg.normalize( chebSet.ChebVec(mSet,0) );
	//The initial moment is then copy into the initial
	//moment of the n-set vector to iterate
	linalg.copy(chebSet.ChebVec(mSet,0),
				chebSet.ChebVec(nSet,0)
				);
	//Multiply the initial vector by the velocity operator
	///Needs to be fixed for chebyshev sets
	_lattice.ApplyVelocity(chebSet.MemSep(), 0, chebSet.ChebVec(mSet,0), chebSet.ChebVec(mSet,1) );

	///Iterate the vector until  the moment corresponding to the process
	for(NumCal::integer m=0;m< myMom ;m++)
		chebSet.ChebyshevIteration( _lattice , mSet, m);
	//for optimization purposes we will use the Jm0 vector
	// as a buffer for this calculation
	Chebt=chebSet.ChebVec(mSet,0);

	//When reaching the correct initial moment for the
	//m set, proceed to initialize the n ster
	if(id == 0)
		cycletime(-1);
	for(NumCal::integer m=myMom;m<myMom+dMom;m++)
	{
		NumCal::scalar chebMom;
		chebSet.ChebyshevIteration( _lattice , mSet, m);
		//The initial moment is copy into the element Jn1
		linalg.copy(chebSet.ChebVec(nSet,0), chebSet.ChebVec(nSet,1) );
		for(NumCal::integer n=0;n< TruncOrder();n++)
		{
			if(id == 0)
				cycletime(dMom*TruncOrder());

			chebSet.ChebyshevIteration( _lattice , nSet, n);
			_lattice.ApplyVelocity( chebSet.MemSep(), 0,chebSet.ChebVec(nSet,1), Chebt );

			chebMom= linalg.dot( chebSet.ChebVec(mSet,1), Chebt );
			outputxx<<m<<" "<<n<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			//for optimization purposes we will use the Jm0 vector
			// as a buffer for this calculation
			_lattice.ApplyVelocity( chebSet.MemSep(),1, chebSet.ChebVec(nSet,1), Chebt );
			chebMom= linalg.dot( chebSet.ChebVec(mSet,1), Chebt );
			outputxy<<m<<" "<<n<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			if(id == 0)
				cycletime(dMom*TruncOrder());

		}
	}
	outputxx.close();
	outputxy.close();
};


void
Kpm::ValleyConductivityMoments(NumCal::integer vindex,NumCal::real RKscal,Lattice _lattice, std::string _condMomXXFile,std::string
_condMomXYFile)
{
	NumCal::scalar *Chebt;
	///Get the id of the process
	const NumCal::integer id = KPM_MPI_GETRANK() ;
	///Get the total number of processes
	const NumCal::integer nproc = KPM_MPI_GETPROC();

	///Distribute the moments between the processes
	const NumCal::integer dMom = kpm_util::GetMomentsPerNode( TruncOrder() );
	const NumCal::integer myMom= id * dMom;

	SetTruncOrder (dMom * nproc);
	//Get the KPM-MPI STATUS
	kpm_util::GetMPIStatus( TruncOrder() );
	//Get the KPM-OMP STATUS
	kpm_util::GetOMPStatus(  );

	_condMomXXFile+= kpm_util::GetNodeLabel();
	_condMomXYFile+= kpm_util::GetNodeLabel();

	std::ofstream outputxx(_condMomXXFile.c_str() );
	std::ofstream outputxy(_condMomXYFile.c_str() );


	if(id==0)
	{
		outputxx<<_lattice.Label()+"valley_conductivityxx"<<" "
				<<TruncOrder()<<" "<<MinBound()<<" "<<MaxBound()<<" "
				<<CutOff()<<" "<<_lattice.OrbitalNumber()<<" "
				<<_lattice.SpinNumber()<<" "<<AREA<<std::endl;

		outputxy<<_lattice.Label()+"valley_conductivityxx"<<" "
				<<TruncOrder()<<" "<<MinBound()<<" "<<MaxBound()<<" "
				<<CutOff()<<" "<<_lattice.OrbitalNumber()<<" "
				<<_lattice.SpinNumber()<<" "<<AREA<<std::endl;
	}

//The define the Set0 as m and the Set 1 as n
	const NumCal::integer mSet=0,nSet=1;

	//Initialize the Vector Jm0 with a random vector
	for(NumCal::integer i=0;i< VecDim() ; i++)
		chebSet.ChebVec(mSet,0)[chebSet.MemSep()*i]= distribution(generator);
	//Normalize the Vector
	linalg.normalize( chebSet.ChebVec(mSet,0) );
	//The initial moment is then copy into the initial
	//moment of the n-set vector to iterate
	linalg.copy(chebSet.ChebVec(mSet,0),
				chebSet.ChebVec(nSet,0)
				);
	//Multiply the initial vector by the velocity operator
	///Needs to be fixed for chebyshev sets
	_lattice.ApplyProjValley_Velocity(chebSet.MemSep(),vindex,RKscal, 0, chebSet.ChebVec(mSet,0), chebSet.ChebVec(mSet,1) );

	///Iterate the vector until  the moment corresponding to the process
	for(NumCal::integer m=0;m< myMom ;m++)
		chebSet.ChebyshevIteration( _lattice , mSet, m);
	//for optimization purposes we will use the Jm0 vector
	// as a buffer for this calculation
	Chebt=chebSet.ChebVec(mSet,0);

	//When reaching the correct initial moment for the
	//m set, proceed to initialize the n ster
	if(id == 0)
		cycletime(-1);

	for(NumCal::integer m=myMom;m<myMom+dMom;m++)
	{
		NumCal::scalar chebMom;
		chebSet.ChebyshevIteration( _lattice , mSet, m);
		//The initial moment is copy into the element Jn1
		linalg.copy(chebSet.ChebVec(nSet,0), chebSet.ChebVec(nSet,1) );
		for(NumCal::integer n=0;n< TruncOrder();n++)
		{
			if(id == 0)
				cycletime(dMom*TruncOrder());

			chebSet.ChebyshevIteration( _lattice , nSet, n);
			_lattice.ApplyVelocity( chebSet.MemSep(), 0,chebSet.ChebVec(nSet,1), Chebt );

			chebMom= linalg.dot( chebSet.ChebVec(mSet,1), Chebt );
			outputxx<<m<<" "<<n<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			//for optimization purposes we will use the Jm0 vector
			// as a buffer for this calculation
			_lattice.ApplyVelocity( chebSet.MemSep(), 1, chebSet.ChebVec(nSet,1), Chebt );

			chebMom= linalg.dot( chebSet.ChebVec(mSet,1), Chebt );
			outputxy<<m<<" "<<n<<" "<<chebMom.real()<<" "<<chebMom.imag()<<std::endl;
			if(id == 0)
				cycletime(dMom*TruncOrder());

		}
	}
	outputxx.close();
	outputxy.close();
};

};
