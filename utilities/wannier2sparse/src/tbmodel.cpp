
#include "tbmodel.hpp"


void tbmodel::readUnitCell(const string inputfile)
{
	lat_vecs = read_unit_cell_file(inputfile);
}



void tbmodel::readOrbitalPositions(const string inputfile)
{
	orbPos_list = read_xyz_file(inputfile);
};



void tbmodel::readWannierModel(const string inputfile)
{
	hl.add_from_wannier(read_wannier_file(inputfile));
};



void tbmodel::readStaticDisorder(const string inputfile)
{
	auto dis_list = read_disorder_file(inputfile);
	if (dis_list.size() == 0 )
	{
		std::cout<<"Pristine calculation"<<std::endl;
		return ;
	}
	
	hl.add_random_hopping_list( dis_list);
};

hopping_list 
tbmodel::createLocalHopping_list( )
{
	assert(volume(lat_vecs) > 0 );
	assert(orbPos_list.size()==hl.WannierBasisSize());
		
	hopping_list chl ;
	chl.SetWannierBasisSize(hl.WannierBasisSize()); //increase the supercell dimension;
	
	for( auto& elem: this->hl.hoppings )
	{
		auto  cellID = get<0>(elem.second);
		auto  value  = get<1>(elem.second);
		auto  edge   = get<2>(elem.second);
		
		bool is_onsite = true;
		for( auto x : (array<int,3>)cellID )
			is_onsite*=(x==0);
	
		if( is_onsite )
			chl.AddHopping( cellID,value,edge) ;
	}

	return chl;
};

hopping_list 
tbmodel::createHoppingDensity_list( oputil::op_matrix OP )
{
	assert(volume(lat_vecs) > 0 );
	assert(orbPos_list.size()==hl.WannierBasisSize());
	assert(OP.size()==hl.WannierBasisSize());
		
	
	hopping_list chl ;
	chl.SetWannierBasisSize(hl.WannierBasisSize()); //increase the supercell dimension;
	
	int i_orb = 0;
	hopping_list::cellID_t cellID={0,0,0};//onsite
	for( auto orb_id_0 : this->orbPos_list )
	{
		int j_orb = 0;
		for( auto orb_id_1 : this->orbPos_list )
		{
			hopping_list::edge_t vertex_edge = {i_orb,j_orb};
			hopping_kind::value_t value = OP(i_orb,j_orb);
			chl.AddHopping( cellID,value,vertex_edge) ;
			j_orb++;
		}
		i_orb++;
	}
	return chl;
};

hopping_list tbmodel::createHoppingCurrents_list(const int dir)
{
	assert( dir <3 && dir >=0 );
	assert(volume(lat_vecs) > 0 );  const auto av_length = pow( volume(lat_vecs), 1.0/3.0 );
	assert_equal( (int)orbPos_list.size(), hl.WannierBasisSize());

	auto curr = this->hl ;
	for( auto& elem: curr.hoppings )
	{
		auto  cellID = get<0>(elem.second);
		auto& value  = get<1>(elem.second);
		auto  edge   = get<2>(elem.second);
		
		//Extract the final orbital position from orbPos_list, which is a tuple
		auto orb_diff= get<1>( orbPos_list[edge[1]] ) ;
		
		//Substract the initial position position from orbPos_list, which is a tuple
		for( int i=0; i < orb_diff.size(); i++)
			orb_diff[i] +=-get<1>(orbPos_list[edge[0]])[i]; 
			
		//Get the displacement vector in cartesian
		auto displ_vec = tag2cartesian(cellID, lat_vecs); 

		//Add the displacement vectors in cartesians to the orbitals
		for( int i=0; i < orb_diff.size(); i++)
			orb_diff[i] += displ_vec[i];    

		//Change the hopping element accordingly
		if( fabs(orb_diff[dir])/av_length > 1E-10 ) 
			value.Rescale( hopping_kind::value_t( 0.0, -orb_diff[dir] ) );
		else
			value.Rescale( hopping_kind::value_t( 0.0, 0.0 ) );
	}
	return curr;
};


// hopping_list 
// tbmodel::createHoppingSpinorialDensity_list(const std::array< std::complex<double>,4 > op)
// {
	//CREATE A SPINORIAL OPERATOR MATRIX
// 	auto sop = this->createSpinorialOp(op);
// 	return this->createHoppingDensity_list(sop);
// };


hopping_list 
tbmodel::createHoppingSpinDensity_list(const double theta, const double phi)
{
	typedef std::complex<double> mycomplex;
	mycomplex I(0,1);
	//CREATE A SPIN MATRIX
	std::array< mycomplex,4 >
	op = { 
			cos(theta)			   ,sin(theta)*exp(-I*phi ),
			sin(theta)*exp( I*phi ),-cos(theta)
		 }; 
	auto sop = this->createSpinorialOp(op);
 	return this->createHoppingDensity_list(sop);
};

hopping_list 
tbmodel::createHoppingSpinProjectionDensity_list(const double theta, const double phi, const int s)
{
    //CREATE A SPIN PROJECTOR MATRIX
    const std::complex<double> I(0,1);
	std::array< std::complex<double>,4 > 
	op = { 
	  0.5 * (1 + s * cos(theta))          ,0.5 * s * sin(theta) * exp(-I * phi),
	  0.5 * s * sin(theta) * exp( I * phi),0.5 * (1 - s * cos(theta))
	};
	auto sop = this->createSpinorialOp(op);
    return this->createHoppingDensity_list(sop);
};

hopping_list 
tbmodel::createHoppingTorqueDensity_list(const double theta, const double phi)
{
	//CREATE A TORQUE MATRIX
	auto Tu = this->createTorqueMatrix(theta,phi );
	return this->createHoppingDensity_list( Tu );
};


hopping_list 
tbmodel::createHoppingSpinCurrents_list(const int dir, const double theta, const double phi )
{
	hopping_kind::value_t I(0.0,1.0);
	
	hopping_list scurr;
	scurr.SetWannierBasisSize(hl.WannierBasisSize()); //increase the supercell dimension;
	scurr.SetBounds(hl.Bounds() ); //increase the supercell dimension;
	
	for( const auto elem: this->createHoppingCurrents_list(dir).hoppings )
	{
		const auto _cellID	= get<0>(elem.second);
		const auto _value	= get<1>(elem.second);
		const auto _edge	= get<2>(elem.second);
		const int  is0		= this->getSpinIndex(_edge[0]);
		const int  is1		= this->getSpinIndex(_edge[1]);
		const auto sz  		= 1.0-2.0*is0;

		if ( is0 == is1 ) //For upup or dndn
		{
			const auto sval_uu = sz*cos(theta);
			const auto sval_ud = 0.5*sin(theta)*exp(-I*phi);

			if( std::norm(sval_uu)>1e-5)
			{
				auto value_uu = _value;  value_uu.Rescale(sval_uu);
				scurr.AddHopping(_cellID,value_uu,ChangeSpinIndex(_edge,{is0,is0}));
			}

			if( std::norm(sval_ud)>1e-5)
			{
				auto value_ud = _value; value_ud.Rescale(sval_ud);
				auto value_du = _value; value_du.Rescale(std::conj(sval_ud) );
				scurr.AddHopping(_cellID,value_ud,ChangeSpinIndex(_edge,{SPIN_UP,SPIN_DW}));
				scurr.AddHopping(_cellID,value_du,ChangeSpinIndex(_edge,{SPIN_DW,SPIN_UP}));
			}
		}
		else // if updw sz ==1  if dwup sz ==-1
		{	
			const auto sval_uu = 0.5*sin(theta)*exp( sz*I*phi);
			if( std::norm(sval_uu)>1e-5)
			{
				auto value_uu = _value; value_uu.Rescale(sval_uu);
				scurr.AddHopping(_cellID,value_uu,ChangeSpinIndex(_edge,{SPIN_UP,SPIN_UP}));
				scurr.AddHopping(_cellID,value_uu,ChangeSpinIndex(_edge,{SPIN_DW,SPIN_DW}));
			}
		}
	}

  return scurr;
};


hopping_list 
tbmodel::createHoppingSpinDensity_list(const std::string op )
{
	char direction = oputil::spin_direction(op);
	if( direction == 'X' )
		return this->tbmodel::createHoppingSpinDensity_list(0.5*M_PI,0.0*M_PI);
	if( direction == 'Y' )
		return this->tbmodel::createHoppingSpinDensity_list(0.5*M_PI,0.5*M_PI);
	if( direction == 'Z' )
		return this->tbmodel::createHoppingSpinDensity_list(0.0*M_PI,0.0*M_PI);

	std::cout<<"Incorrect spin requested in createHoppingSpinDensity_list: "<<direction<<std::endl;
	assert(false);
};


hopping_list 
tbmodel::createHoppingSpinProjectionDensity_list(const std::string op)
{
	char direction = oputil::projection_direction(op);
	int s = +1; //Assume is UP.
	switch(direction)
	{
		case 'U':
			s = +1;
		break;

		case 'D':
			s = -1;
		default:
			s = +1;
	}
	direction = oputil::spin_direction(op);
	return  this->createHoppingSpinProjectionDensity_list(s, direction);
}


hopping_list 
tbmodel::createHoppingSpinProjectionDensity_list(const int s, const char sdir)
{
	switch(sdir)
	{
		case 'X':
		  return createHoppingSpinProjectionDensity_list(0.5*M_PI, 0.0*M_PI, s);
		break;

		case 'Y':
		  return createHoppingSpinProjectionDensity_list(0.5*M_PI,0.5*M_PI, s);
		break;

		case 'Z':
		  return createHoppingSpinProjectionDensity_list(0.0*M_PI,0.0*M_PI, s);
		break;
	}
	std::cout<<"Incorrect spin requested in createHoppingSpinProjectionDensity_list: s"<<s<<"S"<<sdir<<std::endl;
	assert(false);	
};

hopping_list 
tbmodel::createHoppingTorqueDensity_list(const std::string op )
{
	char direction = oputil::spin_direction(op);
	if( direction == 'X' )
		return this->tbmodel::createHoppingTorqueDensity_list(0.5*M_PI,0.0*M_PI);
	if( direction == 'Y' )
		return this->tbmodel::createHoppingTorqueDensity_list(0.5*M_PI,0.5*M_PI);
	if( direction == 'Z' )
		return this->tbmodel::createHoppingTorqueDensity_list(0.0*M_PI,0.0*M_PI);

	std::cout<<"Incorrect spin requested in createHoppingSpinDensity_list: "<<direction<<std::endl;
	assert(false);
};


hopping_list 
tbmodel::createHoppingCurrents_list(const std::string op )
{
	char direction = oputil::velocity_direction(op);
	int dir = 0; //Assume is in X.
	if( direction == 'Y' )
		dir = 1;
	if( direction == 'Z' )
		dir = 2;
	return  this->createHoppingCurrents_list(dir);
};


hopping_list 
tbmodel::createHoppingSpinCurrents_list(string op)
{
	char direction = oputil::velocity_direction(op);
	int dir = 0; //Assume is in X.
	switch(direction)
	{
		case 'X':
			dir = 0;
		break;

		case 'Y':
			dir = 1;
		break;

		case 'Z':
			dir = 2;
		break;

		default:
			dir = -1;
	}
	direction = oputil::spin_direction(op);
	return  this->createHoppingSpinCurrents_list(dir,direction );
}


hopping_list 
tbmodel::createHoppingSpinCurrents_list(const int dir, const char sdir)
{
	std::cout<<"Which mean dir="<<dir<<std::endl;
	switch(sdir)
	{
		case 'X':
			return createHoppingSpinCurrents_list(dir,0.5*M_PI,0.0*M_PI );
		break;

		case 'Y':
			return createHoppingSpinCurrents_list(dir,0.5*M_PI,0.5*M_PI );
		break;

		case 'Z':
			return createHoppingSpinCurrents_list(dir,0.0*M_PI,0.0*M_PI );
		break;
	}
	std::cout<<"Incorrect spin requested in createHoppingSpinCurrents_list: V"<<dir<<"S"<<sdir<<std::endl;
	assert(false);	
};


hopping_list tbmodel::WannierOperator(std::string op_id )
{
	if ( oputil::is_velocity(op_id)  )
	{
		return this->createHoppingCurrents_list(op_id) ;
	}

	if ( oputil::is_spinvelocity(op_id) )
	{
		return this->createHoppingSpinCurrents_list(op_id);
	}

	if ( oputil::is_spin(op_id) )
	{
		return this->createHoppingSpinDensity_list(op_id);
	}

	if ( oputil::is_spinprojection(op_id) )
	{
		return this->createHoppingSpinProjectionDensity_list(op_id);
	}

	if ( oputil::is_torque(op_id) )
	{
		return this->createHoppingTorqueDensity_list(op_id);
	}
	std::ifstream f(op_id+".OP");
	if ( f.good() )
	{
		
		std::cout<<"The file "<<op_id+".OP "<<"was found"<<std::endl;

		std::string type ;
		f>>type;
		
		using namespace std;
		const int num_elems = 4;
		array< complex<double>,num_elems> op_elems;
		
		for( auto& elem : op_elems )
		{
			double re, im;
			f>>re>>im;
			elem = {re, im};
		}
		
		std::cout<<op_id<<"is of the "<<type<<" type"<<std::endl;
		
		auto sop = this->createSpinorialOp(op_elems);
		return this->createHoppingDensity_list(sop);
	}
	
	
	std::cout<<"Incorrect wannier operator requested "<<op_id<<std::endl;
	assert(false);		
};


oputil::op_matrix tbmodel::createSpinorialOp(const std::array< std::complex<double>,4 > op)
{
	const int tot_dim  =  this->orbPos_list.size();
	const int num_spin = 2;
	const int orb_dim  =  tot_dim/num_spin;

	oputil::op_matrix SOP ( tot_dim ) ;
	for(int io = 0 ; io < orb_dim ; io++)
	for(int si = 0 ; si < num_spin; si++)
	for(int sj = 0 ; sj < num_spin; sj++)
	{
		const int i = si*orb_dim + io;
		const int j = sj*orb_dim + io;

		const int opidx = si*num_spin + sj ; //The operator is ordered in row major order
		SOP(i,j) = op[opidx]; 
	}
	return SOP;
}

// oputil::op_matrix tbmodel::createSpinMatrix(const double theta, const double phi)
// {
//	const std::complex<double> I(0,1);
//	std::array< std::complex<double>,4 > 
//	op = { 
//			cos(theta)			   ,sin(theta)*exp(-I*phi ),
//			sin(theta)*exp( I*phi ),-cos(theta)
//		 }; 
//	return tbmodel::createSpinorialOp(op);
//};

// oputil::op_matrix tbmodel::createSpinProjectionMatrix(const double theta, const double phi, const int s)
// {
//	const std::complex<double> I(0,1);
//	std::array< std::complex<double>,4 > 
//	op = { 
//	  0.5 * (1 + s * cos(theta))          ,0.5 * s * sin(theta) * exp(-I * phi),
//	  0.5 * s * sin(theta) * exp( I * phi),0.5 * (1 - s * cos(theta))
//	}; 
//	return tbmodel::createSpinorialOp(op);
//};

oputil::op_matrix tbmodel::createTorqueMatrix(const double theta, const double phi)
{
	const std::complex<double> I(0,1);
	const int tot_dim  =  this->orbPos_list.size();
	const int num_spin = 2;
	const int orb_dim  =  tot_dim/num_spin;

	//Gather all constant onsites into the onsite matrix
	oputil::op_matrix H00 ( tot_dim ) ;
	auto onsites = this->createLocalHopping_list();
	for( const auto elem: onsites.hoppings )
	{
		auto cellID= get<0>(elem.second);
		auto value = get<1>(elem.second);
		auto edge  = get<2>(elem.second);
		
		if( value.is_constant() )	
			H00( get<0>(edge), get<1>(edge) ) += value(); 
	};
	std::cout<<"Creating torque matrix "<<std::endl;
	//Compute the torque matrix
	oputil::op_matrix Tu ( tot_dim ) ;
	double avg_Jex = 0.0;
	for(int io = 0 ; io < orb_dim ; io++)
	{
		auto i = 0*orb_dim + io;
		auto j = 1*orb_dim + io;
		double J[3] = {0.0,0.0,0.0};

		std::cout<<H00(i,i)<<" "<<H00(i,j)<<std::endl<<H00(i,j)<<" "<<H00(j,j)<<std::endl;

		auto J00 = 0.5*( H00(i,i) - H00(j,j) );
		auto J01 = 0.5*( H00(i,j) + std::conj(H00(j,i)) );

		J[0] = J01.real(); J[1] =-J01.imag(); J[2] = J00.real();
		double Jex = sqrt( J[0]*J[0] + J[1]*J[1]+ J[2]*J[2] );
		avg_Jex += Jex; 
		std::cout<<"For orb "<<io+1<<" unit_m  = ("<<-J[0]/Jex<<" "<<-J[1]/Jex<<" "<<-J[2]/Jex<<") and |J|"<<Jex<<" "<<std::endl;

		Tu(i,i) = 1.0*sin(theta)*( J[1]*cos(phi) - J[0]*sin(phi) ); 
		Tu(j,j) =-1.0*sin(theta)*( J[1]*cos(phi) - J[0]*sin(phi) ); 

		Tu(i,j) = (-I*J[0] - J[1])*cos(theta) + J[2]*sin(theta)*( I*cos(phi)+sin(phi) ) ; 
		Tu(j,i) = std::conj( Tu(i,j) ) ; 
	}
	std::cout<<"Average exchange Jex = "<<avg_Jex/orb_dim<<" "<<std::endl;
	return Tu;
};


map<int,string>
tbmodel::map_id2orbs()
{
	map<int,string> id2orbs;
	int id = 0;
	for( auto orb_id : this->orbPos_list )
	{
		auto orb_label = get<0>(orb_id);
		string orbNAME = orb_label ;
		int end_name_pos = orbNAME.find("_");
		orbNAME = orbNAME.substr (0,end_name_pos);
		id2orbs.insert( {id,orbNAME} );
		id++;
	}
	return id2orbs;
}


map<int,int>
tbmodel::map_id2spin()
{
	map<int,int> id2spin;
	int id = 0;
	for( auto orb_id : this->orbPos_list )
	{
		auto orb_label = get<0>(orb_id);
		int sz = 0 ;
		if( orb_label.find("_s+_")!=std::string::npos )
			sz= 1;
		if( orb_label.find("_s-_")!=std::string::npos )
			sz=-1;
		if( sz != 0)
			id2spin.insert( {id,sz} );
		id++;
	}
	return id2spin;
}



array<double,3> tag2cartesian(const array<int,3>& tag,const array < array<double,3> , 3 >& lat_vecs  )
{
	array<double,3> cart_vect = {0,0,0};
	for(int i =0 ; i < cart_vect.size(); i++ )
	for(int j =0 ; j < tag.size(); j++ )
	  cart_vect[i] += tag[j]*lat_vecs[j][i];
	return cart_vect;
};

double volume( const array < array<double,3> , 3 >&  uc )
{
	return  std::fabs(
			( uc[0][1]*uc[1][2]-uc[1][1]*uc[0][2])*uc[2][0]+
			( uc[0][2]*uc[1][0]-uc[0][0]*uc[1][2])*uc[2][1]+
			( uc[0][0]*uc[1][1]-uc[0][1]*uc[1][0])*uc[2][2]);
};


