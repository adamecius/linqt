#ifndef TBMODELS 
#define TBMODELS
#include <string>
#include <sstream>
#include <array>
#include <vector>
#include <tuple>
#include <complex>
#include <map>
#include <iostream>
#include <limits>
#include <cassert>
#include <functional>
#include "wannier_parser.hpp"
#include "hopping_list.hpp"

using namespace std;


template<typename T>
void assert_equal(const T x,const T y )
{
    if( x != y )
        std::cerr<<"ASSERT_EQUAL FAILED: "<<x<<" != "<<y<<std::endl;
    assert(x==y);
    return ;
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
}

class tbmodel
{
    public:
    typedef tuple<string, array<double, 3> > orbPos_t;
    typedef array < array<double,3> , 3 > unitCell_t;


    inline void readUnitCell(const string inputfile)
    {
        lat_vecs= read_unit_cell_file(inputfile);

    }

    inline void readOrbitalPositions(const string inputfile)
    {
        orbPos_list= read_xyz_file(inputfile);
    };

    inline void readWannierModel(const string inputfile)
    {
        hl = create_hopping_list(read_wannier_file(inputfile));
    };

    hopping_list createHoppingCurrents_list(const int dir)
    {
        assert( dir <3 && dir >=0 ); 
        assert(volume(lat_vecs) > 0 );
        assert_equal( (int)orbPos_list.size(), hl.WannierBasisSize());

        auto curr = this->hl ;
        for( auto& elem: curr.hoppings )
        {
            auto  cellID   = get<0>(elem.second);
            auto  edge  = get<2>(elem.second);

            auto orb_diff = get<1>(orbPos_list[edge[1]]) ;//save final orbital position
            
            for( int i=0; i < orb_diff.size(); i++)
                orb_diff[i] +=-get<1>(orbPos_list[edge[0]])[i];//substract initial position

            auto displ_vec = tag2cartesian(cellID, lat_vecs); //get the displacement vector in cartesian
            for( int i=0; i < orb_diff.size(); i++)
                orb_diff[i] += displ_vec[i];    //add it to the orbital_difference

            //Change the hopping element accordingly
            get<1>(elem.second) *=  hopping_list::value_t( 0.0, -orb_diff[dir] );           
        }
        return curr;
    };

    hopping_list createHoppingDensity_list()
    {
        assert(volume(lat_vecs) > 0 );
        assert(orbPos_list.size()==hl.WannierBasisSize());
        
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
				hopping_list::value_t hop_value(1.0,0.0);
				chl.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop_value,vertex_edge) } );
				j_orb++;
			}
			i_orb++;
		}

        return chl;
    };

    map<int,string>
    map_id2orbs()
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
    map_id2spin()
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


    hopping_list createHoppingSpinDensity_list(const char dir)
    {
		std::cout<<"Creating the  SpinDensity operator in "<<dir<<" direction"<<std::endl;
		auto id2spin = this->map_id2spin();
		auto id2orbs = this->map_id2orbs();
        
        auto dens = this->createHoppingDensity_list();

        for( auto& elem: dens.hoppings )
        {
            auto edge  = get<2>(elem.second);
            auto value =&get<1>(elem.second); 
            auto s1=id2spin[edge[0]], s2= id2spin[edge[1]];
            auto o1=id2orbs[edge[0]], o2= id2orbs[edge[1]];
			
            if( s1!= 0 && s2!= 0 && o1==o2 ) //When no spinless and diagonal in orbital index
            {
				switch(dir)
				{
					case 'x':
					*value = (s1 + s2 == 0 ? 1.0 : 0.0   );
					break;

					case 'y':
					*value = (s1 + s2 == 0 ? hopping_list::value_t(0.0,s2) : 0.0   );
					break;

					case 'z':
					*value = (s1 == s2 ? s2 : 0.0   );
					break;

					default:
					*value = 0 ;
				}
			}
			else 
				*value = 0 ;
        }
        return dens;
    }


    hopping_list createHoppingSpinCurrents_list(const int dir, const double theta, const double phi )
    {
		const std::complex<double> I(0,1);
		std::complex<double> sz,sval;
		
        hopping_list scurr ;
		scurr.SetWannierBasisSize(hl.WannierBasisSize()); //increase the supercell dimension;
		scurr.SetBounds(hl.Bounds() ); //increase the supercell dimension;

        for( const auto elem: this->createHoppingCurrents_list(dir).hoppings )
        {
            auto cellID= get<0>(elem.second);
            auto value = get<1>(elem.second);
            auto edge  = get<2>(elem.second);
	        int is0; getSpinIndex(edge[0],is0); 
            int is1; getSpinIndex(edge[1],is1); 
			auto sz  = 1.0-2.0*is0;

			// spin conserving hopping. There will always be an up and a down, So two elements are going to be defined
			if ( is0 == is1 ) 
			{
				//Create a spin conserving hopping. 
				sval= sz*cos(theta);
				if( std::norm(sval)>1e-5)
					scurr.AddHopping(cellID,sval*value,edge);	
			
				//Create a spin nonconserving hopping 
				sval= 0.5*sin(theta)*exp(-I*phi);
				if( std::norm(sval)>1e-5)
				{
					ChangeSpinIndex(edge[0], 0  );
					ChangeSpinIndex(edge[1], 1  ); //up->dn
					scurr.AddHopping(cellID,sval*value,edge); 
					//complex conjugate
					ChangeSpinIndex(edge[0], 1  );
					ChangeSpinIndex(edge[1], 0  ); //up->dn
					scurr.AddHopping(cellID,std::conj(sval)*value,edge);				
				}
			}
			else // spin nonconserving hopping
			{
				//Creates a spin conserving hopping. 
				sval=  0.5*sin(theta)*exp(sz*I*phi);
				if( std::norm(sval)>1e-5)
				for( int is =0 ; is<MAX_SPIN; is++)
				{
					ChangeSpinIndex(edge[0], is );
					ChangeSpinIndex(edge[1], is );//up->up
					scurr.AddHopping(cellID,sval*value,edge);				
				}
			}
		}

        return scurr;
    };

    hopping_list createHoppingSpinCurrents_list(const int dir, const char sdir)
    {
		switch(sdir)
		{
			case 'x':
				return createHoppingSpinCurrents_list(dir,0.5*M_PI,0.0*M_PI );
			break;
			
			case 'y':
				return createHoppingSpinCurrents_list(dir,0.5*M_PI,0.5*M_PI );
			break;
			
			case 'z':
				return createHoppingSpinCurrents_list(dir,0.0*M_PI,0.0*M_PI );
			break;

		}
    };


	inline 
	hopping_list add_onsite_disorder(const hopping_list  hl , double W=0.0 )
	{	
		auto id2spin = this->map_id2spin();
	//	for( auto& elem: dens.hoppings )
	//	{
	//		auto edge  = get<2>(elem.second);
	//		auto value =&get<1>(elem.second); 
	//		auto s1=id2spin[edge[0]], s2= id2spin[edge[1]];
	//		auto o1=id2orbs[edge[0]], o2= id2orbs[edge[1]];
			
	//		if( s1==s2 && edge[0]==edge[1] ) //When no spinless and diagonal in orbital index
		//		*value = W( rand()/RAND_MAX -0.5 );
		//}
		return hl;
	};

    int num_orbs;
    unitCell_t lat_vecs;
    vector< orbPos_t > orbPos_list;
    hopping_list hl;



	private :
	const int MAX_SPIN = 2; 
	inline //orbPos_list.size() gives you the number of positions. Including the repeated one due to MAX_SPIN
	const int NumberOfSites() const { return orbPos_list.size()/MAX_SPIN; } 
	inline 
	void getSpinIndex(const int n, int& s)const {  s = n/NumberOfSites() ; };
	inline
	void ChangeSpinIndex(int &n, const int s) const {  n = n%NumberOfSites() + s*NumberOfSites(); };


};

#endif
