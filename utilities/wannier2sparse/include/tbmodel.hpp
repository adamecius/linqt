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
#include<iostream>
#include<limits>
#include<algorithm>
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
        std::cout<<"Creating the Velocity matrix V"<<dir<<std::endl;
        assert( dir <3 && dir >=0 ); 
        assert(volume(lat_vecs) > 0 );
        assert_equal( (int)orbPos_list.size(), hl.WannierBasisSize());

        hopping_list chl = this->hl ;
        for( auto& elem: chl.hoppings )
        {
            auto  tag   = get<0>(elem.second);
            auto  edge  = get<2>(elem.second);

            auto orb_diff = get<1>(orbPos_list[edge[1]]) ;//save final orbital position
            
            for( int i=0; i < orb_diff.size(); i++)
                orb_diff[i] +=-get<1>(orbPos_list[edge[0]])[i];//substract initial position

            auto displ_vec = tag2cartesian(tag, lat_vecs); //get the displacement vector in cartesian
            for( int i=0; i < orb_diff.size(); i++)
                orb_diff[i] += displ_vec[i];    //add it to the orbital_difference

            //Change the hopping element accordingly
            get<1>(elem.second) *=  hopping_list::value_t( 0.0, -orb_diff[dir] );
        }
        return chl;
    };


    hopping_list createHoppingDensity_list()
    {
        assert(volume(lat_vecs) > 0 );
        assert(orbPos_list.size()==hl.WannierBasisSize());
        
        hopping_list chl ;
        chl.num_wann= hl.WannierBasisSize(); //number of wannier functions
             
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
					std::cout<<"("<<edge[0]<<","<<s1<<","<<o1<<") "<<"-->("<<edge[1]<<","<<s2<<","<<o2<<") ="<<*value<<std::endl;
					break;

					case 'y':
					*value = (s1 + s2 == 0 ? hopping_list::value_t(0.0,s2) : 0.0   );
					std::cout<<"("<<s1<<","<<o1<<") "<<"-->("<<s2<<","<<o2<<") ="<<*value<<std::endl;
					break;

					case 'z':
					*value = (s1 == s2 ? s2 : 0.0   );
					std::cout<<"("<<s1<<","<<o1<<") "<<"-->("<<s2<<","<<o2<<") ="<<*value<<std::endl;
					break;

					default:
					*value = 0 ;
				}
			}
			else 
				*value = 0 ;

        }

        std::cout<<"Sucess."<<std::endl;
        return dens;
    }


    hopping_list createHoppingSpinCurrents_list(const int dir, const char sdir)
    {
        std::cout<<"Creating the spin Velocity matrix V"<<dir<<"S"<<sdir<<std::endl;
        auto id2spin = this->map_id2spin();
        auto curr = this->createHoppingCurrents_list(dir);
        for( auto& elem: curr.hoppings )
        {
            auto edge  = get<2>(elem.second);
            auto value =&get<1>(elem.second);
            auto s1=id2spin[edge[0]], s2= id2spin[edge[1]];
            if( s1!= 0 && s2!= 0 ) //When no spinless
            {
				switch(sdir)
				{
					case 'x':
						std::cout<<"NOT IMPLEMENTED"<<std::endl;
					break;

					case 'y':
						std::cout<<"NOT IMPLEMENTED"<<std::endl;
					break;

					case 'z':
					*value *= (s1 == s2 ? s2 : 0.0   );
					break;

					default:
					*value *= 0 ;
				}
			}

        }
        std::cout<<"Sucess."<<std::endl;
        return curr;
    }


    int num_orbs;
    unitCell_t lat_vecs;
    vector< orbPos_t > orbPos_list;
    hopping_list hl;
};

#endif
