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
        std::cout<<"Creating the Current matrix J"<<dir<<std::endl;
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
        
        hopping_list chl = this->hl ;
        for( auto& elem: chl.hoppings )
        {
            auto  tag  = get<0>(elem.second);
            if( tag != hopping_list::cellID_t({0,0,0} ) )//Send to zero  non diagonal elements
                get<1>(elem.second) *=  0.0 ;
        }

        return chl;
    };

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
        std::cout<<"Creating the spin Density matrix S"<<dir<<std::endl;
        auto id2spin = this->map_id2spin();
        auto dens = this->createHoppingDensity_list();
        for( auto& elem: dens.hoppings )
        {
            auto edge  = get<2>(elem.second);
            auto value =&get<1>(elem.second);
            auto s1=id2spin[edge[0]], s2= id2spin[edge[1]];
            if( s1!= 0 && s2!= 0 )
            switch(dir)
            {
                case 'x':
                if( s1 + s2 == 0  )
                    *value=1.0;
              	break;
            
                case 'y':
                if( s1 + s2 == 0  )
                    *value= hopping_list::value_t(0.0,s2);
        		break;

                case 'z':
                if( s1 == s2 ) 
                    *value= s2;
    		    break;

                default:
                    *value = 0 ;
            }
        }
        std::cout<<"Sucess."<<std::endl;
        return dens;
    }


    hopping_list createHoppingSpinCurrents_list(const int dir, const char sdir)
    {
        std::cout<<"Creating the spin Current matrix J"<<dir<<"S"<<sdir<<std::endl;
        auto id2spin = this->map_id2spin();
        auto curr = this->createHoppingCurrents_list(dir);
        for( auto& elem: curr.hoppings )
        {
            auto edge  = get<2>(elem.second);
            auto value =&get<1>(elem.second);
            auto s1=id2spin[edge[0]], s2= id2spin[edge[1]];
            if( s1!= 0 && s2!= 0 )
            switch(sdir)
            {
                case 'x':
                if( s1 + s2 == 0  )
                    *value *=1.0;
              	break;
            
                case 'y':
                if( s1 + s2 == 0  )
                    *value *= hopping_list::value_t(0.0,s2);
        		break;

                case 'z':
                if( s1 == s2 ) 
                    *value *= s2;
    		    break;

                default:
                    *value *= 0 ;
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
