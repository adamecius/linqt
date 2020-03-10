#include<limits>
#include<iomanip>
#include<cstdlib>
#include<cassert>
#include<fstream>
#include<tuple>
#include<vector>

#include "wannier_parser.hpp"
using namespace std ;

int main( int argc, char* argv[]){
    typedef tuple<string, array<double, 3> > xyz_elem;
    vector< xyz_elem > xyz_data;

    string  xyz_filename = "example.xyz";
    ofstream xyz_file(xyz_filename); 
    xyz_file.precision( std::numeric_limits<double>::digits10+2 );

    const int orbs_per_type = 10;
    const array<string,3> orb_type = {"A","B","C"};
    array<double,3> elem ;
    xyz_file<<orbs_per_type*orb_type.size()<<std::endl;

    for(auto type : orb_type)
    for(int i = 0 ; i < orbs_per_type; i++ )
    {
        elem = {10.*rand()/RAND_MAX,rand()/(double)RAND_MAX,(double)rand()};
        xyz_data.push_back(xyz_elem(type,elem ) );
        xyz_file<<type<<" "<<elem[0]<<" "<<elem[1]<<" "<<elem[2]<<std::endl;
    }xyz_file.close();

    assert( read_xyz_file(xyz_filename)  == xyz_data   );
    

return 0;}
