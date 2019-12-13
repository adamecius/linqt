#include<cassert>
#include<fstream>
#include<tuple>
#include<vector>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
using namespace std ;

int main( int argc, char* argv[]){    

    hopping_list hl_t;
    hopping_list::cellID_t cellID;
    hopping_list::edge_t  vertex_edge;
    hopping_list::value_t  hop;
    string tag ; 

    std::cout<<"Performing GRAPHENE HR TEST"<<std::endl;
    {
        const hopping_list::cellID_t SCDIM({1,1,1});
        const int unit_cell_basis_size = 2;

        //GRAPHENE TEST
        hopping_list hl_i = create_hopping_list(read_wannier_file("graphene_hr.dat"));
        hl_t = hopping_list();
        hl_t.SetWannierBasisSize(unit_cell_basis_size);
        cellID={ 0, 0, 0}; vertex_edge={0,1}; hop = hopping_list::value_t(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={0,1}; hop = hopping_list::value_t(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={0,1}; hop = hopping_list::value_t(1.0,0.12);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 0, 0}; vertex_edge={1,0}; hop = hopping_list::value_t(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={1,0}; hop = hopping_list::value_t(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={1,0}; hop = hopping_list::value_t(1.0,0.21);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={0,0}; hop = 1.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={0,0}; hop = 2.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={0,0}; hop = 3.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={0,0}; hop = hopping_list::value_t(4.123456,0.100000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 1, 0}; vertex_edge={0,0}; hop = hopping_list::value_t(5,-0.200000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1,-1, 0}; vertex_edge={0,0}; hop = hopping_list::value_t(1.0,1.0);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={1,1}; hop = 1.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={1,1}; hop = 2.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={1,1}; hop = 3.0;
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={1,1}; hop = hopping_list::value_t(4.123456,-0.100000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 1, 0}; vertex_edge={1,1}; hop = hopping_list::value_t(4.123456,-0.100000);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1,-1, 0}; vertex_edge={1,1}; hop = hopping_list::value_t(1.0,-1.0); 
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        assert( hl_t == hl_i );
    }
    std::cout<<"GRAPHENE HR TEST PASSED"<<std::endl<<std::endl;


    std::cout<<"Performing  WRAP IN SUPERCELL FOR UNIT_CELL TEST"<<std::endl;
    {
        hopping_list hl_i = create_hopping_list(read_wannier_file("graphene_hr.dat"));
        const hopping_list::cellID_t SCDIM({1,1,1});
        const int unit_cell_basis_size = 2;

        hl_t = hopping_list();
        hl_t.SetWannierBasisSize(unit_cell_basis_size);
        hl_t.cellSizes = SCDIM;

        cellID={0 , 0, 0}; vertex_edge={0,1};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={-1, 0, 0}; vertex_edge={0,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0,-1, 0}; vertex_edge={0,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={1,0};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 1, 0, 0}; vertex_edge={1,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0, 1, 0}; vertex_edge={1,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0,0,0}; vertex_edge={0,0};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={-1, 0, 0}; vertex_edge={0,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0,-1, 0}; vertex_edge={0,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 1, 0, 0}; vertex_edge={0,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0, 1, 0}; vertex_edge={0,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 1,-1, 0}; vertex_edge={0,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={-1, 1, 0}; vertex_edge={0,0};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={ 0, 0, 0}; vertex_edge={1,1};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={-1, 0, 0}; vertex_edge={1,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0,-1, 0}; vertex_edge={1,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 1, 0, 0}; vertex_edge={1,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0, 1, 0}; vertex_edge={1,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 1,-1, 0}; vertex_edge={1,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={-1, 1, 0}; vertex_edge={1,1};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        hl_i = wrap_in_supercell(SCDIM, hl_i );
        assert( hl_t == hl_i );
    }
    std::cout<<"WRAP IN SUPERCELL FOR UNIT_CELL TEST"<<std::endl<<std::endl;

    std::cout<<"PERFORMING WRAP IN SUPERCELL FINAL TEST"<<std::endl;
    ///WRAP IN SUPERCELL FINAL TEST: Insert the hopping elements of a one dimensional chain 
    /// and construct a supercell using wrap_in_supercell function.
    /// compare this with a supercell constructed inserting the hoppings directly
    /// 
    {
        const hopping_list::value_t hop_val( 1.0, 1.0 ); 
        const hopping_list::cellID_t SCDIM({2,2,1});
        const int unit_cell_basis_size = 2;
        hopping_list hl_uc;
        hl_uc.SetWannierBasisSize(unit_cell_basis_size);
        cellID={ 0, 0, 0}; vertex_edge={0,1}; hop = hop_val;
        hl_uc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={-1, 0, 0}; vertex_edge={0,1}; hop = hop_val;
        hl_uc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0,-1, 0}; vertex_edge={0,1}; hop = hop_val;
        hl_uc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 0, 0}; vertex_edge={1,0}; hop = std::conj(hop_val);
        hl_uc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 1, 0, 0}; vertex_edge={1,0}; hop = std::conj(hop_val);
        hl_uc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={ 0, 1, 0}; vertex_edge={1,0}; hop = std::conj(hop_val);
        hl_uc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        hl_uc = wrap_in_supercell(SCDIM, hl_uc );
       
        hopping_list hl_sc;
        hl_sc.SetWannierBasisSize( unit_cell_basis_size*SCDIM[0]*SCDIM[1]*SCDIM[2] );
        cellID={0,0,0}; vertex_edge={0,1}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={0,3}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={0,5}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,0}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,4}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={1,2}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={2,1}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={2,7}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={2,3}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={3,2}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={3,6}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={3,0}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={4,5}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={4,1}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={4,7}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={5,4}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={5,0}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={5,6}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={6,5}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={6,7}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={6,3}; hop = hop_val;
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={7,6}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={7,2}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        cellID={0,0,0}; vertex_edge={7,4}; hop = std::conj(hop_val);
        hl_sc.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );
        
        assert( hl_uc == hl_sc );
    }
    std::cout<<"WRAP IN SUPERCELL FINAL TEST PASSED"<<std::endl<<std::endl;
 

return 0;}
