#include<cassert>
#include<fstream>
#include<tuple>
#include<vector>
#include "hopping_list.hpp"
#include "wannier_parser.hpp"
#include "tbmodel.hpp"
using namespace std ;

int main( int argc, char* argv[]){    

    tbmodel model;

    model.readOrbitalPositions("spin_graphene.xyz");
    model.readUnitCell("spin_graphene.uc");
    model.readWannierModel("spin_graphene_hr.dat");    

    hopping_list JxSx= model.createHoppingSpinCurrents_list(0,'x');
    hopping_list JxSy= model.createHoppingSpinCurrents_list(0,'y');
    hopping_list JxSz= model.createHoppingSpinCurrents_list(0,'z');

    hopping_list Sx = model.createHoppingSpinDensity_list('x');
    hopping_list Sy = model.createHoppingSpinDensity_list('y');
    hopping_list Sz = model.createHoppingSpinDensity_list('z');

    hopping_list hl_i = model.hl;
    hopping_list hl_t;
    hopping_list::cellID_t cellID;
    hopping_list::edge_t  vertex_edge;
    hopping_list::value_t  hop;
    string tag ; 

    std::cout<<"Performing  WRAP IN SUPERCELL FOR UNIT_CELL TEST"<<std::endl;
    {
        const hopping_list::cellID_t SCDIM({1,1,1});
        const int unit_cell_basis_size = 4;

        hl_t = hopping_list();
        hl_t.SetWannierBasisSize(unit_cell_basis_size);
        hl_t.cellSizes = SCDIM;


        cellID={0 , 0, 0}; vertex_edge={0,2};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={2,0};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );


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

        cellID={0 , 0, 0}; vertex_edge={2,3};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={-1, 0, 0}; vertex_edge={2,3};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0,-1, 0}; vertex_edge={2,3};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={0 , 0, 0}; vertex_edge={3,2};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 1, 0, 0}; vertex_edge={3,2};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={ 0, 1, 0}; vertex_edge={3,2};
        hop += get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        cellID={0,0,0};
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );


        cellID={0,0,0}; vertex_edge={0,0};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={ 0, 0, 0}; vertex_edge={1,1};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={ 0, 0, 0}; vertex_edge={2,2};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );

        cellID={ 0, 0, 0}; vertex_edge={3,3};
        hop  = get<1>(hl_i.hoppings[get_tag(cellID,vertex_edge)]);
        hl_t.hoppings.insert( {get_tag(cellID,vertex_edge), hopping_list::hopping_t(cellID,hop,vertex_edge) } );


        hl_i = wrap_in_supercell(SCDIM, hl_i );
        assert( hl_t == hl_i );
    }
    std::cout<<"WRAP IN SUPERCELL FOR UNIT_CELL TEST"<<std::endl<<std::endl;

return 0;}
