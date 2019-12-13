#include <deque>
#include <string>
#include <iostream>

#include "tbmodel.hpp"
#include "hopping_list.hpp"

using namespace std;

int main( int argc, char* argv[]){

deque< string > arguments(argv,argv+argc);
const string program_name = arguments[0]; arguments.pop_front();  

if( arguments.empty() )
{
    cerr<<"ERROR: The program: "<<program_name <<" should be called with arguments (LABEL, Dim0, Dim1, Dim2). "<<endl;
    return -1;
}

const string  label = arguments[0]; arguments.pop_front();  
cout<<"Using "<<label<<" as the system's identification label"<<endl
    <<"This label will be used to detect the label.xyz, label_hr.dat, and label.uc files"<<endl;

array<int, 3> cellDim={1,1,1};
for(int i = 0 ; i < 3 ; i++ )
{
    assert( !arguments.empty() );
    cellDim[i] = stoi(arguments[0]); arguments.pop_front();  
}

tbmodel model;

model.readOrbitalPositions(label+".xyz");
model.readUnitCell(label+".uc");
model.readWannierModel(label+"_hr.dat");    

std::cout<<"Creating the supercell ("<<cellDim[0]<<","<<cellDim[1]<<","<<cellDim[2]<<")"<<std::endl;
save_hopping_list_as_csr(label+".HAM.CSR"  , wrap_in_supercell(cellDim, model.hl) );
save_hopping_list_as_csr(label+".VX.CSR"   , wrap_in_supercell(cellDim, model.createHoppingCurrents_list(0)) );
save_hopping_list_as_csr(label+".VYSZ.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(1,'z')) );
save_hopping_list_as_csr(label+".SX.CSR"   , wrap_in_supercell(cellDim, model.createHoppingSpinDensity_list('x')) );
save_hopping_list_as_csr(label+".SY.CSR"   , wrap_in_supercell(cellDim, model.createHoppingSpinDensity_list('y')) );
save_hopping_list_as_csr(label+".SZ.CSR"   , wrap_in_supercell(cellDim, model.createHoppingSpinDensity_list('z')) ); 
std::cout<<"Supercells created successfully"<<std::endl;

cout<<"The programa finished"<<std::endl;
return 0;}
