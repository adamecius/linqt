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

model.readOrbitalPositions(label+".xyz"); std::cout<<" finished"<<std::endl;
model.readUnitCell(label+".uc"); std::cout<<" finished"<<std::endl;
model.readWannierModel(label+"_hr.dat"); 

std::cout<<"Creating the supercell ("<<cellDim[0]<<","<<cellDim[1]<<","<<cellDim[2]<<")"<<std::endl;
save_hopping_list_as_csr(label+".HAM.CSR"  ,model.add_onsite_disorder( wrap_in_supercell(cellDim, model.hl) ) );
std::cout<<"Creating VX"<<std::endl;
	save_hopping_list_as_csr(label+".VX.CSR"   , wrap_in_supercell(cellDim, model.createHoppingCurrents_list(0)) );
std::cout<<"Creating VY"<<std::endl;
	save_hopping_list_as_csr(label+".VY.CSR"   , wrap_in_supercell(cellDim, model.createHoppingCurrents_list(1)) );
std::cout<<"Creating VXSX"<<std::endl;
	save_hopping_list_as_csr(label+".VXSX.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(0,'x')) );
std::cout<<"Creating VYSX"<<std::endl;
	save_hopping_list_as_csr(label+".VYSX.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(1,'x')) );
std::cout<<"Creating VXSY"<<std::endl;
	save_hopping_list_as_csr(label+".VXSY.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(0,'y')) );
std::cout<<"Creating VYSY"<<std::endl;
	save_hopping_list_as_csr(label+".VYSY.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(1,'y')) );
std::cout<<"Creating VXSZ"<<std::endl;
	save_hopping_list_as_csr(label+".VXSZ.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(0,'z')) );
std::cout<<"Creating VYSZ"<<std::endl;
	save_hopping_list_as_csr(label+".VYSZ.CSR" , wrap_in_supercell(cellDim, model.createHoppingSpinCurrents_list(1,'z')) );
std::cout<<"Creating SX"<<std::endl;
	save_hopping_list_as_csr(label+".SX.CSR"   , wrap_in_supercell(cellDim, model.createHoppingSpinDensity_list('x')) );
std::cout<<"Creating SY"<<std::endl;
	save_hopping_list_as_csr(label+".SY.CSR"   , wrap_in_supercell(cellDim, model.createHoppingSpinDensity_list('y')) );
std::cout<<"Creating SZ"<<std::endl;
	save_hopping_list_as_csr(label+".SZ.CSR"   , wrap_in_supercell(cellDim, model.createHoppingSpinDensity_list('z')) ); 
std::cout<<"Supercells created successfully"<<std::endl;

cout<<"The programa finished"<<std::endl;
return 0;}
