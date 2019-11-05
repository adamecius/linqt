#!/usr/bin/python3

import time
import matplotlib.pyplot as plt
import numpy as np
np.show_config();
print("\n");
import scipy as sp
from numpy.linalg import eigvals,eigvalsh,multi_dot, norm
import copy 
import sys
from operators import bloch_op

from scipy.sparse import csr_matrix
def write_mycsr( outputname, Amat ):
    Amat = csr_matrix(Amat);

    with open(outputname , 'w') as f:
        f.write("%d %d\n" %  (max(Amat.shape), Amat.nnz) )
      
        for val in Amat.data:
            f.write("%1.10f %1.10f " % (np.real(val), np.imag(val)))
        f.write("\n");
    
        for row in Amat.indices:
            f.write("%d " % row)
        f.write("\n");

        for idx in Amat.indptr:
            f.write("%d " % idx)
        f.write("\n");


def get_user_parameters():
    numargs=len(sys.argv);
    if( numargs != 6 ):
        print( "Wannier2Sparse need to be called as: wannier2sparse.py LABEL D0 D1 D2 OP");
        print(  "\tLABEL" ,": string defining the filenames label.Ham and label.Geo");
        print(  "\tboth required to be present in the calling directory" );
        print(  "\t(D0,D1,D2): integers defining supercell dimension.");
        print( "\tIf omitted (1,1,1) is assumed.");
        print( "\tOP: string defining the output operator.");
        print( "\tOP only support the values Ham, Vx, Vy, Jx, Jy, Jxz, Jyz. If ommitted, Ham, is assumed.");
        sys.exit();

    label = sys.argv[1];
    scdim = ( int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]) )
    opname = sys.argv[5];
    hampath= label+".Ham";
    geopath= label+".Geo";
    print( "Wannier2Sparse will use ",hampath," and ",geopath)
    print( "for constructing a ",scdim[0],"x",scdim[1],"x",scdim[2], "unit cell.");
    print( "and returning the ",opname," operator");
    return label, scdim, opname;
    
def main():
    label, scdims, opname = get_user_parameters()

    #LOAD ALL DATA RELATED TO THE HAMILTONIAN OF THE UNIT CELL
    print( "Reading hamiltonian and geometrical information");
    Hk = bloch_op(label);
    #LOAD THE GEOMETRICAL INFORMATION FROM THE GEO FILE. 
    Hk.load_geofile(label);

    #Expand the hamiltonian into a supercell
    Hk.expand_to_supercell(scdims); 
    print(" Printing the ",opname, " operator in CSR format" );
    if( opname == "Ham" ):
        write_mycsr(label+".HAM.CSR", Hk.hamiltonian_operator() );
    elif( opname == "Vx" ):
        write_mycsr(label+".VX.CSR", Hk.current_operator("x") );
    elif( opname == "Jx" ):
        Jx = Hk.current_operator("x");
        Jx.data= 2*np.pi*Jx.data/Hk.volume/Hk.num_cells; 
        write_mycsr(label+".JX.CSR", Jx );    
    elif( opname == "Jy" ):
        Jy = Hk.current_operator("y");
        Jy.data= 2*np.pi*Jy.data/Hk.volume/Hk.num_cells; 
        write_mycsr(label+".JY.CSR", Jy );
    elif( opname == "JySz" ):
        JySz = Hk.spin_current_operator("y");
        JySz.data= 2*np.pi*JySz.data/Hk.volume/Hk.num_cells; 
        write_mycsr(label+".JYSZ.CSR", JySz );
    print(" FINISHED" );

    
if __name__ == "__main__":
    main()
    
    
