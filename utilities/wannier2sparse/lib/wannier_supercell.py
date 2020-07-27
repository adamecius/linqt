import wannier2sparse as w2s
import wannier_tools as wt
import numpy as np

class wannier_supercell(wt.wannier_system):

    def __init__(self, wan_syst, dims ):
        set_supercell_dims(dims);
 
    def __init__(self, label, dims ):
        super().__init__(label);
        self.set_supercell_dims(dims);

    def set_supercell_dims(self, dims ):
        self.dims = dims;
        n1,n2,n3 = dims;
        self.kpoints = ( np.mgrid[0:1:(n1+1)*1j, 0:1:(n2+1)*1j, 0:1:(n3+1)*1j].T)[:-1,:-1,:-1].reshape(n1*n2*n3,3);

    def save_op_asCSR(self, op_list=["HAM",] ):
        for op in op_list:
            print("Building the operator ", op)
            w2s.wannier2sparse(  self.label, self.dims, [op] );

