import numpy as np
from scipy.sparse import coo_matrix


class bloch_op:

    def __init__(self, dim, rows, cols, vals, displ ):
        self.shape= np.array((dim,dim), dtype=int); #an operator is always a square matrix
        self.dim  = max(self.shape); #The dimension of the operator is simply its rank
        self.idxs = (rows,cols);#the idxs variables is what its use in the sparse matrices
        self.vals = np.array(vals, dtype=complex);#Convert to complex for efficiency
        self.displ= np.array(list(zip(displ[0],displ[1]))); #Convert to tuple for efficiency
        self.nnz  = len(rows); #the number of non zero elements

    def help(self ):
        print( "The block operator class posses the following important variables:") 
        print( "shape. Returns the shape of the operator") 
        print( "dim.  Returns the maximum dimension of the operator") 
        print( "idxs. Returs the initial and final orbital indexes") 
        print( "vals. Returs the values contained in the operator in coo format") 
        print( "displ. Returs the displacement vectors") 
        print( "lat_vecs. Returs the lattices vectors in absolute units") 
        print( "orbs_pos. Returs the position of the orbitals in absolute units") 
        print( "nnz. Returs the number of nonzero elements") 
        
    def define_unitcell(self, lat_vecs ):
        self.lat_vecs = np.array(lat_vecs) #storage internally the lattice vector
        self.rec_vecs = 2.0*np.pi*np.linalg.inv(lat_vecs) ; #use it to compute the reciprocal lattice vector
        self.displ = self.displ.dot( lat_vecs );
        self.volume= np.abs( np.linalg.det(lat_vecs) )

    def orbital_positions(self, orbs_pos):
        ox,oy = np.transpose(orbs_pos);
        self.orbs_pos = np.array(list(zip(ox,oy)));#Convert to tuple for efficiency
        #Include the orbital information intothe displacements
        opos = self.orbs_pos.dot(self.lat_vecs);
        self.displ+= [ opos[i_f] - opos[i_o] for (i_o,i_f) in np.transpose(np.array(self.idxs, dtype=int))]
        
    def eval_kop(self, kp ):
        phi_k= self.displ.dot(kp);
        return np.array( coo_matrix((self.vals*np.exp(1.0j*phi_k),self.idxs), shape=self.shape).todense() ); 

    def evalin_kgrid(self,scdim):
        scrv = self.rec_vecs/scdim; #reciprocal lattice vector of the supercell
        kps  = [ scrv.dot((i,j)) for i in range(scdim[0]) for j in range(scdim[1])]
        idx  = [ i*scdim[1] + j for i in range(scdim[0]) for j in range(scdim[1])]
        
        rows = np.empty(0)
        cols = np.empty(0)
        vals = np.empty(0)
        for i,kp in enumerate(kps):
            shift= idx[i]*self.dim;
            Ak = coo_matrix( self.eval_kop(kp) );
            rows = np.append(rows, Ak.row + shift);
            cols = np.append(cols, Ak.col + shift);
            vals = np.append(vals, Ak.data);

        new_shape = self.shape*scdim[0]*scdim[1];
        return coo_matrix((vals, (rows,cols)),shape=new_shape )

