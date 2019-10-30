import numpy as np
from scipy.sparse import coo_matrix
from numpy import exp, pi, ndindex

#TO DO:There are functions in this class which assume the operator is 2Dm search with 2DFunc
class bloch_op:

    def __init__(self, system_label, use_spin = True ):
        spatial_dim = 3;
        self.scdims = (1,1,1); 
        self.num_cells = np.prod( self.scdims );
        self.use_spin = use_spin;

        rows, cols,DX, DY,DZ, ret,imt = np.loadtxt(system_label+".Ham", unpack=True); 
        vals = ret + 1j*imt
        self.from_orbs = np.array(rows, dtype= int);#the idxs variables is what its use in the sparse matrices
        self.to_orbs   = np.array(cols, dtype= int);#the idxs variables is what its use in the sparse matrices
        self.hoppings  = np.array(vals, dtype=complex);#Convert to complex for efficiency
        self.shifts    = np.array([DX,DY,DZ],dtype = int); #Convert to tuple for efficiency   
        self.displs    = np.transpose(self.shifts); #Convert to tuple for efficiency   

        self.num_orbs = int(np.max([rows,cols])) +1 ;
        self.orbs_pos = np.zeros( (self.num_orbs,spatial_dim),dtype=float );            

        #Save the spin transitions
        spin_idx_tran = self.num_orbs/2;
        spin_beg = np.floor( self.from_orbs/spin_idx_tran );
        spin_end = np.floor( self.to_orbs/spin_idx_tran );
        spin_trans = np.array( np.transpose([spin_beg , spin_end]) ,dtype=int );
        self.sz_vals = [  (1 - 2*si*so)*(si==so)  for (si,so) in spin_trans  ]

    def define_unitcell(self, lat_vecs ):
        self.lat_vecs = np.array(lat_vecs) #storage internally the lattice vector
        self.rec_vecs = 2.0*pi*np.linalg.inv(lat_vecs) ; #use it to compute the reciprocal lattice vector
        self.volume   = np.abs( np.linalg.det(lat_vecs) )
        self.displs    = self.displs.dot( self.lat_vecs);

    def load_geofile(self, system_label ):
        Ro = np.loadtxt(system_label+".Geo", dtype=float);
        lat_vecs = Ro[:3]; #extract the lattice vectors;
        Ro = Ro[3:] #keep only the orbitals
        self.num_orbs = len(Ro); 
        self.define_unitcell(lat_vecs );
        self.orbs_pos = Ro.dot(self.lat_vecs);
        self.displs+=self.orbs_pos[self.to_orbs] - self.orbs_pos[self.from_orbs] ;
 
    def expand_to_supercell(self, scdims ):
        self.scdims = scdims; 
        self.num_cells = np.prod( self.scdims );                
        hopps= np.zeros( (self.num_cells, len(self.hoppings)  ), dtype=complex);
        forbs= np.zeros( (self.num_cells, len(self.from_orbs)), dtype=int);
        torbs= np.zeros( (self.num_cells, len(self.to_orbs)  ), dtype=int);
        for idx,kp in enumerate( ndindex(tuple(scdims)) ): #go throw all the sc indexes
            kp  = self.rec_vecs.dot( kp )/scdims;
            hopps[idx]= self.hoppings*np.exp( 1j*self.displs.dot(kp) ) ;
            forbs[idx]= self.from_orbs+ idx*self.num_orbs;
            torbs[idx]= self.to_orbs  + idx*self.num_orbs;
        self.hoppings = np.array(hopps).flatten();
        self.from_orbs= np.array(forbs).flatten();
        self.to_orbs  = np.array(torbs).flatten();

    def hamiltonian_operator(self):
        dim = self.num_cells*self.num_orbs;
        indexes = (self.from_orbs,self.to_orbs);
        vals = self.hoppings;
        Amat = coo_matrix((vals,indexes), shape=(dim,dim) )  ;
        Amat.sum_duplicates();
        Amat.eliminate_zeros();
        if( (Amat != Amat.getH()).nnz==0 ):
            print( "Non hermitian hamiltonian operator")
        return Amat;

    def current_operator(self, cur_dir="x"):
        dir2int = { "x":0 , "y":1, "z":2 };
        i = dir2int[cur_dir];
        dim = self.num_cells*self.num_orbs;
        indexes = (self.from_orbs,self.to_orbs);
        vals= 1j*self.hoppings*np.tile(np.transpose(self.displs)[i], self.num_cells );
        Amat= coo_matrix((vals,indexes ), shape=(dim,dim) )  ;
        Amat.sum_duplicates();
        Amat.eliminate_zeros();
        if( (Amat != Amat.getH()).nnz==0 ):
            print( "Non hermitian hamiltonian operator")
        return Amat;

    def spin_current_operator(self, cur_dir="x"):
        dir2int = { "x":0 , "y":1, "z":2 }; 
        i = dir2int[cur_dir];
        dim = self.num_cells*self.num_orbs;
        indexes = (self.from_orbs,self.to_orbs);
        vals =1j*self.hoppings*np.tile(np.transpose(self.displs)[i], self.num_cells );
        vals*=np.tile(self.sz_vals, self.num_cells);
        Amat= coo_matrix((vals,indexes ), shape=(dim,dim) )  ;
        Amat.sum_duplicates();
        Amat.eliminate_zeros();
        if( (Amat != Amat.getH()).nnz==0 ):
            print( "Non hermitian hamiltonian operator")
        return Amat;
