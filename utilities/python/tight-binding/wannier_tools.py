import numpy as np
import kwant
from kwant.builder import HoppingKind as hop

def wannier2kwant( ham_file, geo_file, fracional=True):
    data    = np.loadtxt( ham_file);     #Get the hamiltonian data
    orbs_pos= np.loadtxt(geo_file);      #Get the lattice constant and orbital position 
    lat_vec = np.transpose(orbs_pos[:3]);#Extract the lattice vector
    orbs_pos=  orbs_pos[3:];
    if( fracional ):
        orbs_pos = orbs_pos.dot(lat_vec);#Rewrite the lattice vectors in the reciprocal cell
    else:
        print("The orbitals are assumed to be in cartesian units");
    print(lat_vec)
    #Use wannier lattice as kwant lattice
    lat = kwant.lattice.general(lat_vec, orbs_pos) 
    orbs= lat.sublattices;
    #Define a periodic system:
    syst = kwant.Builder( kwant.TranslationalSymmetry(lat.vec((1,0,0)), lat.vec((0,1,0)), lat.vec((0,0,1))) )
    #Create a placeholder for all the sites in the system
    syst[(orb(0,0,0) for orb in orbs )] = 0.0;

    #Check if the matrix is onsite
    def is_onsite( data ):
        orb_i, orb_f, dx,dy,dz = np.array(data[:-2],dtype=int);
        return orb_i==orb_f and dx==dy==0;  
    #go through all elements and set onsites and hoppings when correspond
    for elem in data:
        #Include the appropiate onsite energies
        if( is_onsite(elem) ):
            orb_idx = np.array(elem[0],dtype=int);
            val     = np.array(elem[-2],dtype=float); 
            syst[orbs[orb_idx](0,0,0)] = val;
        #Include the appropiate hopping temrs
        else:
            oi,of,dx,dy,dz = np.array(elem[:5],dtype=int);            
            val   = np.array(elem[-2:],dtype=float); 
            syst[hop((dx, dy,dz), orbs[oi], orbs[of])]  = val[0] + 1j*val[1];

    return syst;
