import tbmodels
import numpy as np
from numpy import exp, cos, sin , pi
from numpy.linalg import eigh, eigvalsh, norm, inv

def tbgen2tbmodels( tbg_file, xyz_file, uc_file, hc_exist =True ):

    tbg_data= np.genfromtxt(tbg_file, dtype=None, encoding=None);     #Get the hamiltonian data
    xyz_data= np.genfromtxt(xyz_file, dtype=None, encoding=None);      #Get the lattice constant and orbital position 
    uc_data = np.loadtxt(uc_file);      #Get the lattice constant and orbital position 
   
    # Create a function that adds 100 to something
    on_site= np.array([rv for *R,oi,of,rv,im in tbg_data if oi==of and np.linalg.norm(R)==0 ]  ) 
    #Read the position file from the xyz file. It is assume that the file is sorted following the orbital indexes
    pos= np.array([pos for label, *pos in xyz_data ]  ) 

    #The orbital position should be pass in fractional units
    #but xyz file is in cartesian. Therefore, convert:
    pos = pos.dot(np.linalg.inv(uc_data) );
    
    
    #Construct the skeleton of the model from the previous data.
    model = tbmodels.Model(on_site=on_site, dim=3, occ=1, pos=pos)

    for *R,oi,of,rv,im in tbg_data:
        
        if( np.linalg.norm(R)==0  and oi==of ): #Onsite case already handled
            continue; 
        
        t =( rv + 1j*im );
        if( hc_exist ):	#If the hermitian conjugate exist, add a factor 1/2 because tbmodels include the hc automatically;
            t /= 2.0
        model.add_hop(t,oi-1,of-1, R);#the -1 is due to one-based -> zerp-based indexes

    return model;


def compute_band_structure( Ham , band_paths, fermi_energy = 0.0 ):
    path_labels, path_points, paths = np.transpose(band_paths);
    num_paths = len(paths);

    kpoints   = list(); #List of kpoints used for the band structure calculation
    eigenvals = list(); #the eigenvalues of the bands
    label_index  = list(); #the last element of the path (used for assigning labels)

    #The paths are given in reciprocal lattice units a.b==2pi, but the k-hamiltonian is performed assuming
    # a.b==1. Therefore, when inserting a kp in the hamiltonian one should scale
    def scale(kp):
        return np.array(kp);
    
    #Compute the initial path point
    kp  = np.array(paths[0]);
    eigenvals.append(eigvalsh(Ham(scale(kp))));
    kpoints.append(kp);
    label_index.append(0);
 
    #Compute all other path points
    kp_index = 0
    for p in range(1,num_paths):
        npoint= path_points[p];
        beg=paths[p-1];
        end=paths[p  ];
        for kp in np.linspace(beg,end,npoint)[1:]:
            eigenvals.append( eigvalsh( Ham(scale(kp)) ) );
            kpoints.append( kp ); #notice the change of basis
            kp_index+=1;
        label_index.append( kp_index )

    bands = np.transpose(eigenvals);
    
    return kpoints,bands,label_index
