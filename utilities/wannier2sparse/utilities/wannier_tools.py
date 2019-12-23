import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri #tricontourf
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi
from scipy.sparse import coo_matrix
import time

def load_xyz(filename):
    file = open(filename, "r");
    numOrbs = int(file.readline());

    xyz_data = list();
    for line in file:
        xyz_data.append([ elem for elem in line.split(' ') if elem != ''])
    xyz_data =np.array(xyz_data);
    OrbsID = xyz_data[:,0]
    xyz_coord = xyz_data[:,1:4].astype( float )

    return numOrbs,OrbsID, xyz_coord

def load_wannier_file(filename):

    file = open(filename, "r");
    date = file.readline()
    numOrbs = int(file.readline());
    numKPT  = int(file.readline());
    numKPT_lines = int(np.ceil(numKPT/15));
    file.close();

    #Read and format automatically the data
    wan_data= np.genfromtxt(filename, skip_header=numKPT_lines+3)   
    shift = (wan_data[:,0:3]).astype(int)
    rows  = wan_data[:,3:4].astype(int).flatten()-1
    cols  = wan_data[:,4:5].astype(int).flatten()-1
    values= wan_data[:,5:].astype(float)
    values= values[:,0]+1j*values[:,1]
 
    return shift,rows,cols,values


def _compute_band_structure( numOrbs,hvec,rows,cols,values , band_paths ):

    def Ham(k):
        data = values*np.exp(np.pi*2j*hvec.dot(k));
        return coo_matrix((data, (rows, cols)), shape=(numOrbs,numOrbs)).toarray();


    path_labels, path_points, paths = np.transpose(band_paths);
    num_paths = len(paths);

    kpoints   = list(); #List of kpoints used for the band structure calculation
    eigenvals = list(); #the eigenvalues of the bands
    label_index  = list(); #the last element of the path (used for assigning labels)

    #Compute the initial path point
    kp  = np.array(paths[0]);
    eigenval=  np.linalg.eigvalsh( Ham(kp) );
    eigenvals.append(eigenval);
    kpoints.append(kp);
    label_index.append(0);

    #Compute all other path points
    kp_index = 0
    for p in range(1,num_paths):
        npoint= path_points[p];
        beg=paths[p-1];
        end=paths[p  ];
        for kp in np.linspace(beg,end,npoint)[1:]:
            eigenval=  np.linalg.eigvalsh( Ham(kp) );
            eigenvals.append(eigenval);
            kpoints.append( kp ); #notice the change of basis
            kp_index+=1;
        label_index.append( kp_index )

    return np.array(kpoints),np.transpose(eigenvals),label_index
    

def compute_band_structure( label , band_paths ):

    wan_file =label+"_hr.dat"
    xyz_file =label+".xyz"
    uc_file  =label+".uc"
    lat_vec = np.loadtxt(uc_file);
    numOrbs,OrbsID, xyz_coord = load_xyz(xyz_file)
    hvec,rows,cols,values = load_wannier_file(wan_file)

    return compute_band_structure( numOrbs,hvec,rows,cols,values , band_paths );

def _density_of_states( numOrbs,hvec,rows,cols,values , kgrid, broadening=10 , num_eners=100):

    def Ham(k):
        data = values*np.exp(np.pi*2j*hvec.dot(k));
        return coo_matrix((data, (rows, cols)), shape=(numOrbs,numOrbs)).toarray();

    kxs,kys,kzs = [ np.linspace(0.0,1.0,ksize, endpoint=False) for ksize in kgrid ]; 
    
    eigvals = list();
    for kx in kxs:
        for ky in kys:
            for kz in kzs:
                kp=[kx,ky,kz];
                eigval =  np.linalg.eigvalsh( Ham(kp) );
                eigvals.append(eigval)
    
	
#	energies = np.linspace(min(eigvals), max(eigvals),num_eners);
    
    #broadening = broadening/1000;
   # dos = [(-1/np.pi/kgrid[0]/kgrid[1]/kgrid[2])*np.imag(np.sum( 1/(eigvals -(energy - 1j*broadening) ) ))  for energy in energies]
        
    return energies,dos

def density_of_states( label , band_paths ):

    wan_file =label+"_hr.dat"
    xyz_file =label+".xyz"
    uc_file  =label+".uc"
    lat_vec = np.loadtxt(uc_file);
    numOrbs,OrbsID, xyz_coord = load_xyz(xyz_file)
    hvec,rows,cols,values = load_wannier_file(wan_file)

    return _density_of_states( numOrbs,hvec,rows,cols,values , kgrid, broadening=100 , num_eners=100  );


