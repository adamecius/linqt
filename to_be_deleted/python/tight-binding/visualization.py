import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri #tricontourf
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi
import time


def check_ortogonality(eigvecs):
    U =(np.conj(eigvecs.T) ) .dot(eigvecs);
    U[np.abs(U) < 10*np.finfo(np.float64).eps] = 0;
    diff = np.linalg.norm(U - np.eye(len(U)))/len(U);
    if( diff > 10*np.finfo(np.float64).eps ):
        print( "The eigenvectors are not orthogonal. the difference is  ",diff,np.finfo(np.float64).eps)
        return False;
    return True;
    
def eigval_and_proj( Hk , Op, check_ortogonal = False ):
      
    #When no operator submited, return only the band structure
    if Op is None:
        return eigvalsh(Hk),np.ones( len(Hk) );

    #If not, compute the eigen vectors
    w, v = np.linalg.eigh( Hk );#The column w[:, i] is the normalized eigenvector of v[i] eigenvalue

    #if requested, check for orthogonality
    if( check_ortogonal ):
        check_ortogonality(v);

    w = np.diag( (np.conj(v.T) ) .dot(Hk.dot(v) ) ) ;
    p = np.diag( (np.conj(v.T) ) .dot(Op.dot(v) ) ) ;
    return np.real(w),np.real(p);

def compute_band_structure( Ham , band_paths, fermi_energy = 0.0, proj_op = None ):
    path_labels, path_points, paths = np.transpose(band_paths);
    num_paths = len(paths);

    kpoints   = list(); #List of kpoints used for the band structure calculation
    eigenvals = list(); #the eigenvalues of the bands
    eigenprojs= list(); #the eigenvalues of the bands
    label_index  = list(); #the last element of the path (used for assigning labels)

    #The paths are given in reciprocal lattice units a.b==2pi, but the k-hamiltonian is performed assuming
    # a.b==1. Therefore, when inserting a kp in the hamiltonian one should scale
    def scale(kp):
        return np.array(kp);
    
    #Compute the initial path point
    kp  = np.array(paths[0]);
    eigval, proj = eigval_and_proj( Ham(scale(kp)) , proj_op) ;
    eigenvals.append(eigval);
    eigenprojs.append(proj);
    kpoints.append(kp);
    label_index.append(0);

    #Compute all other path points
    kp_index = 0
    for p in range(1,num_paths):
        npoint= path_points[p];
        beg=paths[p-1];
        end=paths[p  ];
        for kp in np.linspace(beg,end,npoint)[1:]:
            eigval, proj = eigval_and_proj( Ham(scale(kp)) , proj_op ) ;
            eigenvals.append(eigval);
            eigenprojs.append(proj);
            kpoints.append( kp ); #notice the change of basis
            kp_index+=1;
        label_index.append( kp_index )

    bands = np.transpose(eigenvals);
    proj  = np.transpose(eigenprojs);
    
    return kpoints,bands,proj,label_index


def compute_kp_weigths_2D(syst, numkp , fermi_energy, broadening=0.1 ):
    assert len(numkp) == 2, "The kpoints cannot go in the third dimension";
    
    #Compute the eigenvalues in a kpoint grid
    energies= list();
    kpoints = list();     
    for kp in np.ndindex(numkp):
        shift = np.array(numkp)/2;
        kp = (np.array(kp)-shift)*2.0*np.pi/numkp 
        ham= syst.hamiltonian_submatrix(params=dict( k_x=kp[0],k_y=kp[1],k_z=0) )
        eigvals = eigvalsh(ham) ;
        energies.append(eigvals);
        kpoints.append(kp);
    #Use these energies, the Fermi energy, and the broadenings to compute 
    #the weigth of each eigenvalue.
    energies = np.array(energies);
    weigths = np.exp(-0.5*(energies - fermi_energy)**2/broadening**2 );  
    cap_below = 1E-7; #Below this value send to zero
    weigths[ weigths < cap_below] = 0;

    return kpoints,weigths;



def fermi_surface_2D(syst, numkp , fermi_energy, broadening=0.1, do_plot=True ):
    #define an array for the output
    fermi_surf_points = list();

    #Compute the weigth of each k point to the total energy at the fermi level
    kpoints, weigths = compute_kp_weigths_2D(syst, numkp , fermi_energy, broadening);

    #Average the weigths per band and join them with the k point.
    for i, kp in enumerate(kpoints):
        weigth = weigths[i];
        max_weight_pos = np.argmax(weigth);
        max_weigth=weigth[max_weight_pos];
        if ( max_weight_pos != 0 ):
            fermi_surf_points.append( (*(kp/2/np.pi), max_weigth) )

    #plot the final result
    if( len(fermi_surf_points) != 0 and do_plot ):
        print("Creating a plot for Fermi surface");
        X,Y,Z = np.transpose(fermi_surf_points)
        fig, ax, = plt.subplots()
        energy_contours = ax.tricontourf(X, Y, Z, np.arange(0, 1.0, .01),extend='both',antialiased=True, cmap = "Greys" )
        fig.colorbar(energy_contours, ax=ax)
        plt.show()    
    elif(len(fermi_surf_points) == 0):
        print( "System is in a gap, nothing to show")
    return fermi_surf_points;


def spin_texture_2D( syst, numkp , fermi_energy, broadening=0.1, do_plot=True):
    #define an array for the output
    spin_text_points = list();
    
    #IMPORTANT. THE HAMILTONIAN IS ASSUME TO BE in FORM ( ( Huu, Hud),( Hdu, Hdd) )
    #where u=up and d=down the spin direction;
    #Define spin operators
    def su_op( state, theta, phi ):
        dim = len(state)//2;
        new_state=np.zeros(2*dim, dtype=complex);
        new_state[range(dim)]        = cos(theta)*state[ range(dim) ]+ sin(theta)*exp(1j*phi)*state[ range(dim, 2*dim)];
        new_state[range(dim, 2*dim)] = sin(theta)*exp(-1j*phi)*state[ range(dim) ]- cos(theta)*state[ range(dim, 2*dim)];
        return new_state;
 
    def compute_spin_texture( state, spdir ):
        if( spdir =="x"):
            return np.real( np.vdot( state,su_op(state,np.pi/2,0) ) );
        if( spdir =="y"):
            return np.real( np.vdot( state,su_op(state,np.pi/2,np.pi/2) ) );
        if( spdir =="z"):
            return np.real( np.vdot( state,su_op(state,0,0) ) );
        return 0.0;
        
    #Compute the weigth of each k point to the total energy at the fermi level
    kpoints, weigths = compute_kp_weigths_2D(syst, numkp , fermi_energy, broadening);
    
    #Average the weigths per band and join them with the k point.
    for i, kp in enumerate(kpoints):
        max_weight_pos = np.argmax(weigths[i]);
        max_weigth=weigths[i][max_weight_pos];
        if ( max_weigth != 0 ):
            eigvs,eigfs = eigh(syst.hamiltonian_submatrix(params=dict( k_x=kp[0],k_y=kp[1], k_z=0)))
            state = eigfs[:,max_weight_pos];
            dirs = ("x","y","z");
            Sx = compute_spin_texture(state, "x");
            Sy = compute_spin_texture(state, "y");
            Sz = compute_spin_texture(state, "z");
            spin_text_points.append( ((*(kp/2/np.pi), Sx,Sy,Sz ) ) )
            
    if( len(spin_text_points) != 0 and do_plot ):
        print("Creating a plot for spin texture");
        X,Y,U,V,M = np.transpose(spin_text_points);
        fig, ax = plt.subplots()
        ax.set_title("Spin Textures")
        spin_texture =ax.quiver(X, Y, U, V, M, scale=10., cmap = "RdBu" )
        qk = ax.quiverkey(spin_texture, 0.60, 0.85, 1,"max", labelpos='E', coordinates='figure')
        fig.colorbar(spin_texture, ax=ax)
        plt.show()
    elif(len(spin_text_points) == 0):
        print( "System is in a gap, nothing to show")
    
    return spin_text_points;
