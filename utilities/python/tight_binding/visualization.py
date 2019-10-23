import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri #tricontourf
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi
import time


def compute_band_structure(syst, band_paths, do_plot =True ):
    path_labels, path_points, paths = np.transpose(band_paths);
    num_paths = len(paths);

    kvec  = list();
    eigenvals = list();
    path_end = list()
    
    #Compute the initial path point
    seconds = time.time()
    k0,k1,k2  = np.array(paths[0])*2*np.pi;
    ham = syst.hamiltonian_submatrix(params=dict( k_x=k0,k_y=k1, k_z=k2) );
    eigenvals.append( eigvalsh(ham) );
    kvec.append(0);
    path_end.append(0);
    seconds -= time.time()
    print("It took ",-seconds," to compute the initial point. You do the math")

    #Beforme moving, it is convenient to compute the reciprocal lattice vectors:
    #which can be done in a single instruction
    rec =np.linalg.pinv( np.array(syst._wrapped_symmetry.periods).T ).T 
    
    #Compute all other path points
    for p in range(1,num_paths):
        npoint= path_points[p];
        beg=np.array(paths[p-1])*2*np.pi;
        end=np.array(paths[p  ])*2*np.pi;
        kp0 = beg;
        for kp in np.linspace(beg,end,npoint)[1:]:
            ham = syst.hamiltonian_submatrix(params=dict( k_x=kp[0],k_y=kp[1], k_z=kp[2]) )
            eigenvals.append( eigvalsh(ham) );
            kvec.append( kvec[-1]+norm( (kp-kp0).dot(rec) )); #notice the change of basis
            kp0 = kp;
        path_end.append( kvec[-1] )

    bands=  np.transpose(eigenvals);    

    for i,band in enumerate(bands):
        plt.plot( kvec,band);   
    labels = [ band_path[0]  for band_path in band_paths]
    plt.xticks(path_end,labels, color='k', size=20)
    plt.show();
    
    
    return kvec,bands,path_end


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
