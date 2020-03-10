import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.tri as tri #tricontourf
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi, kron
from scipy.sparse import coo_matrix
import time

class wannier_system:
       
    def __init__(self, label ):
        wan_file =label+"_hr.dat"
        xyz_file =label+".xyz"
        uc_file  =label+".uc"
       
        self.lat_vec = np.loadtxt(uc_file);
        self.load_xyz(xyz_file)
        self.load_wannier_file(wan_file)

        #Convert the orbital positions in lattice vector basis
        scaled_coords = np.dot( self.xyz_coord,  np.linalg.inv(self.lat_vec) )
        coord2idx =  dict(zip(np.arange(self.numOrbs), scaled_coords)) 
        #Use the rows and columns to determine the position of the initial and 
        #final element in a hopping
        pos_i = np.array([ coord2idx[i] for i in self.rows]);
        pos_j = np.array([ coord2idx[j] for j in self.cols]);
        #Add both to the shift vectors to obtain a real position
        #in lattice vector units
        self.pos = self.shift+pos_j-pos_i;

        self.bandpath = np.array( ['G',1,[0,0,0]] );
        
    def load_xyz(self,filename):
        file = open(filename, "r");
        self.numOrbs = int(file.readline());

        xyz_data = list();
        for line in file:
            xyz_data.append([ elem for elem in line.split(' ') if elem != ''])
        xyz_data =np.array(xyz_data);
        self.OrbsID = xyz_data[:,0]
        self.xyz_coord = xyz_data[:,1:4].astype( float )


    def load_wannier_file(self,filename):

        file = open(filename, "r");
        date = file.readline()
        numOrbs = int(file.readline());
        numKPT  = int(file.readline());
        numKPT_lines = int(np.ceil(numKPT/15));
        file.close();

        #Read and format automatically the data
        wan_data= np.genfromtxt(filename, skip_header=numKPT_lines+3)   
        self.shift = (wan_data[:,0:3]).astype(int)
        self.rows  = wan_data[:,3:4].astype(int).flatten()-1
        self.cols  = wan_data[:,4:5].astype(int).flatten()-1
        values= wan_data[:,5:].astype(float)
        self.values= values[:,0]+1j*values[:,1]
       
        
    def spin_operator(self , comp=None ):
        s0 = [ [ 0 , 1  ] , [ 1 , 0 ] ];
        sx = [ [ 0 , 1  ] , [ 1 , 0 ] ];
        sy = [ [ 0 ,-1j ] , [ 1j, 0 ] ];
        sz = [ [ 1 , 0  ] , [ 0 ,-1 ] ];
        sop = {"0": s0 ,"x": sx, "y": sy, "z": sz};

        orb_dim = len(self.xyz_coord)//2;
        orb_ID = np.identity( orb_dim  ) ;
        
        if comp is None: #Return all components
            return [ np.kron(si,orb_ID) for si in [s0,sx,sy,sz] ] 
            
        if comp in sop:        
            SOP = np.kron(sop[comp],orb_ID);
            return SOP;
        else:
            print( "Nonexistent direction in spin_operator. Returning Identity" )
            return np.identity( len(self.xyz_coord) ) ;
        
    def ham_operator(self,k):
        data = self.values*np.exp(np.pi*2j*(self.pos).dot(k));
        return coo_matrix((data, (self.rows, self.cols)), shape=(self.numOrbs,self.numOrbs)).toarray();
        
    def projected_eigenvalues(self, k , proj_op  ):
        Hk = self.ham_operator(k);
        #When no operator submited, return only the band structure
        if proj_op is None:
            return eigvalsh( Hk ) ;

        #If not, compute the eigen vectors
        w, v = np.linalg.eigh( Hk );#The column w[:, i] is the normalized eigenvector of v[i] eigenvalue

        w = np.diag( (np.conj(v.T) ) .dot(Hk.dot(v) ) ) ;
        p = np.diag( (np.conj(v.T) ) .dot(proj_op.dot(v) ) ) ;
        return ( np.real(w),np.real(p) );

    def projections(self, k , proj_op  ):
        Hk = self.ham_operator(k);
        #When no operator submited, return only the band structure
        if proj_op is None:
            return eigvalsh( Hk ) ;

        #If not, compute the eigen vectors
        w, v = np.linalg.eigh( Hk );#The column w[:, i] is the normalized eigenvector of v[i] eigenvalue

        p = np.diag( (np.conj(v.T) ) .dot(proj_op.dot(v) ) ) ;
        return np.real(p);
    
    

    def Momentum_Rec2AbsMatrix(self ):
        return 2*np.pi* np.linalg.inv(self.lat_vec);
    
    
    def get_kpoints(self , absolute_coords = False): #Compute the kpoints used for the band structure calculation
        kpoints   = list(); 
        init_k = self.bandpath[0][2];
        for path_label, npoint, end_k in self.bandpath[1:]: #Not consider initial point anymore
                path_kpoints = np.linspace(init_k,end_k,npoint, endpoint=False ); 
                init_k  = end_k
                for kp in path_kpoints:
                    kpoints.append(kp)
        kpoints.append(init_k);

        #If required, rescale to absolute value
        if absolute_coords  is True :
             kpoints = np.dot( kpoints, np.transpose(self.Momentum_Rec2AbsMatrix() ) );
                
        return np.array(kpoints);

    def set_bandpath(self,bandpath):
        self.bandpath=np.array(bandpath);

    def get_XLabels(self ):
        Xaxis = self.bandsXaxis()
        xpos = (np.cumsum(self.bandpath[:,1],dtype=int)-1)
        return (Xaxis[xpos],self.bandpath[:,0])                
    
    
    def bandsXaxis(self ): #The kpoints are assume to be in absolute coordinates 
        kpoints=self.get_kpoints(absolute_coords = True);
        Xaxis=np.cumsum(np.linalg.norm(np.diff( kpoints,axis=0, prepend = 0),axis=1));
        Xaxis-=Xaxis[0];#Remove the initial value.
        return Xaxis;
    
    
    def compute_band_structure(self, fermi_energy = 0.0, proj_op = None, ax=None, plot_proj=False, proj_range=None ):
     
        #Compute the kpoints based on the path
        kpoints= self.get_kpoints();

        #Compute the eigenvalues and the projected values 
        peigenvals = np.array( [self.projected_eigenvalues( kp , proj_op) for kp in kpoints ] );

        if proj_op is None:
            return self.plot_band_structure( bands = np.transpose(peigenvals-fermi_energy), ax=ax );

        bands = np.transpose(peigenvals[:,0]-fermi_energy );
        projs = np.transpose(peigenvals[:,1] );
        
        return self.plot_band_structure( bands = bands , projs = projs , ax=ax, plot_proj = plot_proj, proj_range=proj_range);

    
    def plot_band_structure(self, bands, projs = None, ax=None, plot_proj = False, proj_range=None):
        
        xaxis  = self.bandsXaxis();
            
        #PLOTING
        #plot options
        if  ax is None:
            fig = plt.figure();
            _ax = fig.add_subplot();
        else:
            _ax  = ax;

        _ax.tick_params(axis='both', which='major', labelsize=16);
        
        if plot_proj is True:
            for i,proj in enumerate(projs):
                _ax.plot(xaxis,proj );
        else:
            _ax.set_ylabel("Energy (eV)", fontsize=16);
            cmap_name="seismic";
            vmin,vmax = 0,1;
            if projs is not None:
                vmin = np.min(projs);
                vmax = np.max(projs);
                if proj_range is not None:
                    vmin,vmax = list(proj_range);
            for i,band in enumerate(bands):
                _ax.plot(xaxis,band , color = 'k');
                c = 'k';
                if projs is not None:
                    c = projs[i];
                #Add points to the plot
                im = _ax.scatter(xaxis, band,s=50, c=c,cmap=cmap_name,vmin=vmin, vmax=vmax);
                im.set_facecolor("none");        
            #Create falso plot for color bar
            if projs is not None and ax is None: #Add the color bar whenever you have the a projectedprlt
                im = _ax.scatter(xaxis, bands[0],s=0, c=np.linspace( vmin,vmax,len(bands[0]) ),cmap=cmap_name);
                fig.colorbar(im, ax=_ax); 

        klabels, path_labels = self.get_XLabels()
        _ax.set_xticks(klabels);
        _ax.set_xticklabels(path_labels);
                
        if ax is None:

            return fig,_ax;

        ax = _ax;
        return ax;


    def compute_projections(self,kpoints, proj_op ):
        return  np.array( [self.projections( kp, proj_op=proj_op ) for kp in kpoints ] );
    
    def compute_dispersion(self,kpoints):
        return np.array( [self.projected_eigenvalues( kp, proj_op=None ) for kp in kpoints ] );        
      
    def compute_fermi_surface(self,kpoints, fermi_energy = 0.0, tol = None ,  proj_op = None ):

        #Compute the eigenvalues and the projected values 
        eigvals = np.array( [self.projected_eigenvalues( kp, proj_op=None ) for kp in kpoints ] );
    
        #Determine the important kpoints
        relevant_kpoints = np.any(np.abs(eigvals - fermi_energy) < tol, axis=1);
        #Select the kpoints
        kpoints = kpoints[relevant_kpoints];

        if proj_op is None:
            return kpoints,eigvals[relevant_kpoints];

        #If one requires a projection, select the kpoints and then use them to compute the projections
        
        peigvals = np.array( [self.projected_eigenvalues( kp, proj_op=proj_op ) for kp in kpoints ] );
        return kpoints,peigvals[relevant_kpoints];

    
    def refine_kpoints(self, kpoints):

        nkp= 2*len(kpoints);
        kmin = np.array([ np.min(kpoints[:,0]),np.min(kpoints[:,1]),0 ] );
        kmax = np.array([ np.max(kpoints[:,0]),np.max(kpoints[:,1]),0 ] );
        kpoints= np.array(list(np.ndindex((nkp,nkp,1))))/(nkp-1)*(kmax-kmin) + kmin ;#kpoints in recpricola lattice vectors

        return kpoints;
    
    
    def compute_DOS(self, energies, fermi_energy = 0.0, kpoints =None  ):
        b1,b2,b3 = self.Momentum_Rec2AbsMatrix().T
        nkp= 2*len(kpoints);
        kmin = np.array([ np.min(kpoints[:,0]),np.min(kpoints[:,1]),0 ] );
        kmax = np.array([ np.max(kpoints[:,0]),np.max(kpoints[:,1]),0 ] );
        kpoints= np.array(list(np.ndindex((nkp,nkp,1))))/(nkp-1)*(kmax-kmin) + kmin ;#kpoints in recpricola lattice vectors

        return kpoints;
        

"""
        #Convert kpoints from reciprocal to real coordinates
        kpoints = np.dot( kpoints, np.transpose(self.lat_vec/(2*np.pi))  );
        kaxis   = np.cumsum( np.insert( np.linalg.norm(kpoints[1:]-kpoints[:-1],axis=1),0,0) );
        klabels = [ kaxis[label] for label in label_index]

        #Transposition of bands is different depending 
        #on whereas we computed the projected band structure or not
        if proj_op is None:
            bands = np.transpose(peigenvals);
        else:
            bands = np.transpose(peigenvals[:,0] );
            projs = np.transpose(peigenvals[:,1] );

            
        #PLOTING
        #plot options
        fig = plt.figure();
        ax  = fig.add_subplot();
     
        ax.set_xticks(klabels);
        ax.set_xticklabels(path_labels);
        ax.tick_params(axis='both', which='major', labelsize=18);
        ax.set_ylabel("Energy (eV)", fontsize=18);
        cmap_name="seismic";
        vmin = np.min(projs);
        vmax = np.max(projs);
        
        for i,band in enumerate(bands):
            ax.plot(kaxis,band , color = 'k');
            c = 'k';
            if proj_op is not None:
                c = projs[i];
            #Add points to the plot
            im = ax.scatter(kaxis, band,s=50, c=c,cmap=cmap_name,vmin=-1, vmax=1);
            im.set_facecolor("none");
        
        #Create falso plot for color bar
        im = ax.scatter(kaxis, bands[0],s=0, c=np.linspace( vmin,vmax,len(bands[0]) ),cmap=cmap_name);
        fig.colorbar(im, ax=ax);
                                                                   
        return fig,ax;


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
        

def _density_of_states( numOrbs,
,rows,cols,values , kgrid, broadening=10 , num_eners=100):

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
"""    