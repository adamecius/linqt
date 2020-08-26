import numpy as np
from numpy.linalg import eigh, eigvalsh, norm
from numpy import exp, cos, sin , pi, kron

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.tri as mtri #tricontourf
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)

from scipy.sparse import coo_matrix
import time


proj_s = 50 ;
cmap_name="seismic";
band_lw = 1;
class wannier_system:
       
    def __init__(self, label ):

        self.label =label;

        wan_file =self.label+"_hr.dat"
        xyz_file =self.label+".xyz"
        uc_file  =self.label+".uc"
       
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
        self.numOrbs = int(file.readline());
        self.numKPT  = int(file.readline());
        self.numKPT_lines = int(np.ceil(self.numKPT/15));
        self.KPT =[];
        for n in range(self.numKPT_lines):
            self.KPT.append(file.readline());
        file.close();

        #Read and format automatically the data
        wan_data= np.genfromtxt(filename, skip_header=self.numKPT_lines+3)   
        self.shift = (wan_data[:,0:3]).astype(int)
        self.rows  = wan_data[:,3:4].astype(int).flatten()-1
        self.cols  = wan_data[:,4:5].astype(int).flatten()-1
        values= wan_data[:,5:].astype(float)
        self.values= values[:,0]+1j*values[:,1]


    def save_wannier_file(self,filename):
        from datetime import date
        
        file = open(filename, "w");
        file.write(" written on  %s\n"%date.today())
        file.write("          %d\n"%self.numOrbs)
        file.write("        %d\n"%self.numKPT)
        for line in self.KPT:
            file.write("%s"%line)

        for n,s in enumerate(self.shift):
            v  = self.values[n]
            r,c=self.rows[n]+1,self.cols[n]+1
            file.write("    %d    %d    %d    %d    %d    %f    %f\n"%(s[0],s[1],s[2],r,c,np.real(v),np.imag(v) ) )
        file.close();

    def save_xyz(self,filename):
        file = open(filename, "w");
        file.write("          %d\n"%self.numOrbs)
        for n,xyz in enumerate(self.xyz_coord):
            s=self.OrbsID[n]
            x,y,z = xyz
            file.write("%s %f %f %f\n"%(s,x,y,z) )
        file.close();
        
    def save(self,label):
        self.save_wannier_file(label+"_hr.dat");
        self.save_xyz(label+".xyz")
        np.savetxt(label+".uc",self.lat_vec);

       
    def num_spins(self):
        return 2;

    def orbs_per_spins(self):
        return self.numOrbs//self.num_spins();
        
    #This function returns the indexes of the spin-dependent onsite terms
    def get_onsites_indexes(self):
        num_spins = self.num_spins();
        orbsPerSpin = self.orbs_per_spins();

        spinless_rows = self.rows%orbsPerSpin; #The orbital index without the spins for the rows
        spinless_cols = self.cols%orbsPerSpin; #The orbital index without the spins for the cols
        spinless_onsites = (spinless_rows-spinless_cols)==0;
        intracell_hop= ( np.sum( np.abs(self.shift),axis=1)==0 ); # Indexes of the intracell hoppings

        return np.arange(len(self.rows))[intracell_hop*spinless_onsites];#The position of these indexes in the wannier system

    #This function returns the orbital index, spin index and value of the the spin-dependent onsite terms
    def get_onsites(self):
        num_spins = self.num_spins();
        orbsPerSpin = self.orbs_per_spins();

        onsiteIdxs= self.get_onsites_indexes();
        spinIdx_rows = self.rows[onsiteIdxs]//orbsPerSpin; #The spin index for a particular orbital in the rows
        spinIdx_cols = self.cols[onsiteIdxs]//orbsPerSpin; #The spin index for a particular orbital in the cols
        onsite_values= self.values[onsiteIdxs]; #The values of the zeeman field
 
        spinless_rows = self.rows%orbsPerSpin; #The orbital index without the spins for the rows
        spinless_diag= spinless_rows[onsiteIdxs]; #This is the diagonal index in the spinless orbital space
 
        #Get the indexes associated with the zeeman index
        onsite_indexes = np.transpose([spinless_diag,spinIdx_rows, spinIdx_cols])

        #Sort both the indexes and the values in terms of the orbital spinless indexes
        #sorted_idx = np.argsort(spinless_diag);
        #onsite_indexes= onsite_indexes[sorted_idx];
        #onsite_values = onsite_values[sorted_idx];
        return list(zip(map(tuple,onsite_indexes),onsite_values))
    
    #This function set the orbital index, spin index and value of the the spin-dependent onsite terms
    def set_onsites(self, spinor_diag):
        onsite_indexes = self.get_onsites_indexes();
        to_idx = { x:onsite_indexes[i] for i,(x,v) in enumerate(self.get_onsites()) };
        for n,En in enumerate(spinor_diag):
            for s1,s2 in ((0,0),(0,1),(1,0),(1,1)):
                self.values[to_idx[(n,s1,s2)]] = En[s1,s2] ;
        return 0;

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
            return ( eigvalsh( Hk ) ) ;

        #If not, compute the eigen vectors
        w, v = np.linalg.eigh( Hk );#The column w[:, i] is the normalized eigenvector of v[i] eigenvalue

        w = np.diag( (np.conj(v.T) ) .dot(Hk.dot(v) ) ) ;
        p = np.diag( (np.conj(v.T) ) .dot(proj_op.dot(v) ) ) ;
        return ( np.real(w),np.real(p) );


    def Momentum_Rec2AbsMatrix(self ):
        return 2*np.pi* np.linalg.inv(self.lat_vec);
    
    #If required, rescale to absolute value
    def toAbsoluteCoords(self,x):
             return np.dot( x, np.transpose(self.Momentum_Rec2AbsMatrix() ) );
                
    def band_kpoints(self , absolute_coords = False): #Compute the kpoints used for the band structure calculation
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

    def set_bandpath(self,bandpath, absolute_coords = False):

        if absolute_coords is True: #Convert to reciprocal
            Abs2Rec = np.linalg.inv(np.transpose(self.Momentum_Rec2AbsMatrix() ) );
            bandpath = [ (x[0],x[1],np.dot(x[2],Abs2Rec)) for x in bandpath]

        self.bandpath=np.array(bandpath);

    def get_XLabels(self ):
        Xaxis = self.bandsXaxis()
        xpos = (np.cumsum(self.bandpath[:,1],dtype=int)-1)
        return (Xaxis[xpos],self.bandpath[:,0])                
    
    
    def bandsXaxis(self ): #The kpoints are assume to be in absolute coordinates 
        kpoints=self.band_kpoints(absolute_coords = True);
        Xaxis=np.cumsum(np.linalg.norm(np.diff( kpoints,axis=0, prepend = 0),axis=1));
        Xaxis-=Xaxis[0];#Remove the initial value.
        return Xaxis;
    
    
    def compute_band_structure(self, fermi_energy = 0.0, proj_op = None, ax=None, plot_proj=False, proj_range=None ):
     
        #Compute the kpoints based on the path
        kpoints= self.band_kpoints();

        #Compute the eigenvalues and the projected values 
        peigenvals  = self.compute_dispersion( kpoints , proj_op)  ;

        if proj_op is None:
            bands = peigenvals.T- fermi_energy;
            return self.plot_band_structure( bands = bands, ax=ax );

        bands = peigenvals[:,0].T - fermi_energy;
        projs = peigenvals[:,1].T;
        
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
            vmin,vmax = 0,1;
            if projs is not None:
                vmin = np.min(projs);
                vmax = np.max(projs);
                if proj_range is not None:
                    vmin,vmax = list(proj_range);
            for i,band in enumerate(bands):
                _ax.plot(xaxis,band , color = 'k', lw=band_lw);
                c = 'k';
                if projs is not None:
                    c = projs[i];
                #Add points to the plot
                im = _ax.scatter(xaxis, band,s=proj_s, c=c,cmap=cmap_name,vmin=vmin, vmax=vmax);
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

    def compute_dispersion(self,kpoints, proj_op=None, absolute_coords = False):
        if ( absolute_coords is True ):
            Abs2Rec = np.linalg.inv(np.transpose(self.Momentum_Rec2AbsMatrix() ) );
            kpoints = np.dot( kpoints , Abs2Rec);
        
        return np.array( [self.projected_eigenvalues( kp, proj_op ) for kp in kpoints ] );        
      
    def compute_fermi_surface(self,kpoints, fermi_energy = 0.0, tol = None ,  proj_op = None ):

        #Compute the eigenvalues and the projected values 
        peigenvals  = self.compute_dispersion(kpoints, proj_op=proj_op)
        eigenvalues = np.array(list(map(list,peigenvals[:,0])))
        eigenvalues-= fermi_energy;#Compute the eigenvalues

        #Determine the important kpoints
        relevant_kpoints = np.any(np.abs(eigenvalues - fermi_energy) < tol, axis=1);
        #Select the kpoints
        kpoints = kpoints[relevant_kpoints];

        if proj_op is None:
            return kpoints,eigvals[relevant_kpoints];

        #If one requires a projection, select the kpoints and then use them to compute the projections
        
        peigvals = np.array(list(map(list,peigenvals[:,1])))
        return kpoints,peigvals[relevant_kpoints];

   
    def compute_DOS(self, energies,gridpoints, fermi_energy = 0.0, broadening = 0.1,  kpoints = None ,proj_op = None ):
        
        if kpoints is None:
            n1,n2,n3 = gridpoints
            kpoints = ( np.mgrid[0:1:(n1+1)*1j, 0:1:(n2+1)*1j, 0:1:(n3+1)*1j].T)[:-1,:-1,:-1].reshape(n1*n2*n3,3);

        eta=broadening;
        def gaussian(x):
            return np.exp( -(x/eta)**2/2 )/eta/np.sqrt(2*np.pi)

        dim = len(kpoints);
        peigenvals  = self.compute_dispersion(kpoints, proj_op=proj_op)
        eigenvalues = np.array(list(map(list,peigenvals[:,0]))).flatten()
        eigenvalues -= fermi_energy;#Compute the eigenvalues
        
        
        if proj_op is None:
            return np.array( [np.sum(gaussian(eigenvalues-EF)) for EF in energies] )/dim;

        P = np.array(list(map(list,peigenvals[:,1]))).flatten()
        return np.array( [np.sum(P*gaussian(eigenvalues-EF)) for EF in energies] )/dim;


    def band_spin_texture(self, kpoints, bidx, uop = None, vop = None, zop=None, zlims = None, ax=None ):
        #Compute first the mean values using the model

        Ops= [uop, vop, zop ] ;

        def get_band_data( x, bidx ):
            y = None;
            if x.shape[0] == kpoints.shape[0]:
                y = x;
            else:
                y = self.compute_dispersion(kpoints, proj_op = x )[:,1];
            y = np.array(list(map(list,y))).T;
            return y[bidx];

        #If is a compatible array of data, add it directly to Ops
        Ops= [ get_band_data(op,bidx) for op in Ops ];

        uop, vop, zop = Ops;
        KX,KY,KZ = np.transpose(kpoints);

        # Create triangulation.
        triang = mtri.Triangulation(KX, KY)

        def get_mask( x, y, z ):
            return mtri.CubicTriInterpolator(triang, z, kind='geom')(x, y) > 0.;

        def xypoints( npts ):
            return np.meshgrid(np.linspace(KX.min(), KX.max(), npts), np.linspace(KY.min(), KY.max(), npts));

        # Set up the figure
        if( ax is None ):
            fig, ax = plt.subplots(dpi=400);

        # Interpolate to regularly-spaced quad grid.
        xi, yi = xypoints(npts=100)
        fzi    = mtri.CubicTriInterpolator(triang, zop, kind='geom');
        zi     = np.ma.array( fzi(xi,yi) , mask=get_mask(xi, yi,zop) )

        cs = ax.contourf(xi, yi, zi,levels=100, cmap="bone")
#        fig.colorbar(cs, ax=ax, shrink=0.9 , label=r"$\rm Energy (meV)$")
        # This is the fix for the white lines between contour levels
        for c in cs.collections:
            c.set_edgecolor("face")


        xi, yi = xypoints(npts=20)
        fui= mtri.CubicTriInterpolator(triang, uop, kind='geom');
        fvi= mtri.CubicTriInterpolator(triang, vop, kind='geom');
        ui = np.ma.array(fui(xi, yi), mask=get_mask(xi, yi,zop));
        vi = np.ma.array(fvi(xi, yi), mask=get_mask(xi, yi,zop));
        Q = ax.quiver(xi, yi, ui, vi, color="C1",  scale=15, width=0.022/4);

        return ax;
    
    
    def createRandomPhase(self, dims, output, mask=None):
        def linearize( i , dims ):
            x,y,z =i ;
            nx,ny,nz =dims ;
            return z*nx*ny + y*nx + x;

        nx,ny,nz = dims;
        numstates  = matdim = len(self.ham_operator([0,0,0]))
        dim = np.prod(dims);
        states = np.zeros((numstates,numstates*dim), dtype=complex)
        grid =(np.mgrid[0:nx,0:ny,0:nz].T).reshape( dim ,3);

        for kp in grid:
            if mask is not None and mask(kp/dims):
                continue;
            else:
                k = linearize(kp,dims);
                kp = kp/dims;
                w, vs  = np.linalg.eigh(self.ham_operator(kp))
                vs = vs.T;
                for i,v in enumerate(vs):
                    states[i,k*matdim:(k+1)*matdim]= np.exp(2j*np.pi*np.random.random())*v;

        for s in range(numstates):
            for i in range(matdim):
                S = states[s][i::matdim].reshape(dims).T
                FS = np.fft.fftn(S)
                states[s][i::matdim] = FS.T.flatten();
            states[s] /=np.linalg.norm(states[s]);

        f = open(output,'w')
        f.write("vector\n")
        f.write( str(numstates)+"\n")
        f.write( str(dim*matdim)+"\n")

        for state in states:
            for x in state:
                f.write(str(x.real)+" "+str(x.imag)+"\n")
        f.close()


