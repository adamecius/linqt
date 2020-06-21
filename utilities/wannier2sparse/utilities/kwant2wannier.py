import kwant
import numpy as np
import datetime
import inspect


from  kwant.builder import Site as kwantSite

def evalute_value(obj,value,params):
    if( callable(value) ):
        #Extract the signature and determine the positional arguments
        pos_params= map(str.strip,str(inspect.signature(value)).strip("()").split(","));
        pos_params= [ params[par] for par in pos_params if( par in params ) ]
        
        #If is an onsite, cast the onsite
        if( isinstance( obj, kwantSite ) ):
            return value(obj,*pos_params);
        
        #If is a hopping(tuple of sites), cast the hopping
        if( len(obj)==2 and all(isinstance(x, kwantSite) for x in obj ) ):
            return value(*obj,*pos_params);
        
        I = np.eye(2);
        #If value is not callable we first assume is scalar
        #and return its value times the identity
        if( isinstance(value, (float,complex) ) ):
            return value*I;

    #If is not an scalar, then we assume is a matrix
    #and check its dimension
    try:
        if ( len(value) == 2 ):
            return np.array( value);
    except:
        return I*0;

    print( "Return 0, Invalid object:", obj, " is invalid" )
    return I*0;


def fold_slab( syst, lat_vec, params):

    #Get the system's translational symmetry
    kTS = syst.symmetry;
    
    #Define the sites within the fundamental domain as the new orbitals
    orbs = [ site.pos for site in syst.sites() ];

    #Create an operator that will remove the tag component outside the fundamental domain.
    #There should be a more elegant way of doing thi
    sym_dir  =  np.dot(kTS.periods , np.linalg.inv(lat_vec));
    fold_tag = np.eye(3)*0;
    for idx in np.argmax(sym_dir,axis=1):
        fold_tag[idx,idx]=1;

    
    #The lattice vectors from the system is needed too for creating the new lattice
    names = ["kwantSite_"+str(i) for i,orb in enumerate(orbs) ];
    lat = kwant.lattice.general( prim_vecs= lat_vec, \
                                 basis    = orbs, \
                                 norbs    = len(orbs),\
                                 name     = names);

    #Create a map of sites from the old kwant system to a family
    fam = lat.sublattices;
    fam = { site:fam[i] for i,site in enumerate(syst.sites()) };
    
    #Create new kwant system with the previous information
    nsyst = kwant.Builder( kwant.TranslationalSymmetry(  lat.vec((1,0,0)), lat.vec((0,1,0)) ) ) 


    #Map the sites from old kwant system to families in the new.
    #Evaluate the functions because they depend on the old sites
    for site,val in syst.site_value_pairs():
        nsyst[ fam[site](0,0,0) ] = evalute_value(site,val,params);
    
    for hop,val in syst.hopping_value_pairs():
        fi,fj  = [ fam[ kTS.to_fd(site) ](*fold_tag.dot(site.tag) ) for site in hop ];
        nsyst[ fi,fj ] = evalute_value( hop ,val,params);
    
    return nsyst,lat;


def syst2linQT( syst, systname, lattice_vectors, params=dict() ):

    syst,lat = fold_slab( syst, lattice_vectors,params)
    lattice2txt( lat , systname )    
    
    with open( systname+"_hr.dat" , 'w') as out:
        now = datetime.datetime.now()
        out.write(" written on ");
        out.write(now.strftime("%Y-%m-%d %H:%M:%S \n"))

        elem_list = convert2wannier(syst,params);
        maxOrbId =np.max([ np.max([elem[3],elem[4]]) for elem in elem_list ] );

        out.write("        %d\n"% maxOrbId ) 
        out.write("        %d\n"% 1) 
        out.write("    %d\n"% 1) 

        for elem in elem_list:
            out.write("%d %d %d %d %d %f %f \n"% (elem));
           
           
def lattice2txt( lattice , systname ):

    #Save unit cell
    np.savetxt( systname+".uc", lattice.prim_vecs )
    
    #Save orbital positions in xyz format
    orbs = lattice.sublattices;
    with open( systname+".xyz" , 'w') as out:
        out.write( " {} ".format( 2*len(orbs) )+" \n")
        for spin in ("_up_","_down_"):
            for site in orbs:
                pos = tuple(site(0,0,0).pos);
                out.write(site.name+spin+" {:.8E} {:.8E} {:.8E}".format(*pos)+" \n")
    return True;


def convert2wannier(syst,params):

    #Kwant only storage hoppings without their conjugated complex.
    #if both are submitted only the first one is storage.
    #Wannier require you to give it both.
    #This function compute the conjugate gradient of
    #a hopping.
    def conj(dx,dy,dz,index_i,index_j, re_val, im_val):
        return (-dx,-dy,-dz,index_j,index_i,re_val,-im_val);


    #Returns a dictionary which assignates an ID to each inequivalent site
    #of a given kwant system.
    def getSitesID( siteList ):
        siteID_pairs = [ (site.family,io) for io,site in enumerate(siteList)] ;
        return dict(siteID_pairs)

    #Kwant allows for the hoppings and onsites of a hamiltonian
    #to be defined as a matrix. This function transform the matrix to coo format 
    #better handling. If the matrix is an scalar, return (0,0,value)
    def get_matrix_values( matrix_value ):
        value_list = list()
        for row, row_value in enumerate(matrix_value):
            for col, value in enumerate(row_value):
                value_list.append((row,col,value));
        return value_list;


    #Kwant allows for the hoppings and onsites of a hamiltonian
    #to be defined as a matrix.  Wannier systems require complex
    # hoppings and values.To fix this, we linearize the (i,j)
    #indices of the matrix in a unique way (i,j)->k.
    #These indexes will be used by the function
    #combine_orbital_inner_indexes() for creating an unique
    #site ID.
    #Combine the inner indexes( the indexes obtained from the matrix value)
    #with the site indexes obtained from getSItesID in a single
    def combine_orbital_inner_indexes( orb_idx=0, numOrbs=1, spin_idx=0, numSpins=2 ):
        return spin_idx*numOrbs + orb_idx;


    #Get the onsite values and lattice information and storage in a list
    def getOnsites( orbIndexes, site_value_pairs ):
        onsites = list();
        numOrbs = len( orbIndexes );
        for site, matrix_value in  site_value_pairs :
            for spin_i,spin_j, value in  get_matrix_values(matrix_value) :               
                index_i =combine_orbital_inner_indexes( orb_idx=orbIndexes[site.family], numOrbs=numOrbs, spin_idx=spin_i );
                index_j =combine_orbital_inner_indexes( orb_idx=orbIndexes[site.family], numOrbs=numOrbs, spin_idx=spin_j );
                upper_matrix = bool(spin_j >= spin_i ); #Consider the upper part of the matrix
                if( np.abs(value) > 1e-13 and upper_matrix ):# Remove very small non zero values 
                    hop = (*site.tag,index_i+1,index_j+1, np.real(value), np.imag(value));
                    onsites.append(hop);
                    h_hop = conj(*hop);
                    if ( hop != h_hop  ):
                        onsites.append(h_hop );
        return onsites;
    
    #Get the hoppings values and lattice information and storage in a list
    def getHoppings( orbIndexes, hopping_value_pairs ):
        hoppings = list();
        numOrbs = len( orbIndexes );
        for (site_i, site_j), matrix_value in hopping_value_pairs:
            for spin_i,spin_j, value in  get_matrix_values(matrix_value) :
                index_i =combine_orbital_inner_indexes( orb_idx=orbIndexes[site_i.family], numOrbs=numOrbs, spin_idx=spin_i );
                index_j =combine_orbital_inner_indexes( orb_idx=orbIndexes[site_j.family], numOrbs=numOrbs, spin_idx=spin_j );
                if( np.abs(value) > 1e-13 ): 
                    hop = (*site_j.tag, index_i+1,index_j+1, np.real(value), np.imag(value));
                    hoppings.append(hop);
                    h_hop = conj(*hop);
                    if ( hop != h_hop):
                        hoppings.append(h_hop );

        return hoppings;

    orbIndexes = getSitesID(syst.sites() ); #Extract from syst.
    onsites  = getOnsites( orbIndexes, syst.site_value_pairs() );
    hoppings = getHoppings( orbIndexes,syst.hopping_value_pairs() );
    return onsites+hoppings;
