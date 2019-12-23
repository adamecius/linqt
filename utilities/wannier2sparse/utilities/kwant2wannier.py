import kwant
import numpy as np
import datetime

def syst2txt(syst, systname ,params ):
    filename= systname+"_hr.dat";
    with open(filename, 'w') as out:
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

def geo2txt(syst,lattice_vectors,site_dict ,systname ):
    filename= systname+".uc";
    with open(filename, 'w') as out:
        for vec in lattice_vectors:
            for x in vec:
                if( x>= 0):
                    out.write(" ")
                out.write("{:.8E} ".format(x))
            out.write("\n");

    filename= systname+".xyz";
    orb_pos = [ tuple(site.pos) for site in syst.sites() ] ;
    with open(filename, 'w') as out:
        out.write( " %d \n"%(2*len(orb_pos)))
        for spin in range(2): #GENERALIZE 
            for orb in orb_pos:
#                orb = tuple( np.dot( orb, lattice_vectors) );
                out.write( site_dict[spin]+" ");
                for x in orb:
                    if( x>= 0):
                        out.write(" ")
                    out.write("{:.8E} ".format(x))
                out.write("\n")


def convert2wannier(syst,params):

    #In Kwant, the hoppings and onsites can be functions
    #depending on some set of parameters.
    #This function take these values and evaluate them using
    #the params passed from the main.
    def evalute_value(obj,value,params):
        return value(obj,**params);

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
            matrix_value = evalute_value(site, matrix_value,params);
            for spin_i,spin_j, value in  get_matrix_values(matrix_value) :
                index_i =combine_orbital_inner_indexes( orb_idx=orbIndexes[site.family], numOrbs=numOrbs, spin_idx=spin_i );
                index_j =combine_orbital_inner_indexes( orb_idx=orbIndexes[site.family], numOrbs=numOrbs, spin_idx=spin_j );
                if( np.abs(value) > 1e-13 ):# Remove very small non zero values
                    hop = (*site.tag,index_i+1,index_j+1, np.real(value), np.imag(value));
                    onsites.append(hop);
                    h_hop = conj(*hop);
                    if ( hop != h_hop):
                        onsites.append(h_hop );
        return onsites;
    
    #Get the hoppings values and lattice information and storage in a list
    def getHoppings( orbIndexes, hopping_value_pairs ):
        hoppings = list();
        numOrbs = len( orbIndexes );
        for (site_i, site_j), matrix_value in hopping_value_pairs:
            matrix_value = evalute_value((site_i, site_j), matrix_value,params);
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
