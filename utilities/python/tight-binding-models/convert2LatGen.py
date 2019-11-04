import kwant
import numpy as np


def syst2txt(syst, filename ,params ):
    with open(filename, 'w') as out:
        for elem in convert2LatGen(syst,params):
            out.write("%d %d %d %d %d %f %f \n"% (elem));

def geo2txt(syst,lattice_vectors,filename ):
    orb_pos = [ tuple(site.pos) for site in syst.sites() ] ;
    with open(filename, 'w') as out:
        for vec in lattice_vectors:
            out.write("%f %f %f \n"% (tuple(vec)));
        for spin in range(2): #GENERALIZE 
            for orb in orb_pos:
                out.write("%f %f %f \n"% (orb));


            
def convert2LatGen(syst,params):
    
    #Evalute the functions contained inside value using params.
    #When there is no value simply return value #TODO
    def evalute_value(obj,value,params):
        return value(obj,**params);

    #Evalute the functions contained inside value using params.
    #When there is no value simply return value #TODO
    def conj(index_i,index_j,dx,dy,dz, re_val, im_val):
        return (index_j,index_i,-dx,-dy,-dz,re_val,-im_val);
    
                 
    #Returns a dictionary which assignates an ID to each inequivalent site
    #of a given kwant system.
    def getSitesID( siteList ):
        siteID_pairs = [ (site.family,io) for io,site in enumerate(siteList)] ;
        return dict(siteID_pairs)
        
    #Kwant systems can have a matrix as hoppings and onsites.
    #in this program we linearize this matrix, and assign
    #a unique index to the diagonal elements.
    #this function returns the index pair and value (row, col, val)
    #of the matrix value. 
    #If is not a matrix, it should returns simply (0,0,value); #TODO
    def get_matrix_values( matrix_value ):
        value_list = list()
        for row, row_value in enumerate(matrix_value):
            for col, value in enumerate(row_value):
                value_list.append((row,col,value));
        return value_list;
    
    #Combine the inner indexes( the indexes obtained from the matrix value)
    #with the site indexes obtained from getSItesID in a single 
    def combine_orbital_inner_indexes( orb_idx=0, numOrbs=1, spin_idx=0, numSpins=1 ):
        return orb_idx*numOrbs + spin_idx; 

    #Get the onsite values and lattice information and storage in a list
    def getOnsites( orbIndexes, site_value_pairs ):
        onsites = list();
        numOrbs = len( orbIndexes );
        for site, matrix_value in  site_value_pairs :
            matrix_value = evalute_value(site, matrix_value,params);
            for spin_i,spin_j, value in  get_matrix_values(matrix_value) :
                index_i =combine_orbital_inner_indexes( orb_idx=orbIndexes[site.family], numOrbs=numOrbs, spin_idx=spin_i );
                index_j =combine_orbital_inner_indexes( orb_idx=orbIndexes[site.family], numOrbs=numOrbs, spin_idx=spin_j );
                hop = (index_i,index_j, *site.tag, np.real(value), np.imag(value));
                onsites.append(hop);
                h_hop = conj(index_i,index_j, *site.tag, np.real(value), np.imag(value));
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
                hop = (index_i,index_j, *site_j.tag, np.real(value), np.imag(value));
                hoppings.append(hop);
                h_hop = conj(index_i,index_j, *site_j.tag, np.real(value), np.imag(value));
                if ( hop == h_hop):
                    hoppings.append(h_hop );
        return hoppings;

    orbIndexes = getSitesID(syst.sites() ); #Extract from syst.
    onsites  = getOnsites( orbIndexes, syst.site_value_pairs() );
    hoppings = getHoppings( orbIndexes,syst.hopping_value_pairs() );
    return onsites+hoppings;
