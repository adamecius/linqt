import kwant
import tinyarray as ta
import numpy as np
from numpy import sqrt
from numpy import linalg as LA
from kwant import wraparound
from kwant.builder import HoppingKind as hop

#Spin matrices
S0 = ta.array([[1,  0 ], [0 ,  1]])
Sx = ta.array([[0,  1 ], [1 ,  0]])
Sy = ta.array([[0, -1j], [1j,  0]])
Sz = ta.array([[1,  0 ], [0 , -1]])


#Define graphene lattice 
lat = 1.0; lattice_vectors = lat*ta.array([(sqrt(3)/2,1/2,0), (sqrt(3)/2,-1/2,0),(0,0,1)] ); 
lat = kwant.lattice.general(lattice_vectors, [(-1/3,-1/3,0),(0,0,0)]) 
a, b = lat.sublattices

#Dirac Model Operator
def make_bulk_graphene():
    
    syst = kwant.Builder( kwant.TranslationalSymmetry(lat.vec((1,0,0)), lat.vec((0,1,0)), lat.vec((0,0,1))) ) 

    
    def onsite_potential(site, x ):
        return S0*x;

    def sublattice_potential(site, x  ):
        return S0*(x,-x)[ site.family == a];

    def nearest_hopping(hop, x ):
        return S0*x;

    def rashba_soc(hop, x  ):
        site_i,site_j = hop;
        dij = site_i.pos-site_j.pos;
        dij = dij/LA.norm(dij);
        geo_fact = (2j/3)*(Sx*dij[1] - Sy*dij[0])
        return x*geo_fact;

    def intrinsic_soc(hop, x  ):
        site_i,site_j = hop;
        geo_fact = 1.0j/3.0/np.sqrt(3.);
        lambda_soc=x*Sz*geo_fact;
        return (lambda_soc,-lambda_soc)[ site_j.family == a];

    def valleyZeeman_soc(hop, x ):
        site_i,site_j = hop;
        geo_fact = 1.0j/3.0/np.sqrt(3.);
        lambda_soc=x*Sz*geo_fact;
        return lambda_soc;


    tb_onsite_dict = {
                        "fermi_energy":onsite_potential,\
                        "staggered_potential":sublattice_potential,\
                        };


    tb_nearest_dict = {
                        "hopping_t0":nearest_hopping,\
                        "lambda_rso":rashba_soc,\
                        };

    tb_next_nearest_dict = {
                        "lambda_iso":intrinsic_soc,\
                        "lambda_vz" :valleyZeeman_soc,\
                        };

    
    
    
    
    #INCORPORATE ONSITES
    def sum_functions( x , **kwargs ):
        total = 0;
        if kwargs is not None:
            for key, value in kwargs.items():
                if key  in tb_onsite_dict:
                    total+= (tb_onsite_dict[key])(x,value);
        return  total; 
    for site in lat.sublattices:
        syst[site(0,0,0)] = sum_functions

    #INCORPORATE NEAREST NEIGHOORS
    def sum_functions( x , **kwargs ):
        total = 0;
        if kwargs is not None:
            for key, value in kwargs.items():
                if key  in tb_nearest_dict:
                    total+= (tb_nearest_dict[key])(x,value);
        return  total; 
    nn_hoppings = (((0,0,0), b, a), ((-1,0,0), b, a), ((0,-1,0), b, a)) #in kwant b,a means b=target, a=origin, if using positive lat_vec
    for hopping in nn_hoppings:
        syst[ hop(*hopping) ] = sum_functions;

#    #INCORPORATE NEXT NEAREST NEIGHOORS   
    def sum_functions( x , **kwargs ):
        total = 0;
        if kwargs is not None:
            for key, value in kwargs.items():
                if key  in tb_next_nearest_dict:
                    total+= (tb_next_nearest_dict[key])(x,value);
        return  total; 
    nnn_hoppings = (((-1,0,0), a, a),((0,1,0), a, a),((1,-1,0), a, a), 
                    ((-1,0,0), b, b),((0,1,0), b, b),((1,-1,0), b, b) )
    for hopping in nnn_hoppings:
        syst[ hop(*hopping) ] = sum_functions;   





    
    return syst
