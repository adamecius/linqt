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
lat_const     = 1.0; 
unit_cell_vec = lat_const*np.array([(sqrt(3)/2,1/2,0), (sqrt(3)/2,-1/2,0), (0,0,1)] ); 

lattice_vectors = ta.dot([[1,1,0],[-1,2,0], [0,0,1]],unit_cell_vec);
orbs = ta.dot( [[0,0,0],[-1/3,-1/3,0],[0,-1,0],[2/3,-4/3,0],[1,-1,0],[2/3,-1/3,0]], unit_cell_vec)

print(lattice_vectors)

lat   = kwant.lattice.general(lattice_vectors, orbs )#site convetion opossite to my C program
a0, b0, a1, b1,a2,b2  = lat.sublattices
a_family = [a0,a1,a2];
b_family = [b0,b1,b2];

#FUNCTION NECESSARY TO IDENTIFY THE CHIRALITY OF THE INTRINSIC
def intrinsic_sign(x):
    TOL=0.99
    ucv = unit_cell_vec;
    int_vecs = [ucv[0],-ucv[1],-ucv[0]+ucv[1]]
    for vec in int_vecs:
        result = np.dot(x,vec) 
        if( result >= TOL ):
            return 1;
    return -1;
        

#Dirac Model Operator
def make_bulk_graphene():
    
    syst = kwant.Builder( kwant.TranslationalSymmetry(lat.vec((1,0,0)), lat.vec((0,1,0))) ) 

    
    
    def onsite_potential(site, x ):
        return S0*x;

    def sublattice_potential(site, x  ):
        return S0*(x,-x)[ site.family in a_family ];

    def ferro_exchange(site, x ):
        return ( Sx*x[0]+Sy*x[1]+Sz*x[2]);

    def antif_exchange(site, x ):
        potential = (   Sx*x[0]+Sy*x[1]+Sz*x[2]  ); 
        return (potential,-potential)[ site.family in a_family];

    def valley_exchange(hop, x ):
        site_i,site_j = hop;
        dij = site_i.pos-site_j.pos;
        geo_fact = 1.0j/3.0/np.sqrt(3.)*intrinsic_sign(dij);
        lambda_soc=x*S0*geo_fact;
        return lambda_soc;
                
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
        lambda_soc=np.dot(Sz,valley_exchange(hop, x ));
        return (lambda_soc,-lambda_soc)[ site_j.family in a_family];

    def valleyZeeman_soc(hop, x ):
        lambda_soc=np.dot(Sz,valley_exchange(hop, x ) );
        return lambda_soc;

    
    
    tb_onsite_dict = {
                        "fermi_energy":onsite_potential,\
                        "staggered_potential":sublattice_potential,\
                        "ferro_exchange":ferro_exchange,\
                        "antif_exchange" :antif_exchange
                        };


    tb_nearest_dict = {
                        "hopping":nearest_hopping,\
                        "lambda_rashba":rashba_soc,\
                        };

    tb_next_nearest_dict = {
                        "lambda_intrinsic":intrinsic_soc,\
                        "lambda_valleyZeeman" :valleyZeeman_soc,\
                        "lambda_valleyExchange" :valley_exchange
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
    syst[ lat.neighbors(1) ] = sum_functions;
    
    #INCORPORATE NEXT NEAREST NEIGHOORS   
    def sum_functions( x , **kwargs ):
        total = 0;
        if kwargs is not None:
            for key, value in kwargs.items():
                if key  in tb_next_nearest_dict:
                    total+= (tb_next_nearest_dict[key])(x,value);
        return  total;
    syst[ lat.neighbors(2) ]= sum_functions;
    

    
    return syst
