import kwant
import tinyarray as ta
import numpy as np
from numpy import sqrt
from numpy import linalg as LA
from kwant import wraparound
from kwant.builder import HoppingKind as hop

#Constant values
hbar = 0.65821;

#Spin matrices
S0 = ta.array([[1,  0 ], [0 ,  1]])
Sx = ta.array([[0,  1 ], [1 ,  0]])
Sy = ta.array([[0, -1j], [1j,  0]])
Sz = ta.array([[1,  0 ], [0 , -1]])


def get_lattice_vectors( lattice_constant = 1.0):
    #Define graphene lattice 
    lat = lattice_constant; 
    return lat*ta.array([(sqrt(3)/2,1/2,0), (sqrt(3)/2,-1/2,0),(0,0,1)] ); 

def get_lattice( lattice_constant = 1.0):
    #Define graphene lattice 
    lattice_vectors = get_lattice_vectors( lattice_constant ); 
    return kwant.lattice.general(lattice_vectors, [(-lattice_constant/sqrt(3),0,0),(0,0,0)])#site convetion opossite to my C program

#Tight Binding Model Operator
def make_bulk_graphene(lat_const = 1.0):

    lat = get_lattice( lattice_constant = lat_const);
    a, b = lat.sublattices
    syst = kwant.Builder( kwant.TranslationalSymmetry(lat.vec((1,0,0)), lat.vec((0,1,0)), lat.vec((0,0,1))) ) 

    def onsite_potential(site, x ):
        return S0*x;

    def sublattice_potential(site, x  ):
        return S0*(x,-x)[ site.family == a];

    def ferro_exchange(site, x ):
        return -( Sx*x[0] + Sy*x[1] + Sz*x[2] );

    def antif_exchange(site, x ):
        potential = -(Sx*x[0]+Sy*x[1]+Sz*x[2]); 
        return (potential,-potential)[ site.family == a];

    def valley_exchange(hop, x ):
        geo_fact = 1.0j/3.0/np.sqrt(3.);
        lambda_soc=x*geo_fact;
        return S0*lambda_soc;

    def nearest_hopping(hop, x ):
        return S0*x;

    def rashba_soc(hop, x  ):
        site_i,site_j = hop;   #This has been checked
        dij = site_i.pos-site_j.pos; #This has been checked
        dij = dij/LA.norm(dij);
        geo_fact = (2.0j/3.0)*(Sx*dij[1] - Sy*dij[0])
        return x*geo_fact;

    def valleyZeeman_soc(hop, x ):
        lambda_soc=np.dot(Sz,valley_exchange(hop,x));
        return lambda_soc;

    def intrinsic_soc(hop, x  ):
        site_i,site_j = hop;
        lambda_soc=valleyZeeman_soc(hop, x );
        return (lambda_soc,-lambda_soc)[ site_j.family == a];
    
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
    nn_hoppings = (((0,0,0), b, a), ((-1,0,0), b, a), ((0,-1,0), b, a)) #In kwant b,a means b=target, a=origin(checked)
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
    nnn_hoppings = ((( 1,0,0), a, a),((0,-1,0), a, a),((-1, 1,0), a, a), 
                    (( 1,0,0), b, b),((0,-1,0), b, b),((-1, 1,0), b, b) )#The vector is R_i-R_j (checked)
    for hopping in nnn_hoppings:
        syst[ hop(*hopping) ] = sum_functions;   
    
    return syst

#Tight Binding Model Operator
def dirac_model(kp , hop_list , lat_const = 1.0):
    #Pauli matrices
    s0 = np.array([[1,0  ],[0 ,1 ]],dtype=complex);    
    sx = np.array([[0,1  ],[1 ,0 ]],dtype=complex); 
    sy = np.array([[0,-1j],[1j,0 ]],dtype=complex); 
    sz = np.array([[1,0  ],[0 ,-1]],dtype=complex); 
  
    #Add additional spin-pseudospin
    SigX= np.kron( s0, sx ); SigY= np.kron( s0, sy ); SigZ= np.kron( s0, sz );
    SX  = np.kron( sx, s0 ); SY  = np.kron( sy, s0 ); SZ  = np.kron( sz, s0 ); 
    I0  = np.kron( s0, s0 );
    
    #Add additional valley degree of freedom
    SigX= np.kron( s0  , SigX); SigY= np.kron( s0  , SigY); SigZ= np.kron( s0  , SigZ);
    SX  = np.kron( s0  , SX  ); SY  = np.kron( s0  , SY  ); SZ  = np.kron( s0  , SZ  ); 
    TX  = np.kron( sx  , I0  ); TY  = np.kron( sy  , I0  ); TZ  = np.kron( sz  , I0   ); 
    I0  = np.kron( I0, s0 );

    # Order: ( tau , spin , pseudospin )

    #Orbital functions
    def onsite_potential( kp , x ):
        return I0*x;

    def sublattice_potential( kp , x  ):
        return SigZ*x;

    def nearest_hopping( kp , x ):
        vF = ( np.sqrt(3)*lat_const*x/hbar/2 );
        return hbar*vF*( np.dot(TZ,SigX)*kp[0] + SigY*kp[1] ); 

    #Magnetic functions
    def ferro_exchange( kp , x ):
        return -( SX*x[0] + SY*x[1] + SZ*x[2] );

    def antif_exchange( kp , x ):
        return -np.dot( ( SX*x[0]+SY*x[1]+SZ*x[2] ), SigZ ); 
                
    #Spin-orbit coupling functions
    def rashba_soc(kp, x  ):
        return x * ( np.dot( np.dot(TZ,SigX) ,SY) - np.dot(SigY,SX) );

    def valley_exchange( kp , x ):
        return x*TZ;

    def valleyZeeman_soc( kp , x ):
        return np.dot( valley_exchange( kp , x ), SZ );

    def intrinsic_soc( kp , x  ):
        return np.dot( valleyZeeman_soc( kp, x ), SigZ);

    model_dict = {
                        "fermi_energy":onsite_potential,\
                        "staggered_potential":sublattice_potential,\
                        "ferro_exchange":ferro_exchange,\
                        "antif_exchange" :antif_exchange,\
                        "hopping":nearest_hopping,\
                        "lambda_rashba":rashba_soc,\
                        "lambda_intrinsic":intrinsic_soc,\
                        "lambda_valleyZeeman" :valleyZeeman_soc,\
                        "lambda_valleyExchange" :valley_exchange
                  };
    
    def Hamiltonian( kp , hop_list ):
        total = onsite_potential( 0 , 0 );
        for key  in hop_list:
            if key  in model_dict:
                total+= (model_dict[key])(kp,hop_list[key]) ;
        return total;

    return Hamiltonian( kp , hop_list )
