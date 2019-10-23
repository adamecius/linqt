import numpy as np
from scipy.sparse import coo_matrix

BOUND = 0.95;

def chebG(m, z ):
    fz = ( 1.0 - z**2 )**(1/2)
    return -1j*( z - 1.0j*fz)**m/fz

def spectral_bounds( A, tol=0.1 ):
    from scipy.sparse.linalg import eigsh, svds
    from numpy.linalg import eigvals
 
    dim= max( A.shape )
    if( dim < 1000 ):
        allE=eigvals( coo_matrix(A).todense())
        Emax= max(allE);
        Emin= min(allE);
    else:
        Emax =eigsh( coo_matrix(A), k=1, which='LA',return_eigenvectors=False,maxiter=100, tol=tol)[0]
        Emin =eigsh( coo_matrix(A), k=1, which='SA',return_eigenvectors=False,maxiter=100, tol=tol)[0]

    return Emin,Emax;


def scale_and_shiftOp(A, a,b ):
    A= coo_matrix(A);
    A.setdiag( A.diagonal() + b/a);
    A.data*=a;
    return A;


def scale_spectrum(A, bounds):
    bmin,bmax=bounds; 
    Emin,Emax=spectral_bounds( A );
    scale = (bmax-bmin)/ (Emax-Emin)
    shift = (bmin*Emax-bmax*Emin)/ (Emax-Emin);
    scale_and_shiftOp(  A, scale, shift );   
    return A;

def scale_and_shift( A ):
    Emin,Emax=spectral_bounds( A );
    W  =( Emax-Emin )*0.5;
    Ec =( Emax+Emin )*0.5;
    return W,Ec


def compute_dos_moms(numMom, Hmat, vec, scalpar=0):
    vec = vec/ np.linalg.norm(vec);
    if ( scalpar == 0 ):
        Hmat = scale_spectrum(Hmat, (-BOUND,BOUND)); #only works for sparse matrices
    else:
        W = scalpar[0]; #Set the energy scale;
        Ec= scalpar[1]; #Set the energy shift;
        Hmat = scale_and_shiftOp(Hmat,1.0/W,-Ec/W); #only works for sparse matrices
    #Compute the chebyshev moments
    mu = np.zeros( numMom, dtype=complex)
    chebT0 = vec;
    chebT1 = Hmat.dot(chebT0);
    for m in range(numMom):
        mu[m]  = np.vdot( vec,  chebT0 )
        chebT0 = 2.0*Hmat.dot(chebT1) - chebT0;
        chebTmp=chebT0; chebT0=chebT1; chebT1=chebTmp;
    mu*=2.0;
    mu[0]*=0.5;
    return mu

def dos(mu ,x, y=0 ): 
    from numpy import cos, sqrt, arccos
    m_arr = range(len(mu) );
    GE_arr =chebG(m_arr,x + 1j*y ) ;    
    return -np.imag( mu.dot(GE_arr) )/np.pi;


def compute_neq_2Dmoms(numMom, Hmat, O1mat, O2mat,  vec, scalpar=0):
    vec = vec/ np.linalg.norm(vec);
    if ( scalpar == 0 ):
        Hmat = scale_spectrum(Hmat, (-BOUND,BOUND)); #only works for sparse matrices
    else:
        W = scalpar[0]; #Set the energy scale;
        Ec= scalpar[1]; #Set the energy shift;
        Hmat = scale_and_shiftOp(Hmat,1.0/W,-Ec/W); #only works for sparse matrices
              
    mu = np.zeros( (numMom,numMom), dtype=complex);

    cheb0T0 = O1mat.dot(vec);
    cheb0T1 = Hmat.dot(cheb0T0);
    for m0 in range(numMom):
        cheb1T0 = vec;
        cheb1T1 = Hmat.dot(cheb1T0);
        for m1 in range(numMom):
            mu[m0,m1]  = np.vdot( O2mat.dot(cheb1T0) , cheb0T0 )
            cheb1T0 = 2.0*Hmat.dot(cheb1T1) - cheb1T0;
            chebTmp=cheb1T0; cheb1T0=cheb1T1; cheb1T1=chebTmp;
        cheb0T0 = 2.0*Hmat.dot(cheb0T1) - cheb0T0;
        chebTmp=cheb0T0; cheb0T0=cheb0T1; cheb0T1=chebTmp;
    mu *= 4.0; 
    mu[0,:] *= 0.5; 
    mu[:,0] *= 0.5; 
    return mu 


def save2DChebMom(dim, scalpar, mu2D, filename):
    numMom = mu2D.shape 
    W, Em = scalpar
    f= open(filename,"w")
    f.write("%d %f %f\r\n" % (dim,2.0*W,Em) );
    f.write("%d %d\r\n" % numMom);
    for m0 in range(numMom[0]) :
        for m1 in range(numMom[1]) :
            z = mu2D[m0,m1]; 
            f.write("%.12E %.12E\r\n" % (np.real(z),np.imag(z)));
    print( "finished saving the moments in", filename )
    f.close();
    
    
    
    
    
    
