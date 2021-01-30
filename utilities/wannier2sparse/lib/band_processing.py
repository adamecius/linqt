#!/usr/bin/env python3
import kwant
import sisl
import numpy as np
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import json
import pickle
import wannier_tools

#This function creates a 2D mesh grid based on a kpoint windows kwindow=[kmin,kmax]
# kmin is the position of the origin of the windows while kmax is the position of the
# end of the rectangle
#    ________ kmax
#   |        |
#   |        |
#   |        |
#   |________|
# kmin
def create_2Dkgrid( kwindow=None, npoints=(1,1) ):
    if kwindow is None:
        kwindow = [ [0.0,0.0],[1.0,1.0]]
    kwindow= np.transpose(kwindow);
    kgrid = [ (np.linspace(kmin,kmax,p))  for (kmin,kmax),p in zip(kwindow,npoints) ];  
    return np.meshgrid(*kgrid,indexing='ij') ;

#Transform the kpoint-grid into a list of kpoints
def create_kpoints( kgrid ):
    K1,K2 = [ K.flatten() for K in kgrid ];
    return np.transpose([K1,K2,0*K2.flatten()]);

#Get the kpoints from the band path
def getBandKpoints(bandpath):
    klabels, kpoints, npoints = zip(*bandpath);
    npoints = np.array(npoints)[1:];
    kpoints = np.array(kpoints);
    kpath = np.swapaxes([kpoints[:-1],kpoints[1:]],0,1);
    kpath = np.concatenate([ np.linspace(kmin,kmax,npts,endpoint=False) for npts,(kmin,kmax) in zip(npoints,kpath) ],axis=0);
    kpath = list(kpath)+[kpoints[-1]];
    return kpath;

#Compute the X-axis of a band structure by summing the norm of the k steps
def getBandXaxis(kpath):
    return np.cumsum( np.linalg.norm(np.diff(kpath,axis=0,prepend=[kpath[0]]),axis=1) );

#Get the spin-projected bands using F= sisl, kwant
def path2kpoints(bandpath):
    klabels, kpoints, npoints = zip(*bandpath);
    npoints = np.array(npoints, dtype=int)
    kpoints = getBandKpoints(bandpath);
    return kpoints;

#select the bands in the (Emin,Emax) energy window
def spin_text_inwindow(text, Emin,Emax):
    EC = (Emin+Emax)/2; DE =(Emax-Emin)/2 
    return text[np.any(np.abs(text[:,0,:]-EC) < DE,axis=1)];

def energies_inwindow(bands, energy_window):
    Emin,Emax = energy_window;
    return np.any(bands>Emin,axis=1)*np.any(bands<Emax,axis=1);

def select_bands(bands, energy_window):
    return bands[energies_inwindow(bands,energy_window)];


#Compute the spin texture and bands of a Hamiltonian operator H with spin operators 
#given by S=Sx,Sy,Sz
def numpy_spin_texture(H, S):
    egval , eigvec = np.linalg.eigh(H);
    eigvec = eigvec.T;
    OP = ( H,*S);
    return np.array( [ [ np.real(np.conjugate(ei).dot(op.dot(ei))) for  ei in eigvec ] for op in OP] );


def getTextureFromWannierTools(kpoints, ws, S ):
    bands = [ws.compute_dispersion(kpoints)];
    spin_texture= bands+[np.swapaxes(ws.compute_dispersion(kpoints,proj_op=s),0,1)[1] for s in S ];
    return np.swapaxes( np.transpose(spin_texture),1,2);

def getTextureFromKwant(kpoints, syst, params, S ):
    spin_texture = [];
    for k in kpoints:
        k_x,k_y,k_z= 2*np.pi*k;
        params.update( {"k_x":k_x,"k_y":k_y,"k_z":k_z} );
        H = syst.hamiltonian_submatrix(params=params);
        spin_texture.append( numpy_spin_texture(H,S) );
    return np.swapaxes(spin_texture,2,0);

def getTextureFromSisl(kpoints, sH ):
    bz = sisl.BrillouinZone(sH, kpoints);
    def callback_bz(es, parent, k, weight):
        a = np.concatenate((es.eig.reshape(-1, 1), es.spin_moment()), axis=-1,);
        return a
    return np.swapaxes(bz.apply.array.eigenstate(wrap=callback_bz, eta=False).T,0,1);

