Wannier2Sparse is a package that allows for creating tight-binding operators from Wannier90 outputs. 

From Wannier90, we use:

(i) The wannier hamiltonian. A file created from Wannier90 which usually have the following sufixx _hr.dat, 
and that contains the Hamiltonian of the system in the following format:

Number of Wannier Functions
Number of Wigner-Seitz points
Degeneracy of each  Wigner-Seitz point by 15 elements per line
...
...
Finally, the remaining num_wann 2 × nrpts lines each contain, respectively, the components of the vector R in terms of the lattice vectors {A i }, 
the (R) indices m and n, and the real and imaginary parts of the Hamiltonian matrix element H mn in the WF
basis.

(ii) From the _.out we extract the position of the wannier functions and the unit cell. 

Example:
Linear Chain.
For the typical linear chain with nearest neighbor hoppings T, onsite potential E0, lattice constant a and a single orbital per unit cell, the linear_chain_hr.dat will look as
1
1
1
 0 0 0 1 1 E0 0.0
 1 0 0 1 1 Real(T) Imag(T)
-1 0 0 1 1 Real(T) Imag(T)

while the .uc file will look as
a 0.0 0.0 
0.0 1.0 0.0 
0.0 0.0 1.0

Finally, if we assume that the orbital is located at a position r0, the *.xyz file has this format
r0x r0y r0z 
