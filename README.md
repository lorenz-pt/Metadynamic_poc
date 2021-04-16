# Metadynamics_poc

Implementation of Metadynamics on a simple topological system (A particle confined on a circle).

The program uses an HMC (Hamiltonian Monte Carlo or Hybrid Monte carlo) algorithm to simulate the paths distribution with a probability density given by the action functional. It makes use of the second order Omelyan symplectic integrator.

Metadynamics is needed for healing topological freezing at low temperature (in the code, this is equivalent to choosing a small 'e' or great values of 'N', the length of the paths).

=============================================================================================================================
MATLAB VERSION: (meta_omelyan.m) (it needs to be updated)
=============================================================================================================================

The Matlab version of the code, with N = 300, T = 156, takes about 1:40 minutes.

=============================================================================================================================
FORTRAN VERSION: (metadi.f)
=============================================================================================================================
A new Fortran version has been added. The code is the transposition of the one implemented in Matlab. 

The code require as an external input (file 1 = 'meta_input'):

-measures:number of measures to be taken,
-i_time: integration time ,
-dt: time step for the integrator,
-par: parameter of the action beta* hbar* chi (in the continuum <Q^2> = beta * hbar* chi),
-Qtrh: threshold value of the charge,
-dq: the charge spacing for Metadynamics,
-hgt: heigth of the time dependent potential,
-strength of the external force for |Q|> Qtrh,
-start: if 0 the field is initially set to 0, if 1, the field values are chosen casually between 0-1.

PARAMETERS TO SET DIRECTLY INSIDE THE CODE:

the mean of the time dependent potential is taken starting by a certain step (specified at line 68).
the time length of the path is 'nlattice' 
the length of the time dependent potential is 'length' which must include also the values at '-Qtrh - dq' and 'Qtrh + dq' (these values then, aren't saved as an output)

OUTPUT:

the code saves the values of the charge on an external file (file 2 =  'meta_misure_#').
the code saves the values of the mean of the time dependent potential on (file 3 = 'tdp_#')
the code prints as an output the value of the acceptance rate of the paths.


 
