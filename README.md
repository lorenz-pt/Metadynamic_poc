# Metadynamics_poc

Implementation of Metadynamics on a simple topological system (A particle confined on a circle).

The program uses an HMC (Hamiltonian Monte Carlo or Hybrid Monte carlo) algorithm to simulate the paths distribution with a probability density given by the action functional. It makes use of the second order Omelyan symplectic integrator.

Metadynamics is needed for healing topological freezing at low temperature (in the code, this is equivalent to choosing a small 'e' or great values of 'N', the length of the paths).

=====================================================================================================================
MATLAB VERSION: (meta_omelyan.m) (it needs to be updated)
=====================================================================================================================

The Matlab version of the code, with N = 300, T = 156, takes about 1:40 minutes.

=====================================================================================================================
FORTRAN VERSION: (metadi.f)
=====================================================================================================================
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
-save_value: the sweep after which the code start to take the mean of the tdp
-tau_: after how many sweeps the tdp is updated

PARAMETERS TO SET DIRECTLY INSIDE THE CODE:

the length of the time dependent potential is 'length' which must include also the values at '-Qtrh - dq' and 'Qtrh + dq' (these values then, aren't saved as an output)

OUTPUT:

the code saves the values of the charge on an external file (file 2 =  'meta_misure_#').
the code saves the values of the mean of the time dependent potential on (file 3 = 'tdp_#')
the code prints the value of the acceptance rate of the paths on the terminal.

SUBROUTINES:

***

initialize_lattice: set the intial values of the path at all the time. if start = 0, every point is set to 0, while if start = 1, the points are chosen casually between -1 and 1. It is called once at the start of the main.

***

lattice_grid: prepare the grid for the path with periodic boundary conditions. It is called once at the start of the main.

***

momentum: produce the new fictuous moementa for the HMC algorithm. it also compute the initial kinetic energy. it is called once at the start of the main, for intialising the impulses and the kinetic energy and then, before every new cycle of measurement.

***

hmc: require as an input the time step of the integrator, the integration time and the constat cc that appears in front of the force (computed at the start). It implements the Omelyan integrator and returns the updated values of the field and of the momenta. it is called once for every measurement.

***

metropolis_update: it requires as an input the constat den, that appears in front of the potential, the current value of the charge and the integer step. step is a control variable:If step = 1, the code initialize the starting potential, kinetic energy and field value and update the value of step. the routine implements the metropolis update and returns the new or old values of the field, charge, potential, kinetic energy. It is called once for every measurement.

***

time_dependent_p: requires as an input the charge. it updates the value of the tdp, depending on the value of the charge, after 'tau' sweeps.

***

metad: it is a function that takes as an input the distances between the field values and returns the coefficient of the strength deriving from the tdp.

***

ran2: casual number generator

***

ranstart/ranfinish: subroutine for initializing ran2.


 
