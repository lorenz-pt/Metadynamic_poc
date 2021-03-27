# Metadynamics_poc
Implementation of Metadynamics to a simple topological system (A particle confined on a circle).

The program uses an HMC (Hamiltonian Monte Carlo or Hybrid Monte carlo) algorithm to simulate the paths distribution with a probability density given by the action functional. It makes use of the second order Omelyan symplectic integrator.

Metadynamics is needed for healing topological freezing at low temperature (in the algorithm, this is equivalent to small parameter 'e' or great values of 'N', the length of the paths).

The code, with the preset values of the parameters, takes about 1:40 minutes.
