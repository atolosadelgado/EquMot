This is a dummy project to test different simplectic integrands for the equations of motion. 

The class TVector implements basic routines for 2D vectors

The class TParticle stores position + velocity, and some other simple functions

The class TPlot encapsulates a pipe to GNU plot to show at the end of every step the position of particles

The class TIntegrator store the collection of TParticles, and implement different methods to integrate the equations of motion inside the function DoStep. At the moment only Euler forward and Verlet integrators are working. The stability was checked printing out the total mechanical energy. As expected, total energy in case of Euler fw increases over time, while Verlet integrator keeps it oscilating within the limits. The files EulerFW_mene.txt and Verlet_mene.txt contain the total mechanical energy for 5 orbits of the current implementation.
