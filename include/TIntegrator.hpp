// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TINTEGRATOR_HPP__
#define __TINTEGRATOR_HPP__

#include "TParticle.hpp"

#include "TPlot.hpp"

#include <vector>
#include "TMatrixtriang.hpp"

class TIntegrator
{
public:

    TIntegrator();

    /// Collection of particles
    std::vector<TParticle> particle_v;

    /// Temporal storage of individual forces
    TMatrixtriang forces_m;

    /// number of particles
    int nparticles = {0};

    /// Time step
    double h = {0.1};

    /// Step number
    int step = {0};

    /// Plot every nrefresh number of steps
    int nrefresh = {5};

    /// Calculate potentials and do the actual step
    void DoStep();

    /// Naive way of calculating the sum of foces for the ith particle
    void CalculateForce( std::vector<TVector> & f );

    /// Calculate force only for the i,j, i>j interaction and symmetrize using TMatrixtriang object
    void CalculateForceTriang( std::vector<TVector> & f );

    /// Plot positions of collection of particles
    void PlotPositions(TPlot & plt);

    /// Initialize n particles randomly
    void SetNparticlesRnd(int n);

    /// Calculate Total Mechanical energy. Total is returned, Keve and Vene are optional arguments
    double GetMene(double * Kene = nullptr, double * Vene = nullptr);

    /// Print total Mechanical energy
    void PrintMene();

    /// Raise exception if particles are too close
    void CheckDistances();

private:
    /// Implementation of the Euler FW integrator
    void IntegratorEulerFw( std::vector<TVector> & forces);

    /// Implementation of the Verlet integrator
    void IntegratorVerlet(std::vector< TVector >& force_v);
};

#endif
