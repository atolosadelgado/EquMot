// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TINTEGRATOR_HPP__
#define __TINTEGRATOR_HPP__

#include "TParticle.hpp"

#include "TPlot.hpp"

#include <vector>
#include "TMatrixtriang.hpp"

#include <functional>

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

    std::function<void( std::vector<TVector> & f )> CalculateForce;

    void SetForceCalSimple   () { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceSimple    (f); return;}; }
    void SetForceCalTriang   () { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceTriang    (f); return;}; }
    void SetForceCalTriangPar() { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceTriangPar (f); return;}; }
    void SetForceCalSimplePar() { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceSimplePar (f); return;}; }


    /// Naive way of calculating the sum of foces for the ith particle
    void CalculateForceSimple( std::vector<TVector> & f );
    void CalculateForceSimplePar(std::vector<TVector> & f);


    /// Calculate force only for the i,j, i>j interaction and symmetrize using TMatrixtriang object
    void CalculateForceTriang( std::vector<TVector> & f );
    void CalculateForceTriangPar( std::vector<TVector> & f );

    /// simple way + TBB, linkage disable in cmake termporary
    void CalculateForceParCXX( std::vector<TVector> & f );

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

    static void Test_benchmark_force_calculation(int i);

private:
    /// Implementation of the Euler FW integrator
    void IntegratorEulerFw( std::vector<TVector> & forces);

    /// Implementation of the Verlet integrator
    void IntegratorVerlet(std::vector< TVector >& force_v);
};

#endif
