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

    /// Calculate potentials and do the actual step, equivalent to ForceCXX+VerletCXX
    void DoStepCXX();

    /// Calculate potentials and do the actual step, equivalent to DoStepCXX but using TBB
    void DoStepTBB();

    /// Pointer to function which actually calculate the force
    std::function<void( std::vector<TVector> & f )> CalculateForce;

    /// Naive way of calculating the sum of foces for the ith particle
    void CalculateForceSimple( std::vector<TVector> & f );

    /// Implentation of CalculateForceSimple using OpenMP
    void CalculateForceSimplePar(std::vector<TVector> & f);

    /// Calculate force only for the i,j, i>j interaction and symmetrize using TMatrixtriang object
    void CalculateForceTriang( std::vector<TVector> & f );

    /// Implentation of CalculateForceTriang using OpenMP
    void CalculateForceTriangPar( std::vector<TVector> & f );

    /// Implentation of CalculateForceSimple using CXX17 TBB
    void CalculateForceSimpleParCXX( std::vector<TVector> & f );

    /// Set usage of CalculateForceSimple
    void SetForceCalSimple   () { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceSimple    (f); return;}; }

    /// Set usage of CalculateForceTriang
    void SetForceCalTriang   () { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceTriang    (f); return;}; }

     /// Set usage of CalculateForceTriangPar
    void SetForceCalTriangPar() { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceTriangPar (f); return;}; }

    /// Set usage of CalculateForceSimplePar
    void SetForceCalSimplePar() { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceSimplePar (f); return;}; }

    /// Set usage of CalculateForceSimpleParCXX
    void SetForceCalSimpleParCXX() { CalculateForce = [&](std::vector<TVector> & f){ this->CalculateForceSimpleParCXX(f); return;}; }


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


    /// Pointer to function which actually integrate the trajectories
    std::function<void(std::vector<TVector> & )> TheIntegrator;

    /// Implementation of the Euler FW integrator
    void IntegratorEulerFw( std::vector<TVector> & forces);

    /// Implementation of the Verlet integrator
    void IntegratorVerlet(std::vector< TVector >& force_v);

    /// Implementation of the Verlet integrator using CXX17 TBB
    void IntegratorVerletPar(std::vector< TVector >& force_v);

    /// Set usage of IntegratorEulerFw
    void SetIntegratorEulerFW   () { TheIntegrator = [&](std::vector<TVector> & f){ this->IntegratorEulerFw(f);   return;}; }

    /// Set usage of IntegratorVerlet
    void SetIntegratorVerlet    () { TheIntegrator = [&](std::vector<TVector> & f){ this->IntegratorVerlet(f);    return;}; }

    /// Set usage of IntegratorVerletPar
    void SetIntegratorVerletPar () { TheIntegrator = [&](std::vector<TVector> & f){ this->IntegratorVerletPar(f); return;}; }


    /// CPU demanding test for bencharming the force calculation method and the integrator method
    static void Test_benchmark_force_calculation(int i_forceMethod, int i_i_integratorMethod = 2);



};

#endif
