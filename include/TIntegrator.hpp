// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TINTEGRATOR_HPP__
#define __TINTEGRATOR_HPP__

#include "TParticle.hpp"

#include "TPlot.hpp"

#include <vector>

class TIntegrator
{
public:

    TIntegrator();

    std::vector<TParticle> particle_v;

    double h = {0.1};
    int step = {0};
    int nrefresh = {5};

    void DoStep();
    void PlotPositions(TPlot & plt);

    void SetNparticlesRnd(int n);

    double GetMene(double * Kene = nullptr, double * Vene = nullptr);
    void PrintMene();

    void CheckDistances();

private:


    void IntegratorEulerFw( std::vector<TVector> & forces);

    void IntegratorVerlet(std::vector< TVector >& force_v);
};

#endif
