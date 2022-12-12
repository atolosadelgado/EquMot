// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include <iostream>


#include "TVector.hpp"

#include "TParticle.hpp"

#include "TPlot.hpp"

#include "TIntegrator.hpp"

#include <vector>
#include <fstream>


int main() {


    TPlot plt;

    TIntegrator myIntegrator;

    myIntegrator.SetNparticlesRnd(2);

    myIntegrator.particle_v[0].pos.x = -30;
    myIntegrator.particle_v[0].pos.y = 0;
    myIntegrator.particle_v[0].vel.x = 0;
    myIntegrator.particle_v[0].vel.y = -3;

    myIntegrator.particle_v[1].pos.x = 30;
    myIntegrator.particle_v[1].pos.y = 0;
    myIntegrator.particle_v[1].vel.x = 0;
    myIntegrator.particle_v[1].vel.y = 3;


    myIntegrator.PlotPositions(plt);

    std::ofstream ofile("EulerFW_2.txt");


    for( int i = 0; i < 5000; ++i)
    {
        myIntegrator.DoStep();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
        myIntegrator.PrintMene();
        ofile << myIntegrator.GetMene() << std::endl;

    }

    ofile.close();

}
