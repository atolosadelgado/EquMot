// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



#include "TPlot.hpp"
#include "TIntegrator.hpp"

#include <fstream>
#include <iostream>

void Test_ensemble()
{
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

    std::ofstream ofile("Verlet_2.txt");


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

int main()
{
    Test_ensemble();
}
