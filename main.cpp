// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



#include "TPlot.hpp"
#include "TIntegrator.hpp"

#include "TMatrixtriang.hpp"

#include "Contants.hpp"

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

void Test_ensemble2()
{
    TPlot plt;

    TIntegrator myIntegrator;
    myIntegrator.h = 1e-3;
    myIntegrator.nrefresh = 1;

    TParticle::f_constant = 10;

    myIntegrator.SetNparticlesRnd(2);

    myIntegrator.particle_v[0].pos.x = -1;
    myIntegrator.particle_v[0].pos.y = 0;
    myIntegrator.particle_v[0].vel.x = 0;
    myIntegrator.particle_v[0].vel.y = -1;
    myIntegrator.particle_v[0].mass = 1;
    std::cout << myIntegrator.particle_v[0] << std::endl;

    myIntegrator.particle_v[1].pos.x = 1;
    myIntegrator.particle_v[1].pos.y = 0;
    myIntegrator.particle_v[1].vel.x = 0;
    myIntegrator.particle_v[1].vel.y = 1;

    myIntegrator.particle_v[1].mass =  1;


    std::cout << myIntegrator.particle_v[1] << std::endl;

    myIntegrator.PlotPositions(plt);



    for( int i = 0; i < 50000; ++i)
    {
        myIntegrator.DoStep();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
        myIntegrator.PrintMene();

    }

}

void Test_TMatrix(int n = 10)
{
    TMatrixtriang m(n);
    TVector dummyval( 3.1415, 0);
    m.SetVal(1,2, dummyval);
    std::cout << m << std::endl;
    if( dummyval.x != m.GetVal(2,1).x )
    {
        throw std::runtime_error("TMatrixtriang do not work\n");
    }
    return;
}

int main()
{

//     Test_SS();
//     return 0;



    TPlot plt;

    TIntegrator myIntegrator;

    myIntegrator.h = 3600;
    myIntegrator.nrefresh = 1000;

    TParticle::f_constant = G_UA_MSun;

    myIntegrator.SetNparticlesRnd(500);

    // Sun
    myIntegrator.particle_v[0].pos.x = 0;
    myIntegrator.particle_v[0].pos.y = 0;
    myIntegrator.particle_v[0].vel.x = 0;
    myIntegrator.particle_v[0].vel.y = 0;
    myIntegrator.particle_v[0].mass = 1.0;
    myIntegrator.particle_v[0].isFixed = true;


    std::cout << myIntegrator.particle_v[0] << std::endl;



    // Earth
    myIntegrator.particle_v[1].pos.x = 1.;
    myIntegrator.particle_v[1].pos.y = 0;
    myIntegrator.particle_v[1].vel.x = 0;
    myIntegrator.particle_v[1].vel.y = earth_obital_speed_UA;
    myIntegrator.particle_v[1].mass =  mass_earth / mass_sun;
    std::cout << myIntegrator.particle_v[1] << std::endl;



    // Moon
    myIntegrator.particle_v[2].pos.x = 1 - moon_distance_to_earth/UA_m;
    myIntegrator.particle_v[2].pos.y = 0;
    myIntegrator.particle_v[2].vel.x = 0;
    myIntegrator.particle_v[2].vel.y = earth_obital_speed_UA - moon_obital_speed_UA;
    myIntegrator.particle_v[2].mass =  moon_mass / mass_sun;
    std::cout << myIntegrator.particle_v[2] << std::endl;

    std::default_random_engine generator (0);

/// Random distribution to place the particles
    std::uniform_real_distribution<double> rnd_vel (-earth_obital_speed_UA,earth_obital_speed_UA);
    for( int i=3; i<500; ++i)
    {
        myIntegrator.particle_v[i].mass =  moon_mass / mass_sun /100 ;
        myIntegrator.particle_v[i].vel.x = rnd_vel(generator) ;
        myIntegrator.particle_v[i].vel.y = rnd_vel(generator) ;

    }



//     myIntegrator.PlotPositions(plt);



//     std::ofstream ofile("Verlet_2.txt");

    for( int i = 0; i < 500; ++i)
    {
        myIntegrator.DoStep();
//         myIntegrator.CheckDistances();
//         myIntegrator.PlotPositions(plt);
//         myIntegrator.PrintMene();

//         std::cin.ignore();
//         if( 0 == i % 10 )  ofile << myIntegrator.GetMene() << std::endl;
        if( 0 == i % 100 ) std::cout << "Step " << i << std::endl;

    }

//     ofile.close();


}
