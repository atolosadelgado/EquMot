// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



#include "TPlot.hpp"
#include "TIntegrator.hpp"

#include "TMatrixtriang.hpp"

#include "Contants.hpp"

#include <fstream>
#include <iostream>


void Test_SS_acretion()
{
    TPlot plt;

    TIntegrator myIntegrator;
    myIntegrator.h = 3600;
    myIntegrator.nrefresh = 1000;

    TParticle::f_constant = G_UA_MSun;

    myIntegrator.SetNparticlesRnd(100);

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
    for( int i=3; i<100; ++i)
    {
        myIntegrator.particle_v[i].mass =  moon_mass / mass_sun /100 ;
        myIntegrator.particle_v[i].vel.x = rnd_vel(generator) ;
        myIntegrator.particle_v[i].vel.y = rnd_vel(generator) ;

    }



    myIntegrator.PlotPositions(plt);



    std::ofstream ofile("Total_mechanical_energy_SolarSystem_accretion.txt");

    for( int i = 0; i < 50000; ++i)
    {
        myIntegrator.DoStep();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
        myIntegrator.PrintMene();

//         std::cin.ignore();
        if( 0 == i % 10 )  ofile << myIntegrator.GetMene() << std::endl;

    }

    ofile.close();

}


void Test_SS()
{

    TPlot plt;

    TIntegrator myIntegrator;
    myIntegrator.h = 3600;
    myIntegrator.nrefresh = 100;

    TParticle::f_constant = G_UA_MSun;

    myIntegrator.SetNparticlesRnd(4);


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



    // Jupiter
    myIntegrator.particle_v[2].pos.x = -5.2;
    myIntegrator.particle_v[2].pos.y = 0;
    myIntegrator.particle_v[2].vel.x = 0;
    myIntegrator.particle_v[2].vel.y = -jupiter_orbital_speed_UA;
    myIntegrator.particle_v[2].mass =  mass_jupiter / mass_sun;
    std::cout << myIntegrator.particle_v[2] << std::endl;


    // Saturn
    myIntegrator.particle_v[3].pos.x = -saturn_distance;
    myIntegrator.particle_v[3].pos.y = 0;
    myIntegrator.particle_v[3].vel.x = 0;
    myIntegrator.particle_v[3].vel.y = -saturn_orbital_speed_UA;
    myIntegrator.particle_v[3].mass =  saturn_mass / mass_sun;
    std::cout << myIntegrator.particle_v[3] << std::endl;

    myIntegrator.PlotPositions(plt);



    std::ofstream ofile("Total_mechanical_energy_SolarSystem.txt");

    for( int i = 0; i < 500000; ++i)
    {
//         myIntegrator.DoStep();
        myIntegrator.DoStepCXX();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
//         myIntegrator.PrintMene();

//         std::cin.ignore();
        if( 0 == i % 10 )  ofile << myIntegrator.GetMene() << std::endl;

    }

    ofile.close();


}

int main()
{
    Test_SS();
    return 0;
}
