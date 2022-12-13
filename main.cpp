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

const double mass_sun   = 1.989e30; //kg
const double UA_m       =1.495978707e11; //m

const double G_SI       = 6.6743e-11; // m3 kg-1 s-2
const double G_UA_MSun  = G_SI / pow(1.495978707, 3) * 1.989e-3; // UA3 msun-1 s-2

const double mass_earth = 5.972e24; //kg
const double earth_obital_speed = 29.78e3; // m/s
const double earth_obital_speed_UA = 29.78e3/UA_m; // UA/s


const double moon_distance_to_earth = 384400e3; //m
const double moon_mass =  7.342e22; //kg
const double moon_obital_speed = 1.0e3; // m/s
const double moon_obital_speed_UA = moon_obital_speed/UA_m; // UA/s



const double saturn_distance = 9.5; //UA
const double saturn_mass= 5.683e26; //kg
const double saturn_orbital_speed =9.69e3; // m/s
const double saturn_orbital_speed_UA = saturn_orbital_speed/UA_m; // UA/s


const double distance_jupiter = 5.2; //UA
const double mass_jupiter = 1.898e27; //kg
const double jupiter_orbital_speed = 13.07e3; // m/s
const double jupiter_orbital_speed_UA = jupiter_orbital_speed/UA_m; // UA/s


const double seconds_per_year = 365.26*24*3600;
const double G_UA_MSun_yr = G_SI / pow(1.495978707, 3) * 1.989e-3 * pow(seconds_per_year,2); // UA3 msun-1 s-2
const double earth_obital_speed_UA_yr = earth_obital_speed/UA_m*seconds_per_year; // m/s

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



    std::ofstream ofile("Verlet_2.txt");

    for( int i = 0; i < 500000; ++i)
    {
        myIntegrator.DoStep();
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

    TPlot plt;

    TIntegrator myIntegrator;
    myIntegrator.h = 3600;
    myIntegrator.nrefresh = 100;

    TParticle::f_constant = G_UA_MSun;

    myIntegrator.SetNparticlesRnd(3);

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



    myIntegrator.PlotPositions(plt);



    std::ofstream ofile("Verlet_2.txt");

    for( int i = 0; i < 50000; ++i)
    {
        myIntegrator.DoStep();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
//         myIntegrator.PrintMene();

//         std::cin.ignore();
        if( 0 == i % 10 )  ofile << myIntegrator.GetMene() << std::endl;

    }

    ofile.close();


}
