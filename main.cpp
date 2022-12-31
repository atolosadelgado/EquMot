// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //



#include "TPlot.hpp"
#include "TIntegrator.hpp"

#include "TMatrixtriang.hpp"

#include "Contants.hpp"

#include <fstream>
#include <iostream>
#include <thread>


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

void Test_Crystal()
{
    TPlot plt;

    TIntegrator myIntegrator;

    myIntegrator.h = 1;
    myIntegrator.nrefresh = 1000;

    TParticle::f_constant = G_UA_MSun;
    TParticle::fE_constant = K_SI_e_nm;
    TParticle::fL_constant = 1e-3;
    TParticle::frad_critical = 1.5;


    const int nx = 20;
    const int ny = 20;
    myIntegrator.SetNparticlesRnd(nx*ny);

    std::default_random_engine generator (0);

    /// Random distribution to place the particles
    const double eV = 1.6e-19; //J
    double v = sqrt( 0.025*eV*2/23/UAM_kg);
    v = 1e-5;
    std::uniform_real_distribution<double> rnd_vel (-v, v);

    for( int inx = 0; inx< nx; ++inx)
    {
        for( int iny = 0; iny < ny; ++iny)
        {
            myIntegrator.particle_v[inx*ny+iny].pos.x = inx;
            myIntegrator.particle_v[inx*ny+iny].pos.y = iny;
            myIntegrator.particle_v[inx*ny+iny].vel.x =    0 ; // rnd_vel(generator); //
            myIntegrator.particle_v[inx*ny+iny].vel.y = 0;
            myIntegrator.particle_v[inx*ny+iny].mass = 1; //1e15*UAM_kg;
            myIntegrator.particle_v[inx*ny+iny].charge = 1;
            myIntegrator.particle_v[inx*ny+iny].SetForce(4);;

            if( (inx+iny) % 2 )
            {
                myIntegrator.particle_v[inx*nx+iny].charge *= -1.;
//                 myIntegrator.particle_v[inx*nx+iny].isFixed = true;




            }

            // Keep 2 layers of atoms fixed
            if( 0 == inx*iny || 4 > inx || 4 > iny || nx-4 <= inx || ny-4 <= iny )
                myIntegrator.particle_v[inx*nx+iny].isFixed = true;

            if( nx-5 == inx && false ==  myIntegrator.particle_v[inx*nx+iny].isFixed )
                myIntegrator.particle_v[inx*ny+iny].vel.x = -1e-5; //




        }
    }

    myIntegrator.PlotPositions(plt);


    std::ofstream ofile("Verlet_2_crystal.txt");

    for( int istep = 0; istep < 5000000; ++istep)
    {
        myIntegrator.DoStepTBB();
        myIntegrator.PlotPositions(plt);

//         myIntegrator.CheckDistances();
//         myIntegrator.PrintMene();

        if( 0 == istep % myIntegrator.nrefresh )
            std::cout << "Step " << istep << std::endl;

        if( 0 == istep % myIntegrator.nrefresh )
        {
            double s2 = 0;
            double pos_bar = 0;
            for( int inx = 0; inx< nx; ++inx)
            {
                for( int iny = 0; iny < ny; ++iny)
                {
                    TParticle & p = myIntegrator.particle_v.at(inx*ny+iny);

                    double dx = (p.pos.x - double(inx));
                    double dy = (p.pos.y - double(iny));
                    double r2 = dx*dx + dy*dy;
                    pos_bar+= sqrt( r2 );
                    s2 += r2;
                }
            }
            ofile << myIntegrator.GetMene() << '\t' << pos_bar/TParticle::nparticles << '\t' << sqrt(s2/TParticle::nparticles) << std::endl;
        }

//         if( 100 == istep )
//         {
//                 for( int i = 0; i< ny; ++i)
//             {
//                 myIntegrator.particle_v[(nx-2)*ny+i].vel.x -=0.5e-5;
//             }
//             std::cin.ignore();
//         }

    }

    ofile.close();


}

void Test_Spring()
{
    TPlot plt;

    TIntegrator myIntegrator;

    myIntegrator.h = 0.01;
    myIntegrator.nrefresh = 1;

    TParticle::fL_constant = 1.0;
//     TParticle::frad_critical = 15;


    myIntegrator.SetNparticlesRnd(4);



    myIntegrator.particle_v[0].pos.x = 0;
    myIntegrator.particle_v[0].pos.y = 0;
    myIntegrator.particle_v[0].vel.x = 0;
    myIntegrator.particle_v[0].vel.y = 0;
    myIntegrator.particle_v[0].isFixed = true;

    myIntegrator.particle_v[0].SetForce(4);

    myIntegrator.particle_v[1].pos.x = 20;
    myIntegrator.particle_v[1].pos.y = 0;
    myIntegrator.particle_v[1].vel.x = 0;
    myIntegrator.particle_v[1].vel.y = 0;
    myIntegrator.particle_v[1].isFixed = true;
    myIntegrator.particle_v[1].SetForce(4);


    myIntegrator.particle_v[2].pos.x = 7;
    myIntegrator.particle_v[2].pos.y = 0;
    myIntegrator.particle_v[2].vel.x = -1;
    myIntegrator.particle_v[2].vel.y = 0;
    myIntegrator.particle_v[2].mass  = 1;
    myIntegrator.particle_v[2].isFixed = false;
    myIntegrator.particle_v[2].SetForce(4);

    myIntegrator.particle_v[3].pos.x = 17;
    myIntegrator.particle_v[3].pos.y = 0;
    myIntegrator.particle_v[3].vel.x = 1;
    myIntegrator.particle_v[3].vel.y = 0;
    myIntegrator.particle_v[3].mass  = 1;
    myIntegrator.particle_v[3].isFixed = false;
    myIntegrator.particle_v[3].SetForce(4);

    myIntegrator.PlotPositions(plt);

    for( int istep = 0; istep < 5000000; ++istep)
    {
        myIntegrator.DoStepTBB();
        myIntegrator.PlotPositions(plt);

        if( 0 == istep % myIntegrator.nrefresh )
            std::cout << "Step " << istep << std::endl;

    }

}

void Test_Spring2()
{



    TPlot plt;

    TIntegrator myIntegrator;

    myIntegrator.h = 0.01;
    myIntegrator.nrefresh = 1;

    TParticle::fL_constant = 1.0;
//     TParticle::frad_critical = 15;


    myIntegrator.SetNparticlesRnd(4);



    myIntegrator.particle_v[0].pos.x = 0;
    myIntegrator.particle_v[0].pos.y = 0;
    myIntegrator.particle_v[0].vel.x = 0;
    myIntegrator.particle_v[0].vel.y = 0;
    myIntegrator.particle_v[0].isFixed = true;

    myIntegrator.particle_v[0].SetForce(4);

    myIntegrator.particle_v[1].pos.x = 18;
    myIntegrator.particle_v[1].pos.y = 0;
    myIntegrator.particle_v[1].vel.x = 0;
    myIntegrator.particle_v[1].vel.y = 0;
    myIntegrator.particle_v[1].isFixed = true;
    myIntegrator.particle_v[1].SetForce(4);


    myIntegrator.particle_v[2].pos.x = 6;
    myIntegrator.particle_v[2].pos.y = 0;
    myIntegrator.particle_v[2].vel.x = 0;
    myIntegrator.particle_v[2].vel.y = 1;
    myIntegrator.particle_v[2].mass  = 1;
    myIntegrator.particle_v[2].isFixed = false;
    myIntegrator.particle_v[2].SetForce(4);

    myIntegrator.particle_v[3].pos.x = 12;
    myIntegrator.particle_v[3].pos.y = 0;
    myIntegrator.particle_v[3].vel.x = 0;
    myIntegrator.particle_v[3].vel.y = -1;
    myIntegrator.particle_v[3].mass  = 1;
    myIntegrator.particle_v[3].isFixed = false;
    myIntegrator.particle_v[3].SetForce(4);

    myIntegrator.PlotPositions(plt);
    myIntegrator.SetCriticalRadius(7);

    for( int istep = 0; istep < 5000000; ++istep)
    {
        myIntegrator.DoStepTBB();
        myIntegrator.PlotPositions(plt);

        if( 0 == istep % myIntegrator.nrefresh )
            std::cout << "Step " << istep << std::endl;

    }



    return;
};

int main()
{

    TPlot plt("plt_positions");
    TPlot plt_kene("plt_kene",true);

    TIntegrator myIntegrator;

    myIntegrator.h = 20;
    myIntegrator.nrefresh = 1000;
    myIntegrator.sleep_time_ms = 0;

    TParticle::f_constant = G_UA_MSun;
    TParticle::fE_constant = K_SI_e_nm;
    TParticle::fL_constant = 1e-2;


    const int nx = 200;
    const int ny = 50;
    myIntegrator.SetNparticlesRnd(nx*ny);

    std::default_random_engine generator (0);

    /// Random distribution to place the particles
    const double eV = 1.6e-19; //J
    double v = sqrt( 0.025*eV*2/23/UAM_kg);
    v = 1e-5;
    std::uniform_real_distribution<double> rnd_vel (-v, v);

    for( int inx = 0; inx< nx; ++inx)
    {
        for( int iny = 0; iny < ny; ++iny)
        {
            TParticle & p = myIntegrator.particle_v.at(inx*ny+iny);
            std::cout << inx << " " << iny << std::endl;
            p.pos.x = inx;
            p.pos.y = iny;
            p.vel.x = 0; // rnd_vel(generator);
            p.vel.y = 0; // rnd_vel(generator);
            p.mass = 1e5;
            p.SetForce(4);;

            // Keep particles in the corner fixed
            if(  0 == inx ||  0 == iny || nx-1 == inx || ny-1 == iny)
            {
                p.isFixed = true;
//                 p.mass *=10000;
                p.vel.x = 0;
                p.vel.y = 0;
                if( 0 == inx&& 0 == iny )
                    p.vel.x = 0.1*v;

            }

        }
    }
//     return 0;

    myIntegrator.PlotPositions(plt);
    myIntegrator.SetCriticalRadius(1.1);
//     myIntegrator.SetDamping( -100.0 );
//     return 0;

    std::ofstream ofile("Verlet_2_crystal.txt");
    std::ofstream ofilePos("Verlet_2_crystal_pos.txt");


    for( int istep = 0; istep < 5000000; ++istep)
    {
        myIntegrator.DoStepTBB();
//          std::cin.ignore();
        myIntegrator.PlotPositions(plt);
        myIntegrator.PlotPositionsW(plt_kene, 'K');

//         myIntegrator.CheckDistances();
//         myIntegrator.PrintMene();

        if( 3000 == istep )
        {
            std::cin.ignore();
            myIntegrator.particle_v[2*ny+25].vel.x = -v;
        }

//         for( int i = 3; i < 8; ++i)
//         {
//         myIntegrator.particle_v[i*ny+i].vel.x = 0;
//         myIntegrator.particle_v[i*ny+i].vel.y = 0;
//         }

//         if( 4e4 == istep)
//             myIntegrator.SetDampingZero();



        double kene, vene;
        if( 0 == istep % myIntegrator.nrefresh )
        {
            double s2 = 0;
            double pos_bar = 0;
            for( int inx = 0; inx< nx; ++inx)
            {
                for( int iny = 0; iny < ny; ++iny)
                {
                    TParticle & p = myIntegrator.particle_v.at(inx*ny+iny);

                    double dx = (p.pos.x - double(inx));
                    double dy = (p.pos.y - double(iny));
                    double r2 = dx*dx + dy*dy;
                    pos_bar+= sqrt( r2 );
                    s2 += r2;
                }
            }
            ofile << myIntegrator.GetMene( &kene, &vene ) << '\t' << pos_bar/TParticle::nparticles << '\t' << sqrt(s2/TParticle::nparticles) << std::endl;
            ofilePos << myIntegrator.particle_v[10*nx+10].pos << std::endl;
        }

        if( 0 == istep % myIntegrator.nrefresh )
            std::cout << "Step " << istep << "\tKene " << kene << std::endl;



    }

    ofile.close();
    ofilePos.close();



}
