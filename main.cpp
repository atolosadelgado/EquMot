// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include <iostream>

#include <thread>

#include "TVector.hpp"

#include "TParticle.hpp"

#include "TPlot.hpp"

#include <vector>
#include <iomanip>
#include <fstream>

class TIntegrator
{
public:

    TIntegrator();
    std::vector<TParticle> particle_v;
    void DoStep();
    void PlotPositions(TPlot & plt);
    double h = {0.1};
    int step = {0};
    int nrefresh = {5};
    double GetMene(double * Kene, double * Vene);
    void PrintMene();
    void CheckDistances();
    void IntegratorEulerFw( std::vector<TVector> & forces);

    void IntegratorVerlet(std::vector< TVector >& force_v);
};

TIntegrator::TIntegrator():particle_v(2) {}

double TIntegrator::GetMene(double* Kene = nullptr, double* Vene = nullptr)
{
    double _Kene = 0;
    double _Vene = 0;
    for( auto particle : particle_v)
    {
        const int id = particle.id;
        _Kene += particle.Kene();
        for( auto aux : particle_v)
        {
            if( id == aux.id ) continue;
            _Vene += aux.Vene( particle.pos );
        }

    }
    if( Kene ) *Kene = _Kene;
    if( Vene ) *Vene = _Vene;
    return _Vene+_Kene;
}


void TIntegrator::PrintMene()
{
    double Kene = 0;
    double Vene = 0;
    GetMene( &Kene,  &Vene);
    std::cout << std::scientific << std::setprecision(3)
              << "Step\t" << step
              << "\tMene\t" << Kene+Vene
              << "\tKene\t" << Kene
              << "\tVene\t" << Vene
              << std::endl;
}

void TIntegrator::PlotPositions(TPlot& plt)
{
  if( 0 != step % nrefresh ) return;
    plt.StartPlot();

    for( auto & particle : particle_v )
    {
        plt.AddPoint( particle.pos );
        std::cout << particle << std::endl;
    }

    plt.ShowPlot();

    std::this_thread::sleep_for(std::chrono::milliseconds(10));
}

void TIntegrator::DoStep()
{

    ++step;
    std::vector<TVector> force_v( particle_v.size() );
    // For each particle, sum up the force of the other particles
    {
        auto force_i = force_v.begin();
        auto particle_i = particle_v.begin();

        for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
        {
            const int id = particle_i->id;
            for( auto aux : particle_v)
            {
                if( id == aux.id ) continue;
                *force_i += aux.Force( particle_i->pos );
            }

        }
    }
    // apply the force and do the actual step
    IntegratorEulerFw(force_v);
//     IntegratorVerlet(force_v);


    return;
}

void TIntegrator::IntegratorEulerFw(std::vector<TVector>& force_v)
{
    auto force_i = force_v.begin();
    auto particle_i = particle_v.begin();

    for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
    {
        particle_i->pos += particle_i->vel.Scale(h);;
        particle_i->pos += (*force_i).Scale(0.5*h*h);
        particle_i->vel += (*force_i).Scale(h);
    }
}

void TIntegrator::IntegratorVerlet(std::vector<TVector>& force_v)
{
    auto force_i = force_v.begin();
    auto particle_i = particle_v.begin();
    for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
    {
        particle_i->pos += particle_i->vel.Scale(h);
        particle_i->pos += (*force_i).Scale(0.5*h*h);
    }

    //recalculate forces at new positions, and sum force at the next position to the previous calculated force
    {
        auto force_i = force_v.begin();
        auto particle_i = particle_v.begin();

        for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
        {
            const int id = particle_i->id;
            for( auto aux : particle_v)
            {
                if( id == aux.id ) continue;
                *force_i += aux.Force( particle_i->pos );
            }
            // update Velocity
            particle_i->vel += (( *force_i ).Scale(0.5*h));

        }
    }

}



void TIntegrator::CheckDistances()
{

    for( auto particle : particle_v)
    {
        const int id = particle.id;
        for( auto aux : particle_v)
        {
            if( id == aux.id ) continue;
            if( 1 > aux.pos.distance( particle.pos ) )
                std::runtime_error("Too close particles!!\n");
        }

    }

}


int main() {


    TPlot plt;

    TIntegrator myIntegrator;

    myIntegrator.particle_v[0].pos.x = -30;
    myIntegrator.particle_v[0].pos.y = 0;
    myIntegrator.particle_v[0].vel.x = 0;
    myIntegrator.particle_v[0].vel.y = -3;

    myIntegrator.particle_v[1].pos.x = 30;
    myIntegrator.particle_v[1].pos.y = 0;
    myIntegrator.particle_v[1].vel.x = 0;
    myIntegrator.particle_v[1].vel.y = 3;


    myIntegrator.PlotPositions(plt);

//     std::ofstream ofile("EulerFW.txt");


    for( int i = 0; i < 5000; ++i)
    {
        myIntegrator.DoStep();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
        myIntegrator.PrintMene();
//         ofile << myIntegrator.GetMene() << std::endl;

    }

//     ofile.close();

}
