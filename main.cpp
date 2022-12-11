// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include <iostream>

#include <thread>

#include "TVector.hpp"

#include "TParticle.hpp"

#include "TPlot.hpp"

#include <vector>

class TIntegrator
{
public:

    TIntegrator();
    std::vector<TParticle> particle_v;
    void DoStep();
    void PlotPositions(TPlot & plt);
    double h = {1.};
    int step = {0};
    void PrintMene();
    void CheckDistances();
    void IntegratorEulerFw( std::vector<TVector> & forces);

    void IntegratorVerlet( std::vector<TVector> & forces);
};

TIntegrator::TIntegrator():particle_v(2) {}

void TIntegrator::PrintMene()
{
    double Mene = 0;
    for( auto particle : particle_v)
    {
        const int id = particle.id;
        Mene += particle.Kene();
        for( auto aux : particle_v)
        {
            if( id == aux.id ) continue;
            Mene += aux.Vene( particle.pos );
        }

    }
    std::cout << "Step\t" << step << "\tMene\t" << Mene << std::endl;
}

void TIntegrator::PlotPositions(TPlot& plt)
{
    plt.StartPlot();

    for( auto & particle : particle_v )
        plt.AddPoint( particle.pos );

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


    return;
}

void TIntegrator::IntegratorEulerFw(std::vector<TVector>& force_v)
{
    auto force_i = force_v.begin();
    auto particle_i = particle_v.begin();

    for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
    {
        particle_i->pos += particle_i->vel*h;
        particle_i->pos += (*force_i)*(0.5*h*h);
        particle_i->vel += (*force_i);
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
    myIntegrator.PlotPositions(plt);


    for( int i = 0; i < 1000; ++i)
    {
        myIntegrator.DoStep();
        myIntegrator.CheckDistances();
        myIntegrator.PlotPositions(plt);
        myIntegrator.PrintMene();

    }

}
