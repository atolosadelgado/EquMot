// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TIntegrator.hpp"

#include <iomanip>
#include <iostream>
#include <thread>

TIntegrator::TIntegrator(){}

double TIntegrator::GetMene(double* Kene, double* Vene)
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
            _Vene += aux.Vene( particle );
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
#ifdef NDEBUG
        std::cout << particle << std::endl;
#endif
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
                *force_i += aux.Force( *particle_i );
            }

        }
    }
    // apply the force and do the actual step
//     IntegratorEulerFw(force_v);
    IntegratorVerlet(force_v);


    return;
}

void TIntegrator::IntegratorEulerFw(std::vector<TVector>& force_v)
{
    auto force_i = force_v.begin();
    auto particle_i = particle_v.begin();

    for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
    {
        if( particle_i->isFixed ) continue;
//         particle_i->pos += particle_i->vel.Scale(h);;
//         particle_i->pos += (*force_i).Scale(0.5*h*h).Scale( 1./particle_i->mass );
//         particle_i->vel += (*force_i).Scale(h).Scale( 1./particle_i->mass );
        particle_i->pos.Add( particle_i->vel , h);;
        particle_i->pos.Add( (*force_i), 0.5*h*h/particle_i->mass );
        particle_i->vel.Add( (*force_i), h/particle_i->mass );

    }
}

void TIntegrator::IntegratorVerlet(std::vector<TVector>& force_v)
{
    auto force_i = force_v.begin();
    auto particle_i = particle_v.begin();
    for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
    {
        if( particle_i->isFixed ) continue;
//         particle_i->pos += particle_i->vel.Scale(h);
//         particle_i->pos += (*force_i).Scale(0.5*h*h).Scale( 1./particle_i->mass );
        particle_i->pos.Add( particle_i->vel, h);
        particle_i->pos.Add( (*force_i), 0.5*h*h/particle_i->mass );
    }

    //recalculate forces at new positions, and sum force at the next position to the previous calculated force
    {
        auto force_i = force_v.begin();
        auto particle_i = particle_v.begin();

        for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
        {
            if( particle_i->isFixed ) continue;
            const int id = particle_i->id;
            for( auto aux : particle_v)
            {
                if( id == aux.id ) continue;
                *force_i += aux.Force( *particle_i );
            }
            // update Velocity
            particle_i->vel.Add( ( *force_i ), 0.5*h/particle_i->mass );

        }
    }

}

void TIntegrator::SetNparticlesRnd(int n)
{
    particle_v = std::vector<TParticle>(n);
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
