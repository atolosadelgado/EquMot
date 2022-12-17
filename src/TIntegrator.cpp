// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TIntegrator.hpp"

#include <iomanip>
#include <iostream>
#include <thread>

TIntegrator::TIntegrator() {}

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

void TIntegrator::CalculateForce(std::vector<TVector>& force_v)
{
    // For each particle, sum up the force of the other particles
//     auto force_i = force_v.begin();
//     auto particle_i = particle_v.begin();

    #pragma omp parallel for
    for( int i = 0 ; i < nparticles ; ++i)
    {
        const int id = particle_v[i].id;
        for( auto aux : particle_v)
        {
            if( id == aux.id ) continue;
            force_v[i] += aux.Force( particle_v[i] );
        }

    }
//     for( ; particle_i != particle_v.end(); ++force_i, ++particle_i)
//     {
//         const int id = particle_i->id;
//         for( auto aux : particle_v)
//         {
//             if( id == aux.id ) continue;
//             *force_i += aux.Force( *particle_i );
//         }
//
//     }
    return;
}

void TIntegrator::CalculateForceTriang(std::vector<TVector>& f)
{
    if( 1 >= nparticles ) return;
    /// First reset the matrix
    TVector a(0.,0.);
    std::fill( forces_m.matrix.begin(), forces_m.matrix.end(), a);
//  Paralelizable aproach
    #pragma omp parallel for
    for( int particle_row = 0; particle_row < nparticles; ++particle_row)
    {
        for( int particle_col = particle_row+1; particle_col < nparticles; ++particle_col)
        {
//             forces_m.SetVal( particle_row +1, particle_col +1, particle_v[particle_row].Force( particle_v[particle_col] ));
            forces_m.at( particle_row +1, particle_col +1) = ( particle_v[particle_row].Force( particle_v[particle_col] ));
        }
    }
    #pragma omp parallel for
    for( int particle_row = 0; particle_row < nparticles; ++particle_row)
    {
        for( int particle_col = 0; particle_col < nparticles; ++particle_col)
        {
            if( particle_col == particle_row ) continue;
            f[particle_row] += forces_m.GetVal( particle_row+1, particle_col +1 );
        }
    }

//     for( int particle_row = 0; particle_row < nparticles; ++particle_row)
//     {
//         for( int particle_col = 0; particle_col < particle_row; ++particle_col)
//         {
//             f[particle_row] += forces_m.GetVal( particle_row+1 , particle_col +1 );
//         }
//         for( int particle_col = particle_row+1; particle_col < nparticles; ++particle_col)
//         {
//             auto v = particle_v[particle_row].Force( particle_v[particle_col] );
//             forces_m.SetVal( particle_row +1, particle_col +1, v);
//             f[particle_row] += v;
//         }
//     }


    return;
}


void TIntegrator::DoStep()
{
    ++step;
    std::vector<TVector> force_v( particle_v.size() );
    CalculateForce( force_v );
//     CalculateForceTriang( force_v );

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
        particle_i->pos.Add( particle_i->vel, h);;
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
    forces_m.SetSize(n);
    nparticles = n;
}


void TIntegrator::CheckDistances()
{

    for( auto particle : particle_v)
    {
        const int id = particle.id;
        for( auto aux : particle_v)
        {
            if( id == aux.id ) continue;
            if( 1e-9 > aux.pos.distance( particle.pos ) )
                throw std::runtime_error("Too close particles!!\n");
        }

    }

}
