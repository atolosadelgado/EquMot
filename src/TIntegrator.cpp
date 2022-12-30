// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TIntegrator.hpp"

#include <iomanip>
#include <iostream>
#include <thread>
#include <execution>
#include <algorithm>

#include "Contants.hpp"

#include <tbb/parallel_for.h>

TIntegrator::TIntegrator() {
    SetForceCalSimple();
    SetIntegratorVerlet();
}

void TIntegrator::SetDampingZero()
{
    fDamp=0.0;
    return;
}

void TIntegrator::SetDamping(double d)
{
    if( d<= 0) fDamp = d;
    else       throw std::runtime_error("Damping factor can not be positive!");
    return;
}


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
#ifndef NDEBUG
        std::cout << particle << std::endl;
#endif
    }

    plt.ShowPlot();

    std::this_thread::sleep_for(std::chrono::milliseconds(1));
}

void TIntegrator::PlotPositionsW(TPlot& plt, char )
{
    if( 0 != step % nrefresh ) return;
    auto W_func = [](TParticle & p){ return p.Kene();};
// //     TODO: Rewrite function TParticle::Vene and TParticle::Mene in order to be usable here!
//     switch( W )
//     {
//         case 'K': break;
//         case 'V': break;
//         case 'M': break;
//
//     };



            plt.StartH2D();

            for( auto & particle : particle_v )
            {
                plt.AddPointH2D( particle.pos , W_func(particle) );
#ifndef NDEBUG
                std::cout << particle << std::endl;
#endif
            }

            plt.ShowPlot();

            std::this_thread::sleep_for(std::chrono::milliseconds(1));

}


void TIntegrator::CalculateForceSimpleParCXX(std::vector<TVector> & f)
{
    for_each( std::execution::par, particle_v.begin(), particle_v.end(), [](TParticle & p) {
//     for_each( particle_v.begin(), particle_v.end(), [](TParticle & p) {
        p.force.x =0;
        p.force.y=0;
        return;
    });

    auto ff = [&](TParticle & p) {
        const int id = p.id;
        for( auto aux : particle_v)
        {
            if( id == aux.id ) continue;
            p.force += aux.Force( p );
        }


    };

    for_each(std::execution::par, particle_v.begin(), particle_v.end(), ff);
//     for_each(particle_v.begin(), particle_v.end(), ff);

    std::transform(std::execution::par,  particle_v.begin(), particle_v.end(), f.begin(), [](TParticle & p) {
        return p.force;
    } );


}


void TIntegrator::CalculateForceSimplePar(std::vector<TVector>& force_v)
{
    // For each particle, sum up the force of the other particles

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
    return;
}

void TIntegrator::CalculateForceSimple(std::vector<TVector>& force_v)
{
    // For each particle, sum up the force of the other particles
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
    return;
}

void TIntegrator::CalculateForceTriangPar(std::vector<TVector>& f)
{
    if( 1 >= nparticles ) return;
    /// First reset the matrix
    TVector a(0.,0.);
    std::fill( forces_m.matrix.begin(), forces_m.matrix.end(), a);

    #pragma omp parallel for
    for( int particle_row = 0; particle_row < nparticles; ++particle_row)
    {
        for( int particle_col = particle_row+1; particle_col < nparticles; ++particle_col)
        {
            forces_m.at( particle_row +1, particle_col +1) = ( particle_v[particle_col].Force( particle_v[particle_row] ));
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
    return;
}

void TIntegrator::CalculateForceTriang(std::vector<TVector>& f)
{

// Sequential
    for( int particle_row = 0; particle_row < nparticles; ++particle_row)
    {
        for( int particle_col = 0; particle_col < particle_row; ++particle_col)
        {
            f[particle_row] += forces_m.GetVal( particle_row+1, particle_col +1 );
        }
        for( int particle_col = particle_row+1; particle_col < nparticles; ++particle_col)
        {
            auto v = particle_v[particle_col].Force( particle_v[particle_row] );
            forces_m.SetVal( particle_row +1, particle_col +1, v);
            f[particle_row] += v;
        }
    }


    return;
}


void TIntegrator::DoStep()
{
    ++step;

    std::vector<TVector> force_v( particle_v.size() );

    CalculateForce( force_v );

    TheIntegrator(force_v);

    return;
}

void TIntegrator::DoStepCXX()
{
    ++step;

    // first step of Verlet algorithm

    {
        auto ff = [&](TParticle & particle_i)
        {

            if( particle_i.isFixed ) return;
            const int id = particle_i.id;
            particle_i.force.x = 0;
            particle_i.force.y = 0;
            for( auto aux : particle_v)
            {
                if( id == aux.id ) continue;
                particle_i.force += aux.Force( particle_i );
            }
            particle_i.pos.Add( particle_i.vel, h);
            particle_i.pos.Add( particle_i.force, 0.5*h*h/particle_i.mass );
            return;
        };
        for_each(std::execution::par, particle_v.begin(), particle_v.end(), ff);

    }


    //recalculate forces at new positions, and sum force at the next position to the previous calculated force
    {

        auto ff = [&](TParticle & particle_i) {
            if( particle_i.isFixed ) return;
            const int id = particle_i.id;
            for( auto aux : particle_v)
            {
                if( id == aux.id ) continue;
                particle_i.force += aux.Force( particle_i );
            }
            // update Velocity
            particle_i.vel.Add( particle_i.force, 0.5*h/particle_i.mass );


        };

        for_each(std::execution::par, particle_v.begin(), particle_v.end(), ff);

    }

}
void TIntegrator::DoStepTBB()
{
    ++step;

    // first step of Verlet algorithm

    {
        auto ff = [&](TParticle & particle_i)
        {

            if( particle_i.isFixed ) return;
            const int id = particle_i.id;
//             particle_i.force.x = 0;
//             particle_i.force.y = 0;
            particle_i.force = particle_i.vel.Scale( this->fDamp );
//             for( auto aux : particle_v)
            for( auto & aux : particle_i.related_particles )
            {
                if( id == aux.get().id ) continue;
                particle_i.force += aux.get().Force( particle_i );
            }
            particle_i.pos.Add( particle_i.vel, h);
            particle_i.pos.Add( particle_i.force, 0.5*h*h/particle_i.mass );
            return;
        };

        tbb::parallel_for( tbb::blocked_range<int>(0,particle_v.size()),
                           [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                ff(particle_v[i]);
            }
        });

    }


    //recalculate forces at new positions, and sum force at the next position to the previous calculated force
    {

        auto ff = [&](TParticle & particle_i) {
            if( particle_i.isFixed ) return;
            const int id = particle_i.id;
//             for( auto aux : particle_v)
            for( auto & aux : particle_i.related_particles)
            {
                if( id == aux.get().id ) continue;
                particle_i.force += aux.get().Force( particle_i );
            }
            // update Velocity
            particle_i.vel.Add( particle_i.force, 0.5*h/particle_i.mass );


        };
        tbb::parallel_for( tbb::blocked_range<int>(0,particle_v.size()),
                           [&](tbb::blocked_range<int> r)
        {
            for (int i=r.begin(); i<r.end(); ++i)
            {
                ff(particle_v[i]);
            }
        });

    }

}

void TIntegrator::SetCriticalRadius(double r)
{
    if( r<0 )
        r = std::numeric_limits<double>::max();


    for( auto & p : particle_v )
    {
        for( int i = 0; i < particle_v.size(); ++i)
        {
            if( TVector::distance( p.pos , particle_v[i].pos ) < r && p.id!=particle_v[i].id )
                p.related_particles.push_back( std::reference_wrapper<TParticle>( particle_v[i]));
        }
// #ifndef NDEBUG
        std::cout << "Particle pos: " << p.pos << "\t Number of related particles: " << p.related_particles.size() << std::endl;
// #endif
    }
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
    return;
}

void TIntegrator::IntegratorVerletPar(std::vector<TVector>& force_v)
{

    // first step of Verlet algorithm

    {
        auto ff = [&](TParticle & particle_i, TVector & force_i)
        {

            if( particle_i.isFixed ) return force_i;
            particle_i.pos.Add( particle_i.vel, h);
            particle_i.pos.Add( force_i, 0.5*h*h/particle_i.mass );
            particle_i.force = force_i;
            return force_i;
        };
        std::transform( std::execution::par,
                        particle_v.begin(),
                        particle_v.end(),
                        force_v.begin(),
                        force_v.begin(),
                        ff
                      );
    }







    //recalculate forces at new positions, and sum force at the next position to the previous calculated force
    {

        auto ff = [&](TParticle & particle_i) {
            if( particle_i.isFixed ) return;
            const int id = particle_i.id;
            for( auto aux : particle_v)
            {
                if( id == aux.id ) continue;
                particle_i.force += aux.Force( particle_i );
            }
            // update Velocity
            particle_i.vel.Add( particle_i.force, 0.5*h/particle_i.mass );


        };

        for_each(std::execution::par, particle_v.begin(), particle_v.end(), ff);

    }
    return;
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

void TIntegrator::Test_benchmark_force_calculation(int i_forceMethod, int i_i_integratorMethod )
{
    TIntegrator myIntegrator;

    myIntegrator.h = 3600;

    TParticle::f_constant = G_UA_MSun;

    myIntegrator.SetNparticlesRnd(500);
    myIntegrator.SetForceCalSimple();

    switch(i_forceMethod)
    {
    case 1:
        myIntegrator.SetForceCalSimple();
        break;
    case 2:
        myIntegrator.SetForceCalSimplePar();
        break;
    case 3:
        myIntegrator.SetForceCalTriang();
        break;
    case 4:
        myIntegrator.SetForceCalTriangPar();
        break;
    case 5:
        myIntegrator.SetForceCalSimpleParCXX();
        break;
    default:
        std::cout << "Error while seting force method. Back to default" << std::endl;
    };
    switch(i_i_integratorMethod)
    {
    case 1:
        myIntegrator.SetIntegratorEulerFW();
        break;
    case 2:
        myIntegrator.SetIntegratorVerlet();
        break;
    case 3:
        myIntegrator.SetIntegratorVerletPar();
        break;
    default:
        std::cout << "Error while seting integrator method. Back to default" << std::endl;
    };


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

    if ( -1 == i_forceMethod && -1 == i_i_integratorMethod )
    {
        for( int i = 0; i < 500; ++i)
        {
            myIntegrator.DoStepCXX();
            if( 0 == i % 100 ) std::cout << "Step " << i << std::endl;
        }

    }
    else if( -2 == i_forceMethod && -2 == i_i_integratorMethod )
    {

        for( int i = 0; i < 500; ++i)
        {
            myIntegrator.DoStepTBB();
            if( 0 == i % 100 ) std::cout << "Step " << i << std::endl;
        }
    }
    else
    {

        for( int i = 0; i < 500; ++i)
        {
            myIntegrator.DoStep();
            if( 0 == i % 100 ) std::cout << "Step " << i << std::endl;
        }
    }

    return;
}
