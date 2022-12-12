// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TPARTICLE_HPP__
#define __TPARTICLE_HPP__

#include "TVector.hpp"

/// Class to store properties of a particle
class TParticle
{
public:
    TParticle ();

    /// position
    TVector pos;

    /// Vecocity
    TVector vel;

    /// Mass/charge
    double mass = {1.0};

    /// Id number of the particle
    int id = {0};

    /// Static counter of particles
    static int nparticles;

    /// Constant of the force
    static double f_constant;

    /// Compute Kinetic energy
    double Kene() const ;

    /// Compute force that would be applied to another particle
    TVector Force( TParticle & ipos);

    /// Compute the potential that other particle r feels
    double Vene( TParticle & r );

    /// To be deprecated
    double Mene( TParticle & r );

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};





#endif
