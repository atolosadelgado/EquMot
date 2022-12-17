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


    /// If fixed, position will not be updated
    bool isFixed = {false};

    /// Id number of the particle
    int id = {0};

    /// Static counter of particles
    static int nparticles;

    /// Constant of the force
    static double f_constant;

    /// Compute Kinetic energy
    double Kene() const ;

    /// Compute force that would be applied to another particle
    inline TVector Force(TParticle& ipar);

    /// Compute the potential that other particle r feels
    double Vene( TParticle & r );

    /// To be deprecated
    double Mene( TParticle & r );

    /// Set fix position, do not update the position/velocity
    void SetFixPosition();

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};

inline TVector TParticle::Force( TParticle & ipar ) {
//     auto d = pos.distance(ipar.pos);
//     d=pow(d,3);
//     return TVector( -f_constant*ipar.mass*mass*(ipar.pos.x - pos.x)/d, -f_constant*ipar.mass*mass*(ipar.pos.y - pos.y)/d);
//     auto d = pos.distance(ipar.pos);
//     d=pow(d,3);
    double dx = ipar.pos.x - pos.x;
    double dy = ipar.pos.y - pos.y;
    double d32= pow( dx*dx + dy*dy, 1.5);
    return TVector( -f_constant*ipar.mass*mass*(dx)/d32, -f_constant*ipar.mass*mass*(dy)/d32);
};



#endif
