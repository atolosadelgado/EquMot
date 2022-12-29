// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TPARTICLE_HPP__
#define __TPARTICLE_HPP__

#include "TVector.hpp"

#include <functional>
#include <numeric>

/// Class to store properties of a particle
class TParticle
{
public:
    TParticle ();

    /// position
    TVector pos;

    /// Vecocity
    TVector vel;

    /// Mass
    double mass = {1.0};

    /// Charge
    double charge = {0.0};

    /// force
    TVector force;

    /// If fixed, position will not be updated
    bool isFixed = {false};

    /// Id number of the particle
    int id = {0};

    /// Static counter of particles
    static int nparticles;

    /// Constant of the force
    static double f_constant;

    /// Constant of the electric force
    static double fE_constant;

    /// Constant of the gravitational force
    static double fG_constant;

    /// Constant of the elastic force
    static double fL_constant;

    /// Compute force only from particles inside the radius frad_critical
    static double frad_critical;

    /// Compute Kinetic energy
    double Kene() const ;

    /// Compute force that would be applied to another particle
    std::function< TVector(TParticle& ipar)> Force;

    void SetForce(int n);

    /// Compute the Electric force that would be applied to another particle
    inline TVector ForceE(TParticle& ipar);

    /// Compute the Electric force that would be applied to another particle
    inline TVector ForceG(TParticle& ipar);

    /// Compute the elasticforce that would be applied to another particle
    inline TVector ForceL(TParticle& ipar);

    /// Compute LJ force
    inline static TVector ForceLJ( TParticle & a, TParticle & b);

    /// Compute the potential that other particle r feels
    double Vene( TParticle & r );

    /// To be deprecated
    double Mene( TParticle & r );

    /// Set fix position, do not update the position/velocity
    void SetFixPosition();

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};

inline TVector TParticle::ForceE( TParticle & ipar ) {
//     auto d = pos.distance(ipar.pos);
//     d=pow(d,3);
//     return TVector( -f_constant*ipar.mass*mass*(ipar.pos.x - pos.x)/d, -f_constant*ipar.mass*mass*(ipar.pos.y - pos.y)/d);
//     auto d = pos.distance(ipar.pos);
//     d=pow(d,3);
    double dx = ipar.pos.x - pos.x;
    double dy = ipar.pos.y - pos.y;
    double d2= dx*dx + dy*dy;
    double d = sqrt(d2);
    double fE_mod = 0;
    if( d < frad_critical )
        fE_mod = -fE_constant*ipar.charge*charge/d2;
    return TVector( fE_mod*(dx)/d, fE_mod*(dy)/d);
};

inline TVector TParticle::ForceG( TParticle & ipar ) {
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

inline TVector TParticle::ForceL(TParticle& ipar)
{
    double dx = ipar.pos.x - pos.x;
    double dy = ipar.pos.y - pos.y;
    double d2= dx*dx + dy*dy;
    double d = sqrt(d2);
    double fL_mod = 0;
    if( d < frad_critical )
        fL_mod = -fL_constant*d;
    return TVector( fL_mod*(dx)/d, fL_mod*(dy)/d);
}


#endif
