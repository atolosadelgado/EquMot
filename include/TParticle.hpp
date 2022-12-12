// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TPARTICLE_HPP__
#define __TPARTICLE_HPP__

#include "TVector.hpp"

class TParticle
{
public:
    TParticle ();
    TVector pos;
    TVector vel;
    double mass = {1.0};
    int id = {0};
    static int nparticles;
    static double f_constant;

    double Kene() const ;
    TVector Force( TParticle & ipos);

    double Vene( TParticle & r );
    double Mene( TParticle & r );

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};





#endif
