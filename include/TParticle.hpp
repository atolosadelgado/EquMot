#ifndef __TPARTICLE_HPP__
#define __TPARTICLE_HPP__

#include "TVector.hpp"

class TParticle
{
public:
    TParticle ();
    TVector pos;
    TVector vel;
    int id = {0};
    static int nparticles;
    static double f_constant;

    double Kene() const ;
    TVector Force( TVector & ipos);

    double Vene( TVector & r );
    double Mene( TVector & r );

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};





#endif
