// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TParticle.hpp"
#include "Contants.hpp"

TParticle::TParticle () {
    pos.SetRnd();
    vel.x = 0;
    vel.y = 0;
    id = (++nparticles);
    Force = [&](TParticle & p){return this->ForceG(p); };
};

TVector TParticle::ForceLJ(TParticle& a, TParticle& ipar)
{
    double dx = ipar.pos.x - a.pos.x;
    double dy = ipar.pos.y - a.pos.y;
    double r= sqrt( dx*dx + dy*dy);
    double sigma = 1.0; //nm
    double epsilon = 1.38e-5; // kg*nm2/s2
    double sigma_r = sigma/r;
    double LJ_mod = 0;
    if( r < frad_critical )
        LJ_mod = 24*epsilon/sigma*( 2*pow( sigma_r,13) -  pow( sigma_r,7) );

    return TVector( LJ_mod*(dx)/r, LJ_mod*(dy)/r);
}


void TParticle::SetForce(int n)
{
    if( 1 == n)
            Force = [&](TParticle & p){return this->ForceG(p); };
    if( 2 == n)
            Force = [&](TParticle & p){return this->ForceE(p); };
    if( 12 == n)
            Force = [&](TParticle & p){return (this->ForceE(p)+=this->ForceG(p)) ; };
    if( 3 == n)
            Force = std::bind( TParticle::ForceLJ, *this, std::placeholders::_1);
    if( 4 == n)
            Force = [&](TParticle & p){return this->ForceL(p); };

}


double TParticle::Kene() const {
    return 0.5*mass*vel.Norm2();
}

// inline TVector TParticle::Force( TParticle & ipar ) {
// //     auto d = pos.distance(ipar.pos);
// //     d=pow(d,3);
// //     return TVector( -f_constant*ipar.mass*mass*(ipar.pos.x - pos.x)/d, -f_constant*ipar.mass*mass*(ipar.pos.y - pos.y)/d);
// //     auto d = pos.distance(ipar.pos);
// //     d=pow(d,3);
//     double dx = ipar.pos.x - pos.x;
//     double dy = ipar.pos.y - pos.y;
//     double d32= pow( dx*dx + dy*dy, 1.5);
//     return TVector( -f_constant*ipar.mass*mass*(dx)/d32, -f_constant*ipar.mass*mass*(dy)/d32);
// };

double TParticle::Vene( TParticle & ipar ) {
//     return -f_constant*mass*ipar.mass/pos.distance(ipar.pos);
    return -fE_constant*charge*ipar.charge/pos.distance(ipar.pos);
};
double TParticle::Mene( TParticle & ipar ) {
    return Kene()+Vene(ipar);
};

void TParticle::SetFixPosition()
{
    isFixed = true;
    return;
}



int TParticle::nparticles = 0;
double TParticle::f_constant = 1000.;
double TParticle::fL_constant = 1.;
double TParticle::fE_constant = K_SI_e;
double TParticle::fG_constant = G_SI;
double TParticle::frad_critical = std::numeric_limits<double>::max();



std::ostream& operator<<(std::ostream& os, const TParticle& v)
{
    os << "Part #" << v.id
    << "\tKene " << v.Kene()
    << "\tPos " ;
    os << v.pos
    << "\tVel " ;
    os << v.vel;
    return os;
}

void TParticle_Test1()
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    TParticle a;
    std::cout << a << std::endl;
}

void TParticle_Test2()
{
    std::cout << __PRETTY_FUNCTION__ << std::endl;
    TParticle a, b;
    std::cout << a.Force(b) << std::endl;
    std::cout << (a.pos+=b.pos) << std::endl;

}

