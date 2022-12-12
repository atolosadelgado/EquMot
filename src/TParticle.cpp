// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TParticle.hpp"

TParticle::TParticle () {
    pos.SetRnd();
    vel.x = pos.x;
    vel.y = pos.y;
    id = (++nparticles);
};

double TParticle::Kene() const {
    return 0.5*vel.Norm2();
}

TVector TParticle::Force( TParticle & ipar ) {
    auto d = pos.distance(ipar.pos);
    d=pow(d,3);
    return TVector( -f_constant*ipar.mass*mass*(ipar.pos.x - pos.x)/d, -f_constant*ipar.mass*mass*(ipar.pos.y - pos.y)/d);
};

double TParticle::Vene( TParticle & ipar ) {
    return -f_constant*mass*ipar.mass/pos.distance(ipar.pos);
};
double TParticle::Mene( TParticle & ipar ) {
    return Kene()+Vene(ipar);
};



int TParticle::nparticles = 0;
double TParticle::f_constant = 1000.;




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

