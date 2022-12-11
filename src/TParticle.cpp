// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TParticle.hpp"

TParticle::TParticle () {
    pos.SetRnd();
    vel.x = 0.01* pos.y;
    vel.y = -0.01*pos.x;
    id = (++nparticles);
};

double TParticle::Kene() const {
    return 0.5*vel.Norm2();
}

TVector TParticle::Force( TVector & ipos) {
    auto d = pos.distance(ipos);
    d=pow(d,3);
    return TVector( -f_constant*(ipos.x - pos.x)/d, -f_constant*(ipos.y - pos.y)/d);
};

double TParticle::Vene( TVector & r ) {
    return f_constant/pos.distance(r);
};
double TParticle::Mene( TVector & r ) {
    return Kene()+Vene(r);
};



int TParticle::nparticles = 0;
double TParticle::f_constant = 100.;




std::ostream& operator<<(std::ostream& os, const TParticle& v)
{
    os << "Particle #" << v.id << " has this Kinetic energy " << v.Kene() << " and it is placed at " ;
    os << v.pos;
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
    std::cout << a.Force(b.pos) << std::endl;
    std::cout << (a.pos+=b.pos) << std::endl;

}

