// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#include "TVector.hpp"


/// Fixed seed to keep reproducibility
unsigned seed = 1; // std::chrono::system_clock::now().time_since_epoch().count();

/// Random generator
std::default_random_engine generator (seed);

/// Random distribution to place the particles
std::uniform_real_distribution<double> rnd_position (-1.,1.);


TVector::TVector(double ix, double iy):x(ix),y(iy) {};

void TVector::Set(  TVector & v ) {
    x = v.x;
    y=v.y;
    return;
};


void TVector::SetRnd() {
    x=rnd_position(generator);
    y=rnd_position(generator);
    return;
};


double TVector::Norm() const {
    return std::sqrt( pow(x,2) + pow(y,2));
}
double TVector::Norm2() const {
    return  pow(x,2) + pow(y,2);
}
double TVector::distance( TVector & v) {
    return std::sqrt( pow(x -v.x,2) + pow(y-v.y,2));
}

double TVector::distance( TVector & v,  TVector & vv) {
    return std::sqrt( pow(vv.x -v.x,2) + pow(vv.y-v.y,2));
}

TVector TVector::operator * (double const & h) {
    return TVector( this->x*h, this->y*h );
}

TVector TVector::Scale(double h)
{
    return TVector( this->x*h, this->y*h );
}

void TVector::Add(TVector& v, double h)
{
    this->x += h*v.x;
    this->y += h*v.y;
    return;
}


TVector& TVector::operator+=(const TVector& v) {

    this->x += v.x;
    this->y += v.y;
    return *this;
}

std::ostream& operator<<(std::ostream& os, const TVector& v)
{
    os <<  '\t' << v.x << ' ' << v.y << std::endl; ;
    return os;
}

void TVector_Test1()
{
    TVector a(1.,0.);
    std::cout << (a*3.5).Norm() << std::endl;
    return;
}
