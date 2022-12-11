#include "TVector.hpp"


unsigned seed = 1; // std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> rnd_position (-100.,100.);


TVector::TVector(double ix, double iy):x(ix),y(iy) {};

void TVector::Set(  TVector & v ) {
    x = v.x;
    y=v.y;
    return;
};


void TVector::SetRnd() {
    x=rnd_position(generator);
    y=rnd_position(generator);
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

TVector & TVector::operator * (double const & scale) {
    this->x*=scale;
    this->y*=scale;
    return *this;
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
