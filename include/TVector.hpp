// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TVECTOR_HPP__
#define __TVECTOR_HPP__


#include <cmath>
#include <random>
#include <chrono>
#include <iostream>




class TVector
{
public:
    TVector() {};
    TVector(double ix, double iy);
    double x = {0.0};
    double y = {0.0};

    void Set(  TVector & v );

    void SetRnd();

    double Norm() const;
    double Norm2() const;
    double distance( TVector & v);
    static double distance( TVector & v,  TVector & vv);

    TVector Scale(double h);

//     TVector & operator * (double const & scale);

    friend std::ostream& operator<<(std::ostream& os, const TVector& v);
    TVector& operator+=(const TVector& v);

};

#endif
