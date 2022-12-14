// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TVECTOR_HPP__
#define __TVECTOR_HPP__


#include <cmath>
#include <random>
#include <chrono>
#include <iostream>



/// Implementation of 2D vector in cartesian coordinates
class TVector
{
public:
    /// Default constructor (0,0)
    TVector() {};

    /// Constructor which takes coordinates as input
    TVector(double ix, double iy);

    /// X component of the vector
    double x = {0.0};
    /// Y component of the vector
    double y = {0.0};

    /// Set components from another vector
    void Set(  TVector & v );

    /// Set components randomly
    void SetRnd();

    /// return the euclidean norm, sqrt( x**2 + y**2)
    double Norm() const;

    /// return the square euclidean norm, ( x**2 + y**2)
    double Norm2() const;

    /// return the euclidean distance, sqrt( dx**2 + dy**2)
    double distance( TVector & v);

    /// return the euclidean distance, sqrt( dx**2 + dy**2)
    static double distance( TVector & v,  TVector & vv);

    TVector operator * (double const & h);

    /// Return vector scaled by factor h
    TVector Scale(double h);

    /// Add vector v scaled by factor h
    void Add(TVector & v, double h);

    friend std::ostream& operator<<(std::ostream& os, const TVector& v);

    /// Custom operator to sum another vector to this vector
    TVector& operator+=(const TVector& v);

};

#endif
