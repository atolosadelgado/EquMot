// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TMatrixtriang_HPP__
#define __TMatrixtriang_HPP__

#include <vector>
#include <iostream>

#include "TVector.hpp"

/// Class for storing a symetric matrix, sij = sji, using a vector as container
class TMatrixtriang
{
public:
    TMatrixtriang();
    TMatrixtriang(int n);
    void SetSize(int n);
    TVector & at (int row, int column);
    void SetVal( int row, int column, TVector & value);
    TVector GetVal( int row, int column) const ;
    friend std::ostream& operator<<(std::ostream& os, const TMatrixtriang& v);

    std::vector<TVector> matrix;
    int _n;



};

#endif
