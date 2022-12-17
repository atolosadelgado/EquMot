// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TMatrixtriang_HPP__
#define __TMatrixtriang_HPP__

#include <vector>
#include <iostream>

/// Class for storing a symetric matrix, sij = sji, using a vector as container
class TMatrixtriang
{
public:
    TMatrixtriang();
    TMatrixtriang(int n);
    void SetSize(int n);
    void SetVal( int row, int column, double value);
    double GetVal( int row, int column) const ;
    friend std::ostream& operator<<(std::ostream& os, const TMatrixtriang& v);

    std::vector<double> matrix;
    int _n;



};

#endif
