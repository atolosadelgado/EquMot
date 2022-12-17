// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TMatrixtriang.hpp"

TMatrixtriang::TMatrixtriang():matrix(1, 0),_n(1){}

/// Constructor to initialize the vector that correspond to a matrix nxn
/// The number of non identical elements corresponds to the Gaus summation of the first n numbers
TMatrixtriang::TMatrixtriang(int n):matrix(n*(n+1)/2 +1, 0),_n(n){}

/// Resize the current matrix using std::vector::resize method
void TMatrixtriang::SetSize(int n)
{
    matrix.resize( n*(n+1)/2 +1 );
    _n = n;
    return;
}

/// Get Value from row,column position
double TMatrixtriang::GetVal (int row, int column) const
{
    if( column < row ) std::swap(column, row);
    return matrix[ column + (row-1)*( 2*_n-row)/2];
}

/// Set Value for row, column position
void TMatrixtriang::SetVal(int row, int column, double value)
{
    if( column < row ) std::swap(column, row);
    matrix[ column + (row-1)*( 2*_n-row)/2] = value;
    return;
}

std::ostream& operator<<(std::ostream& os, const TMatrixtriang& v)
{

    for( int row = 1; row <= v._n; ++row )
    {
        std::cout << '(';
        for( int column = 1; column  <= v._n; ++column  )
        {
            std::cout << v.GetVal( row, column ) << '\t';
        }
        std::cout << ')' << std::endl;

    }

    std::cout << std::endl;
    return os;
}
