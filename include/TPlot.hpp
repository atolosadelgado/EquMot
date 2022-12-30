// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TPLOT_HPP__
#define __TPLOT_HPP__

#include "TVector.hpp"

class TPlot
{
public:
    TPlot();
    ~TPlot();
    void StartPlot();
    void AddPoint( TVector & v);
    void ShowPlot();
    FILE * plotHandle;
    void StartH2D();
    void AddPointH2D(  TVector & v, double w );
};


#endif
