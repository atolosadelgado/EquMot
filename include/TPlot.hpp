// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


#ifndef __TPLOT_HPP__
#define __TPLOT_HPP__

#include "TVector.hpp"
#include <string>

class TPlot
{
public:
    TPlot(std::string _s = "", bool _save_png = false);
    ~TPlot();
    void StartPlot(int n=0);
    void AddPoint( TVector & v);
    void ShowPlot();
    FILE * plotHandle;
    bool save_png = {false};
    std::string name = {""};
    void StartH2D(int n=0);
    void AddPointH2D(  TVector & v, double w );
};


#endif
