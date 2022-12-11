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
};


#endif
