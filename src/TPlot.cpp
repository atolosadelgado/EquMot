// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TPlot.hpp"


TPlot::TPlot()
{
    plotHandle = popen("gnuplot -persist", "w");

    if(plotHandle == NULL)
    {
        std::runtime_error("Can not open pipe to GNU plot");
    }

    fprintf(plotHandle,"set terminal qt\n");
    fprintf(plotHandle,"set xlabel \"x\"\n");
    fprintf(plotHandle,"set ylabel \"y\"\n");
    fprintf(plotHandle,"set xrange [-15:15]\n");
    fprintf(plotHandle,"set yrange [-15:15]\n");
    fprintf(plotHandle,"set grid\n");
}

void TPlot::StartPlot()
{
//     fprintf( plotHandle, "plot '-'\n");
    fprintf( plotHandle, "plot '-' w circles lc 1 fs solid\n");
}


void TPlot::AddPoint(TVector& v)
{
    fprintf( plotHandle, "%g %g 0.05\n", v.x, v.y);
}

void TPlot::ShowPlot()
{
    fprintf( plotHandle, "e\n");
    fflush(plotHandle);
}

TPlot::~TPlot()
{
    fprintf( plotHandle, "exit\n");
    pclose(plotHandle);
}

