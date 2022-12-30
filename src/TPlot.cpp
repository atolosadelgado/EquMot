// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include "TPlot.hpp"


TPlot::TPlot()
{

#ifdef GNUPLOT
    plotHandle = popen("gnuplot -persist", "w");

    if(plotHandle == NULL)
    {
        throw std::runtime_error("Can not open pipe to GNU plot");
    }

    fprintf(plotHandle,"set terminal qt\n");
    fprintf(plotHandle,"set xlabel \"x\"\n");
    fprintf(plotHandle,"set ylabel \"y\"\n");
    fprintf(plotHandle,"set xrange [-5:25]\n");
    fprintf(plotHandle,"set yrange [-5:25]\n");
    fprintf(plotHandle,"set grid\n");
    fprintf(plotHandle, "set style fill solid 1.0\n");
#endif

}

void TPlot::StartPlot()
{
#ifdef GNUPLOT
//     fprintf( plotHandle, "plot '-'\n");
    fprintf( plotHandle, "plot '-' w circles lc 1 fs solid\n");
#endif

}


void TPlot::AddPoint(TVector& v)
{
#ifdef GNUPLOT
    fprintf( plotHandle, "%g %g 0.05\n", v.x, v.y);
#endif
}

void TPlot::ShowPlot()
{
#ifdef GNUPLOT
    fprintf( plotHandle, "e\n");
    fflush(plotHandle);
#endif
}


void TPlot::StartH2D()
{
#ifdef GNUPLOT
//     fprintf( plotHandle, "plot '-'\n");
    fprintf( plotHandle, "plot '-' w boxxy fc palette z\n");
#endif

}


void TPlot::AddPointH2D(TVector& v, double w)
{
#ifdef GNUPLOT
    fprintf( plotHandle, "%g %g 0.5 0.5 %g\n", v.x, v.y, w);
#endif
}




TPlot::~TPlot()
{
#ifdef GNUPLOT
    fprintf( plotHandle, "exit\n");
    pclose(plotHandle);
#endif
}

