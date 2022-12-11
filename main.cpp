// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// //  Alvaro Tolosa Delgado 2022
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

#include <iostream>

#include <thread>

#include "TVector.hpp"

#include "TParticle.hpp"

#include "TPlot.hpp"




int main() {

    TParticle a, b;

    double h = 1.;
    TPlot plt;

    for( int i = 0; i < 1000; ++i)
    {
        std::cout << "Step " << i << "\t" << a.Mene(b.pos) + b.Mene(a.pos) << std::endl;


        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        plt.StartPlot();

        plt.AddPoint( a.pos );
        plt.AddPoint( b.pos );


        plt.ShowPlot();


        TVector forcea = a.Force(b.pos);
        TVector forceb = b.Force(a.pos);

        a.pos += a.vel*h;
        a.pos += forceb*(0.5*h*h);

        b.pos += b.vel*h;
        b.pos += forcea*(0.5*h*h);

        a.vel+= forceb;
        b.vel+=forcea;
        if( 1> a.pos.distance(b.pos) )
            break;
    }









}
