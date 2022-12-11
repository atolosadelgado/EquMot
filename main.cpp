// Hello world

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include <thread>

unsigned seed = 1; // std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> rnd_position (-100.,100.);




class TVector
{
public:
    TVector() {};
    TVector(double ix, double iy):x(ix),y(iy) {};
    double x = {0.0};
    double y = {0.0};

    void Set(  TVector & v ) {
        x = v.x;
        y=v.y;
        return;
    };
    void SetRnd() {
        x=rnd_position(generator);
        y=rnd_position(generator);
    };

    double Norm() const {
        return std::sqrt( pow(x,2) + pow(y,2));
    }
    double Norm2() const {
        return  pow(x,2) + pow(y,2);
    }
    double distance( TVector & v) {
        return std::sqrt( pow(x -v.x,2) + pow(y-v.y,2));
    }
    static double distance( TVector & v,  TVector & vv) {
        return std::sqrt( pow(vv.x -v.x,2) + pow(vv.y-v.y,2));
    }

    TVector & operator * (double const & scale) {
        this->x*=scale;
        this->y*=scale;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const TVector& v);
    TVector& operator+=(const TVector& v) {

        this->x += v.x;
        this->y += v.y;
        return *this;
    }

};

std::ostream& operator<<(std::ostream& os, const TVector& v)
{
    os <<  '\t' << v.x << ' ' << v.y << std::endl; ;
    return os;
}

class TParticle
{
public:
    TParticle () {
        pos.SetRnd();
        vel.x = 0.01* pos.y;
        vel.y = -0.01*pos.x;
        id = (++nparticles);
    };
    TVector pos;
    TVector vel;
    int id = {0};
    static int nparticles;
    static double f_constant;

    double Kene() const {
        return 0.5*vel.Norm2();
    }
    TVector Force( TVector & ipos) {
        auto d = pos.distance(ipos);
        d=pow(d,3);
        return TVector( -f_constant*(ipos.x - pos.x)/d, -f_constant*(ipos.y - pos.y)/d);
    };

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};

int TParticle::nparticles = 0;
double TParticle::f_constant = 100.;

class myGNUplotter
{
public:
    myGNUplotter();
    ~myGNUplotter();
    void StartPlot();
    void AddPoint( TVector & v);
    void ShowPlot();
    FILE * plotHandle;
};

myGNUplotter::myGNUplotter()
{
    plotHandle = popen("gnuplot -persist", "w");

    if(plotHandle == NULL)
    {
        std::runtime_error("Can not open pipe to GNU plot");
    }

    fprintf(plotHandle,"set terminal qt\n");
    fprintf(plotHandle,"set xlabel \"x\"\n");
    fprintf(plotHandle,"set ylabel \"y\"\n");
    fprintf(plotHandle,"set xrange [-100:100]\n");
    fprintf(plotHandle,"set yrange [-100:100]\n");
    fprintf(plotHandle,"set grid\n");
}

void myGNUplotter::StartPlot()
{
    fprintf( plotHandle, "plot '-'\n");
}


void myGNUplotter::AddPoint(TVector& v)
{
    fprintf( plotHandle, "%g %g\n", v.x, v.y);
}

void myGNUplotter::ShowPlot()
{
    fprintf( plotHandle, "e\n");
    fflush(plotHandle);
}

myGNUplotter::~myGNUplotter()
{
    fprintf( plotHandle, "exit\n");
    pclose(plotHandle);
}





std::ostream& operator<<(std::ostream& os, const TParticle& v)
{
    os << "Particle #" << v.id << " has this Kinetic energy " << v.Kene() << " and it is placed at " ;
    os << v.pos;
    return os;
}

int main() {
    std::cout << "Hello world!\n";

    TParticle a;
    TParticle b;

    double h = 1;
    myGNUplotter plt;

    for( int i = 0; i < 100; ++i)
    {
        std::cout << a << std::endl;
        std::cout << b << std::endl;


        std::this_thread::sleep_for(std::chrono::milliseconds(100));
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




    std::cout << a << std::endl;
    std::cout << b << std::endl;
    std::cout << a.Force(b.pos) << std::endl;

    std::cout << (a.pos+=b.pos) << std::endl;




}
