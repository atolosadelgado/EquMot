// Hello world

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

unsigned seed = 1; // std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::uniform_real_distribution<double> rnd_position (-100.,100.);

class TVector
{
public:
    TVector(){};
    TVector(double ix, double iy):x(ix),y(iy){};
    double x = {0.0};
    double y = {0.0};

    void Set(  TVector & v ){ x = v.x; y=v.y; return; };
    void SetRnd(){x=rnd_position(generator); y=rnd_position(generator);};

    double Norm() const { return std::sqrt( pow(x,2) + pow(y,2)); }
    double Norm2() const { return  pow(x,2) + pow(y,2); }
    double distance( TVector & v){ return std::sqrt( pow(x -v.x,2) + pow(y-v.y,2)); }
    static double distance( TVector & v,  TVector & vv){ return std::sqrt( pow(vv.x -v.x,2) + pow(vv.y-v.y,2)); }

    TVector & operator * (double const & scale) {
        this->x*=scale;
        this->y*=scale;
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const TVector& v);
    TVector& operator+=(const TVector& v){

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
    TParticle (){ pos.SetRnd(); vel.x = 0.1* pos.y; vel.y = -0.1*pos.x; id = (++nparticles); };
    TVector pos;
    TVector vel;
    int id = {0};
    static int nparticles;

    double Kene() const {return 0.5*vel.Norm2(); }
    TVector Force( TVector & ipos){ auto d = pos.distance(ipos); d=pow(d,3); return TVector( -(ipos.x - pos.x)/d, -(ipos.y - pos.y)/d); };

    friend std::ostream& operator<<(std::ostream& os, const TParticle& v);
};

int TParticle::nparticles = 0;

std::ostream& operator<<(std::ostream& os, const TParticle& v)
{
    os << "Particle #" << v.id << " has this Kinetic energy " << v.Kene() << " and it is placed at " ;
    os << v.pos;
    return os;
}

int main(){
    std::cout << "Hello world!\n";

    TParticle a;
    TParticle b;

    double h = 1;

    for( int i = 0; i < 100; ++i)
    {
        std::cout << a << std::endl;
        std::cout << b << std::endl;

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

  FILE *plotHandle = popen("gnuplot -persist", "w");

  if(plotHandle == NULL){
   return -1;
  }

  fprintf(plotHandle,"set terminal qt\n");
  fprintf(plotHandle,"set xlabel \"frequency\"\n");
  fprintf(plotHandle,"set ylabel \"Amplitude\"\n");
  fprintf(plotHandle,"set xrange [-3:3]\n");
  fprintf(plotHandle,"set grid\n");
  fprintf(plotHandle,"plot (x**2)*sin(x*100)\n");
// plot '-' w p ls 1, '-' w p ls 2, '-' w p ls 3
// 1 2
// e
// 2 1
// e

  fflush(plotHandle);

  pclose(plotHandle);

}
