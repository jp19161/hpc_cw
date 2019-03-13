#ifndef CLASS_BURGERS
#define CLASS_BURGERS
#include <cmath>
#include "Model.h"
using namespace std;

class Burgers{
    
public:
    Burgers(Model& m);
    ~Burgers(){};

    void Initial_velocity();
    void Integrate_velocity();
    void Print_velocity();
    void Energy_Calculation();
    
private:
// parameters will go here
// call to functoins also 
    Model* m;
 

        // Numerics
        double x0;
        double y0;
        double L;
        double Lx;
        double Ly;
        double T ;
        int Nx ;
        int Ny ;
        int Nt ;
        double dx;
        double dy;
        double dt;

        // Physics
        double ax;
        double ay;
        double b;
        double c;
        
        // Parameters needed by burgers.cpp functions 
       double *u;
       double *v;
       double *u_new;
       double *v_new;
       double u_next;
       double v_next;
       double *x;
       double *y;
};




#endif /* Burgers_h */