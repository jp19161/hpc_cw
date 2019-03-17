#ifndef CLASS_BURGERS
#define CLASS_BURGERS
#include <cmath>
#include "Model.h"
#include <mpi.h>

using namespace std;

class Burgers{
    
public:
    Burgers(Model& m);
    ~Burgers();

    void Initial_velocity();
    void Integrate_velocity();
    void Print_velocity();
    void Energy_Calculation();
    void Communication();
    void AssembleGrid();
    void PatchUpU();
    void PatchUpV();
    
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
       
       //MPI parameters
        int myrank, mysize, retval_rank, retval_size;
        int mpicols;
        int mpirows;
        int *width_process_u;
        int *height_process_u;
        int *width_process_v;
        int *height_process_v;
        int local_nx;
        int local_ny;
        int px;
        int py;
        int my_pos_x = 0;
        int my_pos_y= 0;
        double *local_u;
        double *local_v;
        double *local_x;
        double *local_y;
        double *local_u_new;
        double *local_v_new;
        double local_u_next;
        double local_v_next;
        double* bigu;
        double* bigv;
        double* BigU;
        double* BigV;

        
        // communication arrays    
         double *outer_left_u;
         double *outer_right_u;
         double *outer_up_u;
         double *outer_down_u;
         double *inner_left_u;
         double *inner_right_u;
         double *inner_up_u;
         double *inner_down_u;
         
};




#endif /* Burgers_h */