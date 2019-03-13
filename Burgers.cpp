#include <chrono>
#include "Model.h"
#include "Burgers.h"
#include <iostream>
#include <stdio.h>
using namespace std;

Burgers::Burgers(Model& m) 
{
    // need to calculate coefficients for integration
    
    // parameters needed for this sections funcitons 
    ax = m.GetAx();
    ay = m.GetAy();
    b = m.GetB();
    c= m.GetC();
    L = m.GetL();
    Nx = m.GetNx();
    Ny = m.GetNy();
    Nt = m.GetNt();
    T = m.GetT();
    dx = m.GetDx();
    dy = m.GetDy();
    dt = m.GetDt();
}
    



// Initialise velocity
void Burgers::Initial_velocity()
{
    double r, rb;
    auto* u = new double[Nx * Ny];
    auto* v = new double[Nx * Ny];
    double x[Nx];
    double y[Ny];

    for(int i = 0; i < Nx; i++) {
        x[i] = -L / 2 + i * dx;
    }
    for(int i = 0; i < Ny; i++) {
        y[i] = -L / 2 + i * dy;
    }
    for(int col = 0; col < Nx; ++col) {
        for(int row = 0; row < Ny; ++row) {
            r = sqrt(pow(x[col], 2) + pow(y[row], 2));
            if(r <= 1) {
                rb = 2 * pow(1 - r, 4) * (4 * r + 1);
                // adjustBounds(row,col);
            } else {
                rb = 0.0;
            }
            u[col * Ny + row] = rb;
            v[col * Ny + row] = rb;
        }
    }

    for(int i = 0; i < Nx; i++) {
        u[i * Ny + 0] = 0.0;
        u[i * Ny + Ny - 1] = 0.0;
        v[i * Ny + 0] = 0.0;
        v[i * Ny + Ny - 1] = 0.0;
    }
    for(int i = 0; i < Ny; i++) {
        u[0 * Ny + i] = 0.0;
        u[(Nx - 1) * Ny + i] = 0.0;
        v[0 * Ny + i] = 0.0;
        v[(Nx - 1) * Ny + i] = 0.0;
    }
}

// Integrate Velocity
void Burgers::Integrate_velocity()
{
    auto* u_new = new double[Nx * Ny];
    auto* v_new = new double[Nx * Ny];
    // auto* u_next = new double[Nx * Ny];
    // auto* v_next = new double[Nx * Ny];
    const double c1 = c / dx / dx + ax / dx;
    const double c2 = c / dx / dx;
    const double c3 = -2.0 * c * (1 / dx / dx + 1 / dy / dy) - ax / dx - ay / dy + 1 / dt;
    const double c4 = c / dy / dy + ay / dy;
    const double c5 = c / dy / dy;

    //for(int k = 1; k < Nt; k++) 
        for(int i = 1; i < Nx - 1; i++) {
            for(int j = 1; j < Ny - 1; j++) {
                u_next = c1 * u[(i - 1) * Ny + j] + c2 * u[(i + 1) * Ny + j] + c3 * u[i * Ny + j] + c4 * u[i * Ny + j - 1] + c5 * u[i * Ny + j + 1] +b / dx * u[i * Ny + j] * (-u[i * Ny + j] + u[(i - 1) * Ny + j]) +  b / dy * v[i * Ny + j] * (u[i * Ny + j - 1] - u[i * Ny + j]);
                v_next = c1 * v[(i - 1) * Ny + j] + c2 * v[(i + 1) * Ny + j] + c3 * v[i * Ny + j] +  c4 * v[i * Ny + j - 1] + c5 * v[i * Ny + j + 1] + b / dx * v[i * Ny + j] * (-v[i * Ny + j] + v[(i - 1) * Ny + j]) +  b / dy * u[i * Ny + j] * (v[i * Ny + j - 1] - v[i * Ny + j]);
                u_new[i * Ny + j] = dt * u_next;
                v_new[i * Ny + j] = dt * v_next;
            }
        }

        // Reassigning values to the velocity arrays
        for(int in = 1; in < Nx - 1; in++) {
            for(int jn = 1; jn < Ny - 1; jn++) {
                u[in * Ny + jn] = u_new[in * Ny + jn];
                v[in * Ny + jn] = v_new[in * Ny + jn];
            }
        }
    
}


//writing velocity field to file 






void Burgers::Energy_Calculation()
{
    double energy;
    for(int i = 0; i < Nx * Ny; i++) {
        energy += (pow(u[i],2) + pow(v[i],2));
    }
    energy *= 0.5 * dx * dy;
    cout << energy << endl;
}



// Model m(argc, argv);
// Burgers b(&m);
//b.Initial_velocity();

    
/*    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Call code to perform time integration here

    hrc::time_point end = hrc::now();

    // Calculate final energy and write output

    return 0;*/



