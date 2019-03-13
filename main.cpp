#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdio.h>
//#include <mpi.h>
#include "Model.h"
#include "Burgers.h"
using namespace std;

int main(int argc, char* argv[])
{
    //MPI_Init(&argc, &argv);
    
    Model m(argc, argv);
    Burgers b(m);  
    /*
     b.Initial_Velocity();
      b.Integrate_Velocity();
       b.Energy_Calculation();
    */
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    hrc::time_point end = hrc::now();

   // MPI_Finalize();
    return 0;
}
