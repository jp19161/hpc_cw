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
    MPI_Init(&argc, &argv);
    
    Model m(argc, argv);
    //cout<<"Created model"<<endl;
    Burgers b(m);  
    //cout<<"Created burgers"<<endl;

    b.Initial_velocity();
    MPI_Barrier(MPI_COMM_WORLD);
    //cout<<"Got here"<<endl;
    //b.~Burgers();
    //cout<<"Here as wee"<<endl;
   b.PatchUpU();
   b.PatchUpV();
//    b.Integrate_velocity();
//    b.Energy_Calculation();

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    hrc::time_point end = hrc::now();


    MPI_Barrier(MPI_COMM_WORLD);
    cout<<"Finished"<<endl;

    MPI_Finalize();
    return 0;
}
