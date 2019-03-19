#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdio.h>
//#include <mpi.h>
#include "Burgers.h"
#include "Model.h"
using namespace std;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    Model m(argc, argv);
//cout << " model m done" << endl;
    Burgers b(m);
//cout << " burger b done" << endl;

    b.Initial_velocity();
    cout << " initial velocity done" << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    b.Integrate_velocity();
    cout << " integrate velcoity done" << endl;

    hrc::time_point end = hrc::now();

    b.PatchUpU();
cout << " U patch done" << endl;

    b.PatchUpV();
cout << " V patch done" << endl;

    b.Print_velocity();
cout << " Printed verlocityy done" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    b.Energy_Calculation();
cout << " calculated energy done" << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    auto diff = end - start;
    cout << "\n\nexexution time: " << chrono::duration<double, milli>(diff).count() << " ms" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
