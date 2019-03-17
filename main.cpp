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
    //
    //cout<<"Here as wee"<<endl;
    
  // b.PatchUpU();
   //cout << "done patch U" << endl;
  // b.PatchUpV();
   //cout << "got here" << endl;
//   b.Print_velocity();
   //cout << "i printed velocity" << endl;
   b.Communication();
    // b.Integrate_velocity();
    cout << "COMMUNICATION FINISHED" << endl;
//    b.Energy_Calculation();

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    hrc::time_point end = hrc::now();

    MPI_Barrier(MPI_COMM_WORLD);
   cout<< "Passed the BARRIER" << endl;
 b.~Burgers();
 cout << "did i get here though???????" << endl;
    MPI_Finalize();
    cout<<"Finished"<<endl;
    return 0;
}
