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
    //cout << "DONE INITIAL VELOCITY" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
//    b.Communication();
//    cout << "DONE COMMUNICATION" << endl;
//    b.MakeBigger();
//    cout << "DONE MAKE BIGGER" << endl;
    b.Integrate_velocity();
//    MPI_Barrier(MPI_COMM_WORLD);
//    cout << "DONE INTEGRATE VELOCITY" << endl;
//    MPI_Barrier(MPI_COMM_WORLD);
    b.PatchUpU();
  //  cout << " PATCHED UP U" << endl;
//     MPI_Barrier(MPI_COMM_WORLD);
    b.PatchUpV();
   //cout << " PATCHED UP V" << endl;
//     
   b.Print_velocity();
 //   cout << " I PRINTED VELOCITY" << endl;

    // cout << "i printed velocity" << endl;

 //   cout << "COMMUNICATION FINISHED" << endl;

 //   cout << "MADE BIGGER" << endl;
    // 
MPI_Barrier(MPI_COMM_WORLD);
    b.Energy_Calculation();

MPI_Barrier(MPI_COMM_WORLD);

    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();
    hrc::time_point end = hrc::now();
    auto diff = end - start;
    cout <<"\n\nexexution time: "<< chrono::duration <double, milli> (diff).count() << " ms" << endl;
//cout << "getting to the barrier" << endl;   
    //MPI_Barrier(MPI_COMM_WORLD);
 //  cout<< "Passed the BARRIER" << endl;
// b.~Burgers();
// cout<< "Passed the destructor" << endl;
 MPI_Barrier(MPI_COMM_WORLD);
 //cout << "DESTRUCTED" << endl;
  MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
   // cout<<"Finished"<<endl;
    return 0;
}
