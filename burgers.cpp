#include <chrono>

#include "Model.h"
#include "Burgers.h"
#include <iostream>
#include <stdio.h>



int main(int argc, char* argv[]) {

    Model m(argc, argv);
    Burgers b(&m);

    // Call code to initialise the problem here
//    int Burgers::Initial_velocity(Model* m){
//        return ax;
//    }

 b.Initial_velocity();

//cout << u<< v << endl;
    
    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Call code to perform time integration here

    hrc::time_point end = hrc::now();

    // Calculate final energy and write output

    return 0;
}