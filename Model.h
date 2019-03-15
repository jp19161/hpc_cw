#ifndef CLASS_MODEL
#define CLASS_MODEL

#include <regex>
#include <iostream>
#include <fstream>
//#include <mpi.h>
using namespace std;

class Model {
	public:
            Model(int argc, char* argv[]);
            ~Model(){};
        // Getters
        bool IsVerbose() const{ return verbose; }
        bool   IsHelp()    const { return help; }
        double GetX0()     const { return x0; }
        double GetY0()     const { return y0; }
        double GetL()       const {return L; }
        //double GetLx()     const { return Lx; }
        //double GetLy()     const { return Ly; }
        double GetT()      const { return T; }
        int    GetNx()     const { return Nx; }
        int    GetNy()     const { return Ny; }
        int    GetNt()     const { return Nt; }
        double GetDx()     const { return dx; }
        double GetDy()     const { return dy; }
        double GetDt()     const { return dt; }
        double GetAx()     const { return ax; }
        double GetAy()     const { return ay; }
        double GetB()      const { return b; }
        double GetC()      const { return c; }
        int GetPx()      const { return px; }
        int GetPy()      const { return py; }

        

        // Add any other getters here...

private:
    
    // Insert parameters here..
    void ParseParameters(int argc, char* argv[]);

    bool verbose;
    bool help;

    // Prints the parameter list if the user calls for it
    void PrintParameters();
    // Calculates the rest of the 'D' parameters
    void FillInput();

    // Displays a help menu if the wrong number of inputs are entered or the user calls for aid
    void DispHelp();

    void IsValid(int argc, char* argv[]);

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
        int px;
        int py;

        // Physics
        double ax;
        double ay;
        double b;
        double c;

        // Add any additional parameters here...
};

#endif  