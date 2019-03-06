#ifndef CLASS_MODEL
#define CLASS_MODEL

//#include "mpi.h"
#include <iostream>
using namespace std;

class Model {
	public:
            Model(int argc, char* argv[])
            {

                for(int i = 0; i < argc; i++) {
                    if(string(argv[i]) == "-h") {
                        DispHelp();
                        help = true;
                    }
                    if(string(argv[i]) == "-v") {
                        PrintParameters();
                        verbose = true;
                    }
                }

                ParseParameters(argc, argv);
            }

            //~Model();

        // Getters
        bool IsVerbose() const{ return verbose; }
        bool   IsHelp()    const { return help; }
        double GetX0()     const { return x0; }
        double GetY0()     const { return y0; }
        double GetL()       const {return L; }
        double GetLx()     const { return Lx; }
        double GetLy()     const { return Ly; }
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

        // Add any other getters here...

private:
    
    // Insert parameters here..
    void ParseParameters(int argc, char* argv[])
    {
        if(argc != 6 && argc != 5) {
            cout << " Please refer to the *HELP* section below" << endl;
            DispHelp();
        }

        ax = stod(argv[1]);
        ay = stod(argv[2]);
        b = stod(argv[3]);
        c = stod(argv[4]);
        L = 10;
        Lx = 10;
        Ly = 10;
        T = 0;
    }

    bool verbose;
    bool help;

    // Prints the parameter list if the user calls for it
    void PrintParameters()
    {
        cout << "Physics constants: ax, ay, b, c" << endl;
        cout << "Numerical constants: x0, y0, Lx, Ly, T, Nx, Ny, Nt, dx, dy, dt" << endl;
    }

    // Displays a help menu if the wrong number of inputs are entered or the user calls for aid
    void DispHelp()
    {
        cout << "HELP" << endl;
        cout << " if you would like help, please enter '-h' as the final input" << endl;
        cout << " if you would to know the parameter list, please enter '-v' as the final input" << endl;
        cout << " this program does blah blah" << endl;
        cout << " the following parameters are required" << endl;
        cout << " ax = constant" << endl;
        cout << " ay = constant" << endl;
        cout << " b = constant" << endl;
        cout << " c = constant" << endl;
        }

//        void ValidateParameters(argv);
//        for (int i = 0; i<4;i++ )
//        if(isalpha(argv[i]) == true) {
//            bool valid = false;
//            return valid
//        }

        bool IsValid(); // returns true or false. true if all parameters are valid

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

        // Add any additional parameters here...

};

#endif  