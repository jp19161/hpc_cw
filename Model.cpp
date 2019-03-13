#include "Model.h"
using namespace  std;




void Model::ParseParameters(int argc, char* argv[])
{
    if(argc !=10) {
        cout << " Please refer to the *HELP* section below" << endl;
        DispHelp();
    }
    ax = stod(argv[1]);
    ay = stod(argv[2]);
    b = stod(argv[3]);
    c = stod(argv[4]);
    L = stod(argv[5]);
    Nx = stoi(argv[6], NULL, 10);
    Ny = stoi(argv[7], NULL,10);
    Nt = stoi(argv[8],NULL, 10);
    T = stod(argv[9]);
}

void Model::PrintParameters()
{
    cout << "Physics constants: ax, ay, b, c" << endl;
    cout << "Numerical constants: x0, y0, L, T, Nx, Ny, Nt, dx, dy, dt" << endl;
    // print properly
        }

        void Model::DispHelp()
    {
        cout << "HELP" << endl;
        cout << " if you would like help, please enter '-h' as the final input" << endl;
        cout << " if you would to know the parameter list, please enter '-v' as the final input" << endl;
        cout << " this program does blah blah" << endl;
        cout << " the following parameters are required:" << endl;
        cout << " ax = constant" << endl;
        cout << " ay = constant" << endl;
        cout << " b = constant" << endl;
        cout << " c = constant" << endl;
        cout << " L = length" << endl;
        cout << " Nx = Number of x points" << endl;
        cout << " Ny = Number of y points" << endl;
        cout << " Nt = Number of time points" << endl;
        cout << " T = time" << endl;
        }

        void Model::IsValid(int argc, char* argv[])
        {
            regex integer("(\\+|-)?[[:digit:]]+");
            // As long as the input is correct ask for another number
            for(int i = 0; i < 3; i++) {
                if(!regex_match(argv[i + 6], integer)) {
                    cout << "Nx, Ny, Nt must be integers" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }

        void Model::FillInput()
        {
            dx = L / (Nx - 1);
            dy = L / (Ny - 1);
            dx = T / Nt;
        }

        Model::Model(int argc, char* argv[])
        {
            // Reading the parameters
            ParseParameters(argc, argv);
            FillInput();
            // ValidateParameters()
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

            
            
        }