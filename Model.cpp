#include "Model.h"
using namespace  std;

void Model::ParseParameters(int argc, char* argv[])
{
    // Get rank of process 

    
    if(argc != 12) {
        cout << " Please refer to the *HELP* section below" << endl;
        DispHelp();
    }
    ax = stod(argv[1]);
    ay = stod(argv[2]);
    b = stod(argv[3]);
    c = stod(argv[4]);
    L = stod(argv[5]);
    Nx = stoi(argv[6], NULL, 10);
    Ny = stoi(argv[7], NULL, 10);
    Nt = stoi(argv[8], NULL, 10);
    T = stod(argv[9]);
    px = stoi(argv[10], NULL, 10);
    py = stoi(argv[11], NULL, 10);
}

void Model::PrintParameters()
{
    cout << "ax:" << ax << endl;
    cout << "ay:" << ay << endl;
    cout << "b:" << b << endl;
    cout << "c:" << c << endl;
    cout << "L:" << L << endl;
    cout << "Nx:" << Nx << endl;
    cout << "Ny:" << Ny << endl;
    cout << "Nt:" << Nt << endl;
    cout << "T:" << T << endl;
    cout << "px:" << px << endl;
    cout << "py:" << py << endl;

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
   // for(int i = 0; i < 3; i++) {
        if(!regex_match(argv[6], integer)) {
            cout << "Nx must be an integer" << endl;
            exit(EXIT_FAILURE); 
        }
                if(!regex_match(argv[7], integer)) {
            cout << "Ny must be an integer" << endl;
            exit(EXIT_FAILURE); 
        }
                if(!regex_match(argv[8], integer)) {
            cout << "Nt must be an integer" << endl;
            exit(EXIT_FAILURE); 
        }
                if(!regex_match(argv[10], integer)) {
            cout << "px must be an integer" << endl;
            exit(EXIT_FAILURE); 
        }
                if(!regex_match(argv[11], integer)) {
            cout << "py must be an integer" << endl;
            exit(EXIT_FAILURE); 
        }
        
   // }
}

void Model::FillInput()
{
    dx = L / (Nx - 1);
    dy = L / (Ny - 1);
    dt = T / Nt;
}

Model::Model(int argc, char* argv[])
{
    // Reading the parameters
    //cout<<"Creating mode"<<endl;
    ParseParameters(argc, argv);
    //cout<<"Parsed parameters"<<endl;
    FillInput();
    //cout<<"Filled"<<endl;
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