#ifndef Burgers_h
#define Burgers_h
#include <cmath>
using namespace std;
class Burgers
{
public:
    Burgers(Model* m);

    // Initial_velocity()
    // I
    // cout<<"ax = "<<ax<<endl;
    void Initial_velocity();
    void Integrate_velocity();
    void Print_velocity();
    void Energy_Calculation();
private:
    Model* m
    double* x, y u, v;
    int rows, cols;

};

Burgers::Burgers(Model* m_) : m(m_)
{
    rows = m->getNy();
    cols  =  m->getNx();
    double arr_size = rows * cols;
    u = new double[arr_size];
    v = new double[arr_size];

    x = new double[cols];
    y = new double[rows];
    double L = m->getLx());
    double dx = m->getDx();
    double dy = m->getDy();
    for(int i = 0; i < cols; ++i) {
        x[i] = i*dx - L / 2.0;
    }
    for(int i = 0; i < rows; ++i) {
        y[i] = i*dy - L / 2.0;
    }
}


void Burgers::Initial_velocity(){
    // Initial Conditions
    double r;
    for  (int i = 0; i < cols; i++){
        for  (int j = 0; j < rows; j++){
            r = sqrt(pow(x[i],2)+pow(y[j],2));
            v[i*rows+j] = u[i*rows+j] = (r<1? 2*pow((1-r),4)*((4*r+1)): 0)
        }
    }
}




void Integrate_velocity()
{
    
}

void Print_velocity()
{
}

void Energy_Calculation()
{
}

#endif /* Burgers_h */