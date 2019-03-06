#ifndef Burgers_h
#define Burgers_h
#include <cmath>
class Burgers
{
public:
    Burgers(Model* m);

    // Initial_velocity()
    // I
    // cout<<"ax = "<<ax<<endl;
    void Initial_velocity();
    double Integrate_velocity();
    double Print_velocity();
    double Energy_Calculation();

private:
    //Numerical
    double x0;
    double y0;
    double Lx;
    double Ly;
    double L;
    double T = 1;
    int Nx = 20;
    int Ny = 20;
    int Nt = 30;
    double dx = L / (Nx - 1);
    double dy = L / (Ny - 1);
    double dt = T / Nt;

    // Physics
    double ax;
    double ay;
    double b;
    double c;
    // arrays needed for the functions to be called 
    double uArray[Nx][Ny];
    double unArray[Nx][Ny];
    double vArray[Nx][Ny];
    double vnArray[Nx][Ny];
    double rhsxArray[Nx][Ny];
    double rhsyArray[Nx][Ny];
    double xArray[Nx][Ny];
    double x[Nx];
    double y[Ny];
    double u[Nx][Ny];
    double v[Nx][Ny];
    double un[Nx][Ny];
    double vn[Nx][Ny];



};

Burgers::Burgers(Model* m)
{
    // taking in all the inputs
    ax = m->GetAx();
    ay = m->GetAy();
    b = m->GetB();
    c = m->GetC();
    x0 = m->GetX0();
    y0 = m->GetY0();
    Lx = m->GetLx();
    Ly = m->GetLy();
    L = m->GetL();
    T = m->GetT();
    Nx = m->GetNx();
    Ny = m->GetNy();
    Nt = m->GetNt();
    dx = m->GetDx();
    dy = m->GetDy();
    dt = m->GetDt();

    
    
    // variables needed for the array multiplactions and stuff


}
void Burgers::Initial_velocity(){
    double uArray[Nx][Ny];
    double unArray[Nx][Ny];
    double vArray[Nx][Ny];
    double vnArray[Nx][Ny];
    double rhsxArray[Nx][Ny];
    double rhsyArray[Nx][Ny];
    double xArray[Nx][Ny];
    double x[Nx];
    double y[Ny];
    double u[Nx][Ny];
    double v[Nx][Ny];
    double un[Nx][Ny];
    double vn[Nx][Ny];


x[0]=-L/2;
y[0]=-L/2;

for (int temp=1; x[Nx]>0; temp++){
    x[temp]=x[temp-1]+dx;
}

for (int temp=1;  y[Ny]>0;temp++){
    y[temp]=y[temp-1]+dy;
}

// Initial Conditions

for  (int i = 0; i < Nx; i++){
    for  (int j = 0; j < Ny; j++){
        double r = pow(pow(x[i],2)+pow(y[j],2),2);
        if (r<=1){
            u[i][j]= 2 * pow(1-r,4) * (4 *r +1);
            v[i][j]= u[i][j];
        }
            else {
                u[i][j]=0;
                v[i][j]=0;
        }
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