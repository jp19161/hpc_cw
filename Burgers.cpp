#include "Burgers.h"
#include "Model.h"
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <stdio.h>
#include <iomanip>
using namespace std;

Burgers::Burgers(Model& m)
{
    // need to calculate coefficients for integration

    // parameters needed for this sections funcitons
    ax = m.GetAx();
    ay = m.GetAy();
    b = m.GetB();
    c = m.GetC();
    L = m.GetL();
    Nx = m.GetNx();
    Ny = m.GetNy();
    Nt = m.GetNt();
    T = m.GetT();
    dx = m.GetDx();
    dy = m.GetDy();
    dt = m.GetDt();
    px = m.GetPx();
    py = m.GetPy();

    // MPI stuffs

    retval_rank = MPI_Comm_rank(MPI_COMM_WORLD, &myrank); // zero-based
    retval_size = MPI_Comm_size(MPI_COMM_WORLD, &mysize);
    // cout << "Current rank is:" << myrank << endl;
    // cout << " Current size number is" << mysize << endl;
    if(retval_rank == MPI_ERR_COMM || retval_size == MPI_ERR_COMM) {
        std::cout << "Invalid communicator" << std::endl;
    }
}

Burgers::~Burgers()
{
    // deleteing all the paointed to variables
    //        delete[]  u;
    //        delete[]  v;
    //        delete[]  u_new;
    //        delete[]  v_new;
    //        delete[]  x;
    //        delete[]  y;
    //        delete[]  width_process_u;
    //        delete[]  height_process_u;
    //        delete[] width_process_v;
    //        delete[] height_process_v;
    //        delete[]  outer_left_u;
    //        delete[] outer_right_u;
    //        delete[] inner_left_u;
    //        delete[] inner_right_u;
    //        delete[] inner_down_u;
    //        delete[] inner_up_u;
    //        delete[] outer_down_u;
    //        delete[] outer_up_u;
    //        delete[] local_u;
    //        delete[]  local_v;
    //        delete[]  local_x;
    //        delete[]  local_y;
    //        delete[]  local_u_new;
    //        delete[]  local_v_new;
    //        delete[] BigU;
    //        delete[] BigV;
    //        delete[] bigu;
    //        delete[] bigv;
    //       delete[] medium_u;
}

// Initialise velocity
void Burgers::Initial_velocity()
{

    // dividing up sub grids
    mpirows = myrank / px;
    mpicols = myrank % px;

    // width of local Nx
    // cout << "Creating no rows" << endl;
    width_process_u = new int[px];
    for(int i = 0; i <= px; i++)
        width_process_u[i] = Nx / px;
    // cout << "Got here" << endl;
    if(Nx % px != 0) {
        for(int i = 0; i < (Nx % px); i++)
            width_process_u[i] = width_process_u[i] + 1;
    }
    // Height of local Ny
    height_process_u = new int[py];
    for(int i = 0; i <= py; i++)
        height_process_u[i] = Ny / py;
    // width_process_u = Ny/py;
    if(Ny % py != 0) {
        for(int i = 0; i < (Ny % py); i++)
            height_process_u[i] = height_process_u[i] + 1;
    }

    local_nx = width_process_u[mpicols];
    local_ny = height_process_u[mpirows];
    // cout<<"Created rows"<<endl;

    // finding global cordinates
    for(int i = 0; i < mpicols; i++) {
        my_pos_x += width_process_u[i];
    }
    for(int i = 0; i < mpirows; i++) {
        my_pos_y += height_process_u[i];
    }

    // checking everything has been definined correctly
//    cout << "rank: " << myrank << endl;
//
//    cout << "local nx: " << local_nx << endl;
//
//    cout << "local ny: " << local_ny << endl;
//
//    cout << "my x position: " << my_pos_x << endl;
//
//    cout << "my y position: " << my_pos_y << endl;

    // Actually doing the MPI stuff hah
    double r, rb;
    local_u = new double[local_nx * local_ny];
    local_v = new double[local_nx * local_ny];
    local_x = new double[local_nx];
    local_y = new double[local_ny];
    /*Left_B = new double
    Right_B
    Upper_B
    Lower_B*/

    //    double r, rb;
    //    u = new double[Nx * Ny];
    //    v = new double[Nx * Ny];
    //    x = new double[Nx];
    //    y = new double[Ny];

    for(int i = 0; i < local_nx; i++) {
        local_x[i] = -L / 2 + (my_pos_x + i) * dx;
    }
    for(int i = 0; i < local_ny; i++) {
        local_y[i] = -L / 2 + (my_pos_y + i) * dy;
    }
    // cout<<local_x[0]<<" "<<local_y[0]<<endl;

    for(int col = 0; col < local_nx; ++col) {
        for(int row = 0; row < local_ny; ++row) {
            r = sqrt(pow(local_x[col], 2) + pow(local_y[row], 2));
            if(r <= 1) {
                rb = 2 * pow(1 - r, 4) * (4 * r + 1);
               // cout << "my rank:" << myrank << " global i " << (col + my_pos_x) << " global j " << (my_pos_y + row)
                     //<< endl;
            } else {
                rb = 0.0;
            }
            local_u[row * local_nx + col] = rb;
            local_v[row * local_nx + col] = rb;
        }
    }

    inner_left_u = new double[local_ny];
    inner_right_u = new double[local_ny];
    outer_left_u = new double[local_ny];
    outer_right_u = new double[local_ny];
    inner_up_u = new double[local_nx];
    inner_down_u = new double[local_nx];
    outer_up_u = new double[local_nx];
    outer_down_u = new double[local_nx];

    inner_left_v = new double[local_ny];
    inner_right_v = new double[local_ny];
    outer_left_v = new double[local_ny];
    outer_right_v = new double[local_ny];
    inner_up_v = new double[local_nx];
    inner_down_v = new double[local_nx];
    outer_up_v = new double[local_nx];
    outer_down_v = new double[local_nx];

    for(int row = 0; row < local_ny; row++) {
        inner_left_u[row] = local_u[row * local_nx]; // first element of the top is local u [local_ny*local_nx-local_nx]
        inner_left_v[row] = local_v[row * local_nx];
        inner_right_u[row] = local_u[row * local_nx + local_nx - 1];
        inner_right_v[row] = local_v[row * local_nx + local_nx - 1];
        // do v
    }
    // different loop for top bottom
    for(int col = 0; col < local_nx; col++) {
        inner_up_u[col] = local_u[col +local_nx*(local_ny-1)]; // first element of the top is local u [local_ny*local_nx-local_nx]
        inner_up_v[col] = local_v[col +local_nx*(local_ny-1)];
        inner_down_u[col] = local_u[col];
        inner_down_v[col] = local_v[col];
        // do v
    }
//    if (myrank==3)
//        for(int col = 0; col < local_nx; col++)
//            cout<<inner_down_u[col]<<" ";
//        cout<<endl;

    MPI_Barrier(MPI_COMM_WORLD);
    //    if (myrank==2)
    //    for(int row = local_ny-1; row >=0; --row) {
    //        /*for(int col = 0; col < local_nx; ++col) {
    //
    //            //cout <<local_u[row*local_nx+col]<<" ";
    //        }*/
    //       // cout <<endl;
    //    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*    for(int i = 0; i < local_nx; i++) {
            u[i * local_ny + 0] = 0.0;
            u[i * local_ny + local_ny - 1] = 0.0;
            v[i * local_ny + 0] = 0.0;
            v[i * local_ny + local_ny - 1] = 0.0;
        }
        for(int i = 0; i < local_nx; i++) {
            u[0 * local_ny + i] = 0.0;
            u[(local_nx - 1) * local_ny + i] = 0.0;
            v[0 * local_ny + i] = 0.0;
            v[(local_nx - 1) * local_ny + i] = 0.0;
        }*/
    // cout << "sup " << endl;
}

/*void Burgers::AssembleGrid()
{
double full_grid_x;
double full_grid_y;



}*/

/*Patching U field by first patching the fields onto the proceses on the leftmost column
 * and then by patchig the wide fields onto process 0 at the bottom
 */



// MPI references
// int MPI_Recv(void *buf, int cnt, MPI_Datatype type, int src,int tag, MPI_Comm comm, MPI_Status *stat);
// int MPI_Send(void *buf, int cnt, MPI_Datatype type, int dest, int tag, MPI_Comm comm);
void Burgers::Communication()
{
    // sending inner right to outer left of prcess to the right of me
    // cout << " I ENTERED COMMUNCATION()" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    if(mpicols != 0) {
        MPI_Recv(outer_left_u, local_ny, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(outer_left_v, local_ny, MPI_DOUBLE, myrank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        cout << " my rank is: " << myrank << " and i am receiving from: " << myrank - 1 << endl;
//        if (myrank==1)
//        for(int i = local_ny-1; i>=0; i--) {
//            cout << outer_left_u[i] << endl;
//        }
//        cout << endl; 
//        cout << endl;
        // Sending left straight afterwards
    }
    if(mpicols != (px - 1)) {
        MPI_Send(inner_right_u, local_ny, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD);
        MPI_Send(inner_right_v, local_ny, MPI_DOUBLE, myrank + 1, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    // sending inner left to outer right of processs to the left of me
    if(mpicols != (px - 1)) {
        MPI_Recv(outer_right_u, local_ny, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(outer_right_v, local_ny, MPI_DOUBLE, myrank + 1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       cout << " my rank is: " << myrank << " and i am receiving from: " << myrank + 1 << endl;
//       for(int i = 0; i < local_ny; i++) {
//           cout << outer_right_u[i] << endl;
//       } 
//       cout << endl;
//       cout << endl;
    }
    //cout << " I FINIHSED THE RECIEVINBG BIT" << endl;
    if(mpicols != 0) {
        MPI_Send(inner_left_u, local_ny, MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD);
        MPI_Send(inner_left_v, local_ny, MPI_DOUBLE, myrank - 1, 1, MPI_COMM_WORLD);
    }

    // seding inner upper to outer down of process above me
    if(mpirows != 0) {
        MPI_Recv(outer_down_u, local_nx, MPI_DOUBLE, myrank - px, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(outer_down_v, local_nx, MPI_DOUBLE, myrank - px, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//        cout << " my rank is: " << myrank << " and i am receiving from: " << myrank - px << endl;
//        for(int i = 0; i < local_nx; i++) {
//            cout << outer_down_u[i] << " ";
//        }
//        cout << endl;
    }
   // MPI_Barrier(MPI_COMM_WORLD);
    if(mpirows != (py - 1)) {
        MPI_Ssend(inner_up_u, local_nx, MPI_DOUBLE, myrank + px, 2, MPI_COMM_WORLD);
        MPI_Ssend(inner_up_v, local_nx, MPI_DOUBLE, myrank + px, 2, MPI_COMM_WORLD);
        //cout << " my rank is: " << myrank << " and i am sendin to: " << myrank + px << endl;
//        for(int i = 0; i < local_nx; i++) {
//            cout << inner_up_u[i] << " ";
//        }
//        cout << endl;
    }

         // sending inner down u to outer upper down u of process below me
        if(mpirows != (py - 1)) {
            MPI_Recv(outer_up_u, local_nx, MPI_DOUBLE, myrank + px, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(outer_up_v, local_nx, MPI_DOUBLE, myrank + px, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//            cout << " my rank is: " << myrank << " and i am receiving from: " << myrank + px << endl;
//            for(int i = 0; i < local_ny; i++) {
//                cout << outer_up_u[i] << " ";
 //           }
        }

        if(mpirows != 0) {
          //  cout << myrank << " RECIEVING from " << (myrank - px) << endl;
            MPI_Ssend(inner_down_u, local_nx, MPI_DOUBLE, myrank - px, 2, MPI_COMM_WORLD);
            MPI_Ssend(inner_down_v, local_nx, MPI_DOUBLE, myrank - px, 2, MPI_COMM_WORLD);

            }
            
            MPI_Barrier(MPI_COMM_WORLD);
//             if (myrank==2)
//        for(int i = 0;i<local_nx;i++) {
//            cout << outer_down_u[i] << "  ";
//       }
//        cout<<endl;
            MPI_Barrier(MPI_COMM_WORLD);

            // checking if it all works
            // cout << "hereherherheherhehre" << endl;
            // for (int i=1; i<local_ny ;i++)
            //{
            //    cout<< "the outer left U is:" << endl;
            //     cout<< outer_left_u[i] << endl;
            //}
            // cout << "LEAVING COMMS" << endl;
        }
    


void Burgers::MakeBigger()
{

    MPI_Barrier(MPI_COMM_WORLD);
   // cout << myrank << "is finished 1" << endl;
    for(int j = 1; j <= local_ny; j++)
        for(int i = 1; i <= local_nx; i++) {
            medium_u[j * (local_nx + 2) + i] = local_u[(j - 1) * local_nx + (i - 1)];
            medium_v[j * (local_nx + 2) + i] = local_v[(j - 1) * local_nx + (i - 1)];
        }

   // cout << myrank << "is finished 2" << endl;
    for(int i = 1; i <= local_ny; i++) {
        medium_u[i * (local_nx + 2)] = outer_left_u[i - 1];
        medium_u[i * (local_nx + 2) + local_nx + 1] = outer_right_u[i - 1];
        medium_v[i * (local_nx + 2)] = outer_left_v[i - 1];
        medium_v[i * (local_nx + 2) + local_nx + 1] = outer_right_v[i - 1];
    }
   // cout << myrank << "is finished 3" << endl;
    for(int i = 1; i <= local_nx; i++) {
        medium_u[i] = outer_down_u[i - 1];
        medium_u[local_ny * (local_nx + 2) + local_nx +2+ i] = outer_up_u[i - 1];
        medium_v[i] = outer_down_v[i - 1];
        medium_v[local_ny * (local_nx + 2) + local_nx +2+ i] = outer_up_v[i - 1];
    }
    //cout << myrank << "is finished 4" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
//    if (myrank==1)
//        for (int j=local_ny+1;j>=0;j--)
//        {
//            for (int i=0;i<medium_nx;i++)
//                cout<<" "<<medium_u[j * medium_nx + i];
//                cout<<endl;
//        }
    MPI_Barrier(MPI_COMM_WORLD);
}

void Burgers::Integrate_velocity()
{
    medium_u = new double[(local_nx + 2) * (local_ny + 2)];
    medium_v = new double[(local_nx + 2) * (local_ny + 2)];
 medium_nx = local_nx + 2;
 for (int i=0; i<local_ny; i++){
    outer_left_u[i] =0;
    outer_right_u[i]=0;
    outer_left_v[i]=0;
    outer_right_v[i]=0;
 }
  for (int i=0; i<local_nx; i++){
    outer_up_v[i]=0;
    outer_down_v[i]=0;
    outer_up_u[i]=0;
    outer_down_u[i]=0;
    }
    
    //    //ORDER
    //    // 1. set outer boundaries that arent involved in integration
    //    // 2. communication boundaries between processes
    //    // 3. Integrate for 1 time steps
    //    // 4. communcate
    //    // 5. repeat for all time steps.
    
            const double c1 = -2.0 * c * (1 / dx / dx + 1 / dy / dy) - ax / dx - ay / dy + 1 / dt;
        const double c2 = c / dx / dx;
        const double c3 = c / dy / dy;
        const double c4 = c / dx / dx + ax / dx;
        const double c5 = c / dy / dy + ay / dy;
    double t=0;
//        cout << "-----------------------------------------------------------------------------------------------------------" << endl;
//    cout << "dt: " << dt << endl;
//    cout << "Nt: " << Nt << endl;
//        cout << "c1/dt: " << c1*dt  << endl;
//    cout << "c2/dt: " << c2 * dt<< endl;
//        cout << "c3/dt: " << c3 *dt << endl;
//    cout << "c4/dt: " << c4 *dt  << endl;
//        cout << "c5/dt: " << c5 *dt << endl;
//    cout << "-----------------------------------------------------------------------------------------------------------" << endl;
    

    while(t <(T))   {
        t+=dt;
        for(int i = 0; i < local_nx+2; i++) {
            for(int j = 0; j < local_ny+2; j++) {
                medium_u[j * medium_nx + i] = 0;
                medium_v[j * medium_nx + i] = 0;
            }
        }
//        if (myrank==1)
//        for (int j=local_ny+1;j>=0;j--)
//        {
//            for (int i=0;i<medium_nx;i++)
//                cout<<" "<<medium_u[j * medium_nx + i];
//                cout<<endl;
//        }
        MPI_Barrier(MPI_COMM_WORLD);
        Communication();
        MPI_Barrier(MPI_COMM_WORLD);
        MakeBigger();
      //  cout<<"got here "<<myrank<<endl;
//        
//        cout<<"got here2 "<<myrank<<endl;
//        if (myrank==1)
//        for (int j=local_ny+1;j>=0;j--)
//        {
//            for (int i=0;i<medium_nx;i++)
//                cout<<" "<<medium_u[j * medium_nx + i];
//                cout<<endl;
//        }
//        cout<<"Feeling fine"<<endl;
         MPI_Barrier(MPI_COMM_WORLD);
//        if(myrank == 3) {
//        for(int j =  local_ny + 1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx+2; i++) {
//                cout << medium_u[j * (local_nx+2) + i]  << " ";
//                
//            }
//            cout << endl;
//        }
//    }

//        
//        
        for(int i = 1; i < (local_nx + 1); i++) {
            for(int j = 1; j < (local_ny + 1); j++) {
                local_u[(j - 1) * local_nx + i - 1] =
                    dt * (c1 * medium_u[j * medium_nx + i] -
                             (b / dy) * medium_v[j * medium_nx + i] * medium_u[j * medium_nx + i] -
                             (b / dx) * medium_u[j * medium_nx + i] * medium_u[j * medium_nx + i] +
                             c2 * medium_u[j * medium_nx + (i + 1)] + c3 * medium_u[(j + 1) * medium_nx + i] +
                             c4 * medium_u[j * medium_nx + (i - 1)] +
                             (b / dx) * medium_u[j * medium_nx + i] * medium_u[j * medium_nx + (i - 1)] +
                             c5 * medium_u[(j - 1) * medium_nx + i] +
                             (b / dy) * medium_v[j * medium_nx + i] * medium_u[(j - 1) * medium_nx + i]);
                             
//
                local_v[(j - 1) * local_nx + i - 1] =
                    dt * (c1 * medium_v[j * medium_nx + i] -
                             (b / dy) * medium_u[j * medium_nx + i] * medium_v[j * medium_nx + i] -
                             (b / dx) * medium_v[j * medium_nx + i] * medium_v[j * medium_nx + i] +
                             c2 * medium_v[j * medium_nx + (i + 1)] + c3 * medium_v[(j + 1) * medium_nx + i] +
                             c4 * medium_v[j * medium_nx + (i - 1)] +
                             (b / dx) * medium_v[j * medium_nx + i] * medium_v[j * medium_nx + (i - 1)] +
                             c5 * medium_v[(j - 1) * medium_nx + i] +
                             (b / dy) * medium_u[j * medium_nx + i] * medium_v[(j - 1) * medium_nx + i]);
            }
        }
//        
        for(int row = 0; row < local_ny; row++) {
            inner_left_u[row] =
            local_u[row * local_nx]; // first element of the top is local u [local_ny*local_nx-local_nx]
            inner_left_v[row] = local_v[row * local_nx];
            inner_right_u[row] = local_u[row * local_nx + local_nx - 1];
            inner_right_v[row] = local_v[row * local_nx + local_nx - 1];
        }
//        // different loop for top bottom
        for(int col = 0; col < local_nx; col++) {
            inner_up_u[col] =
            local_u[col +local_nx*(local_ny-1)]; // first element of the top is local u [local_ny*local_nx-local_nx]
            inner_up_v[col] = local_v[col +local_nx*(local_ny-1)];
            inner_down_u[col] = local_u[col];
            inner_down_v[col] = local_v[col];
        }
//        MPI_Barrier(MPI_COMM_WORLD);
//cout<<"Feeling fin here as wele"<<endl;
//         MPI_Barrier(MPI_COMM_WORLD);
////        MPI_Barrier(MPI_COMM_WORLD);
////        if (myrank==0){
        cout <<"time step is:   "<< t <<endl;
      t+=dt;
////        cout<<"Dt after"<<dt<<endl;
////        }
//
//cout << " the  medium u velocity" << endl;
//if (myrank == 3){
//    cout << "my rank is :" << myrank << endl;
//for(int j =  local_ny + 1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx+2; i++) {
//                cout << medium_u[j * (local_nx+2) + i]  << " ";
//                
//            }
//            cout << endl;
// }
        }






//
//MPI_Barrier(MPI_COMM_WORLD);
//
//cout << " the  local u velocity for rank 3" << endl;
//if (myrank == 3){
//    cout << "my rank is :" << myrank << endl;
//for(int j =  local_ny -1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx; i++) {
//                cout << local_u[j * (local_nx) + i]  << " ";
//                
//            }
//            cout << endl;
//        }
//}
//MPI_Barrier(MPI_COMM_WORLD);
//cout << " the  local u velocity for rank 0" << endl;
//if (myrank == 0){
//    cout << "my rank is :" << myrank << endl;
//for(int j =  local_ny -1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx; i++) {
//                cout << local_u[j * (local_nx) + i]  << " ";
//                
//            }
//            cout << endl;
//        }
//}
}
//}
//MPI_Barrier(MPI_COMM_WORLD);
//if (myrank == 1){
//    cout << "my rank is :" << myrank << endl;
//for(int j =  local_ny + 1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx+2; i++) {
//                cout << local_u[j * (local_nx+2) + i]  << " ";
//                
//            }
//            cout << endl;
//        }
//}
//MPI_Barrier(MPI_COMM_WORLD);
//
//if (myrank == 2){
//    cout << "my rank is :" << myrank << endl;
//for(int j =  local_ny + 1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx+2; i++) {
//                cout << local_u[j * (local_nx+2) + i]  << " ";
//                
//            }
//            cout << endl;
//        }
//}
//MPI_Barrier(MPI_COMM_WORLD);
//
//if (myrank == 3){
//    cout << "my rank is :" << myrank << endl;
//for(int j =  local_ny + 1; j>=0 ; j--) {
//            for(int i = 0; i < local_nx+2; i++) {
//                cout << local_u[j * (local_nx+2) + i]  << " ";
//                
//            }
//            cout << endl;
//        }
//}



    
//     for(int i = 1; i < local_nx - 1; i++) {
//            for(int j = 1; j < local_ny - 1; j++) {
//                local_u_next = c1 * u[(i - 1) * local_ny + j] + c2 * u[(i + 1) * local_ny + j] + c3 * u[i * local_ny +
//   j] + c4 * u[i * local_ny + j - 1] + c5 * u[i * local_ny + j + 1] +b / dx * u[i * local_ny + j] * (-u[i * local_ny +
//   j] + u[(i - 1) * local_ny + j]) +  b / dy * v[i * local_ny + j] * (u[i * local_ny + j - 1] - u[i * local_ny + j]);
//                local_v_next = c1 * v[(i - 1) * local_ny + j] + c2 * v[(i + 1) * local_ny + j] + c3 * v[i * local_ny +
//   j] +  c4 * v[i * local_ny + j - 1] + c5 * v[i * local_ny + j + 1] + b / dx * v[i * local_ny + j] * (-v[i * local_ny +
//   j] + v[(i - 1) * local_ny + j]) +  b / dy * u[i * local_ny + j] * (v[i * local_ny + j - 1] - v[i * local_ny + j]);
//                //cout << "sup 3" << endl;
//                local_u_new[i * local_ny + j] = dt * local_u_next;
//                local_v_new[i * local_ny + j] = dt * local_v_next;
//            }
//        }*/
// Reassigning values to the velocity arrays
/*        for(int in = 1; in < local_nx - 1; in++) {
            for(int jn = 1; jn < local_ny - 1; jn++) {
                local_u[in * local_ny + jn] = local_u_new[in * local_ny + jn];
                local_v[in * local_ny + jn] = local_v_new[in * local_ny + jn];
            }
        }*/
    //    if(myrank == 1)
    //        for(int j = (local_ny + 1); j >= 0; j--) {
    //            for(int i = 0; i < (local_nx + 2); i++)
    //                cout << setw(8) << fixed << medium_u[j * (local_nx + 2) + i];
    //            cout << endl;


//                local_u_next = c1 * local_u[(i - 1) * local_ny + j] + c2 * local_u[(i + 1) * local_ny + j] + c3 *
//                local_u[i * local_ny + j] + c4 * local_u[i * local_ny + j - 1] + c5 * local_u[i * local_ny + j + 1] +b
//                / dx * local_u[i * local_ny + j] * (-local_u[i * local_ny + j] + local_u[(i - 1) * local_ny + j]) +  b
//                / dy * local_u[i * local_ny + j] * (local_u[i * local_ny + j - 1] - local_u[i * local_ny + j]);
//                local_v_next = c1 * local_v[(i - 1) * local_ny + j] + c2 * local_v[(i + 1) * local_ny + j] + c3 *
//                local_v[i * local_ny + j] +  c4 * local_v[i * local_ny + j - 1] + c5 * local_v[i * local_ny + j + 1] +
//                b / dx * local_v[i * local_ny + j] * (-local_v[i * local_ny + j] + local_v[(i - 1) * local_ny + j]) +
//                b / dy * local_v[i * local_ny + j] * (local_v[i * local_ny + j - 1] - local_v[i * local_ny + j]);
//                //cout << "sup 3" << endl;
//                local_u_new[i * local_ny + j] = dt * local_u_next;
//                local_v_new[i * local_ny + j] = dt * local_v_next;

// auto* u_next = new double[Nx * Ny];
// auto* v_next = 0new double[Nx * Ny];

// cout << "sup2 " << endl;

// dont forget about layer of padding all the way the outside of the big sqaure

/*//bottom left corner
    if(mpicols==0 && mpirows==0){
    }

// top right corner
     if(mpicols==px && mpirows==py){
    }

// top left corner
     if(mpicols==0 && mpirows == py){
    }

// top left corner
     if(mpicols==px && mpirows == 0){
    }

// left rows
    if(mpicols ==0 && mpirows !=0 && mpirows !=py){
    }

//right rows
    if(mpicols==px && mpirows !=0 && mpirows !=py){
    }

// bottom rows
    if(mpicols !=0 && mpicols != px && mpirows ==0){
    }

// top rows
    if (mpicols !=0 && mpicols != px && mpirows ==py){
    }*/

// MPI DATA TYPE METHOD
/*   MPI_Datatype row;
   MPI_Type_vector top_row(1, local_nx, local_nx, MPI_DOUBLE, &top_row)

   MPI_Datatype col;
   MPI_Type_vector right_col(local_ny, 1, local_nx, MPI_DOUBLE, &right_col)

   MPI_Send(local_u, 1, top_col, int destination, int tag, MPI_Comm communicator)
   MPI_Recv(void* data, local_nx, top_col, int source, int tag, MPI_Comm communicator, MPI_Status* status)*/

/*        for(int i = 1; i < local_nx - 1; i++) {
            for(int j = 1; j < local_ny - 1; j++) {
                local_u_next = c1 * u[(i - 1) * local_ny + j] + c2 * u[(i + 1) * local_ny + j] + c3 * u[i * local_ny +
   j] + c4 * u[i * local_ny + j - 1] + c5 * u[i * local_ny + j + 1] +b / dx * u[i * local_ny + j] * (-u[i * local_ny +
   j] + u[(i - 1) * local_ny + j]) +  b / dy * v[i * local_ny + j] * (u[i * local_ny + j - 1] - u[i * local_ny + j]);
                local_v_next = c1 * v[(i - 1) * local_ny + j] + c2 * v[(i + 1) * local_ny + j] + c3 * v[i * local_ny +
   j] +  c4 * v[i * local_ny + j - 1] + c5 * v[i * local_ny + j + 1] + b / dx * v[i * local_ny + j] * (-v[i * local_ny +
   j] + v[(i - 1) * local_ny + j]) +  b / dy * u[i * local_ny + j] * (v[i * local_ny + j - 1] - v[i * local_ny + j]);
                //cout << "sup 3" << endl;
                local_u_new[i * local_ny + j] = dt * local_u_next;
                local_v_new[i * local_ny + j] = dt * local_v_next;
            }
        }*/
// Reassigning values to the velocity arrays
/*        for(int in = 1; in < local_nx - 1; in++) {
            for(int jn = 1; jn < local_ny - 1; jn++) {
                local_u[in * local_ny + jn] = local_u_new[in * local_ny + jn];
                local_v[in * local_ny + jn] = local_v_new[in * local_ny + jn];
            }
        }*/

 void Burgers::Energy_Calculation()
{
    
    if (myrank ==0){
        cout << "did i get here" << endl; 
    double energy = 0.0;
    for(int i = 0; i < Nx * Ny; i++) {
        energy += BigU[i]*BigU[i] + BigV[i]*BigV[i];
    }
    energy *= 0.5 * dx * dy;
    cout << "the energy is: " << energy << endl;
    }
}

// Model m(argc, argv);
// Burgers b(&m);
// b.Initial_velocity();

/*    typedef std::chrono::high_resolution_clock hrc;
    typedef std::chrono::milliseconds ms;
    hrc::time_point start = hrc::now();

    // Call code to perform time integration here

    hrc::time_point end = hrc::now();

    // Calculate final energy and write output

    return 0;*/
    void Burgers::PatchUpU()
{
    // Patching in the row
    if(mpicols == 0) {
        cout << "my rank" << myrank << "entered merging" << endl;
        bigu = new double[Nx * local_ny+5];
        int pos = 0;
        //        bigu = new double [sizeu*local_ny];
        int posnow;
         cout << "my rank" << myrank << "Started merging" << endl;
        // from other in the same row
        for(int i = 1; i < px; i++) {
            int recv_size = width_process_u[i];
            pos += width_process_u[i - 1];

            for(int j = 0; j < local_ny; j++) {
                posnow = j * Nx + pos;
                MPI_Recv(&bigu[posnow], recv_size, MPI_DOUBLE, myrank + i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        // from own
        for(int j = 0; j < local_ny; j++) {
            posnow = j * Nx;
            for(int i = 0; i < local_nx; i++)
                bigu[posnow + i] = local_u[j * local_nx + i];
        }
    } else {
        for(int j = 0; j < local_ny; j++)
            MPI_Ssend(&local_u[j * local_nx], local_nx, MPI_DOUBLE, myrank - mpicols, mpicols, MPI_COMM_WORLD);
    }
    cout << "my rank" << myrank << "Done column merging" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    // Patching the wide field onto process 0
    if(mpicols == 0) {
        if(myrank == 0) {
            int posj = 0;

            BigU = new double[Nx * Ny];
            // from other in the 0 column
            for(int J = 1; J < py; J++) {
                posj += height_process_u[J - 1];
                //cout << "Receiving from " << (myrank + J * px) << endl;
                for(int j = 0; j < height_process_u[J]; j++)
                    MPI_Recv(
                        &BigU[(posj + j) * Nx], Nx, MPI_DOUBLE, myrank + J * px, J, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                //cout << "Received from " << (myrank + J * px) << endl;
            }
            // from own widefield
            for(int j = 0; j < local_ny; j++)
                for(int i = 0; i < Nx; i++)
                    BigU[j * Nx + i] = bigu[j * Nx + i];
        } else {
            //cout << "my rank" << myrank << " Sending to 0" << endl;

            for(int j = 0; j < local_ny; j++)
                MPI_Ssend(&bigu[j * Nx], Nx, MPI_DOUBLE, 0, mpirows, MPI_COMM_WORLD);
        }
    }
    //cout << "my rank finshed patching U " << myrank << endl;
}

void Burgers::PatchUpV()
{
    // Patching in the row
    if(mpicols == 0) {
        bigv = new double[Nx * local_ny];
        int pos = 0;
        //        bigu = new double [sizeu*local_ny];
        int posnow;
        // from other in the same row
        for(int i = 1; i < px; i++) {
            int recv_size = width_process_u[i];
            pos += width_process_u[i - 1];

            for(int j = 0; j < local_ny; j++) {
                posnow = j * Nx + pos;
                MPI_Recv(&bigv[posnow], recv_size, MPI_DOUBLE, myrank + i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
        // from own
        for(int j = 0; j < local_ny; j++) {
            posnow = j * Nx;
            for(int i = 0; i < local_nx; i++)
                bigv[posnow + i] = local_v[j * local_nx + i];
        }
    } else {
        for(int j = 0; j < local_ny; j++)
            MPI_Send(&local_v[j * local_nx], local_nx, MPI_DOUBLE, myrank - mpicols, mpicols, MPI_COMM_WORLD);
    }
    // cout<<"Done column merging"<<endl;
    // Patching the wide field onto process 0
    if(mpicols == 0) {
        if(myrank == 0) {
            int posj = 0;

            BigV = new double[Nx * Ny];
            // from other in the 0 column
            for(int J = 1; J < py; J++) {
                posj += height_process_u[J - 1];
                for(int j = 0; j < height_process_u[J]; j++)
                    MPI_Recv(
                        &BigV[(posj + j) * Nx], Nx, MPI_DOUBLE, myrank + J * px, J, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
            // from own widefield
            for(int j = 0; j < local_ny; j++)
                for(int i = 0; i < Nx; i++)
                    BigV[j * Nx + i] = bigv[j * Nx + i];
        } else {
            for(int j = 0; j < local_ny; j++)
                MPI_Send(&bigv[j * Nx], Nx, MPI_DOUBLE, 0, mpirows, MPI_COMM_WORLD);
        }
    }
}

//void Burgers::Energy_Calculation(){
//    double energy = 0.0;
//    
//    if(myrank == 0){
//    for (int = 0; )
//}

void Burgers::Print_velocity()
{
    if(myrank == 0) {
        // Open file
        ofstream f_out;

        f_out.open("velocity.txt"); //,ios::trunc | ios::out);
        if(!f_out.good()) {
            cout << "Error: unable to open output file: sine.txt" << endl;
        } else {
            f_out << "the U velocity field" << endl;

            f_out.precision(6);
            for(int j = Ny - 1; j >= 0; j--) {
                for(int i = 0; i < Nx; i++)
                    f_out << fixed << BigU[j * Nx + i] << " ";
                f_out << endl;
            }
        }

        f_out << "the V velocity field" << endl;

        if(myrank == 0) {   
            f_out.precision(6);
            for(int j = Ny - 1; j >= 0; j--) {
                for(int i = 0; i < Nx; i++)
                    f_out << fixed << BigV[j * Nx + i] << " ";
                f_out << endl;
            }
        }

        f_out.close();
    }
}
