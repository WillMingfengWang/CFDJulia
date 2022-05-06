#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
#define PI 3.1415926


int main() {

    double x_l = -1.0; // left border of compute domain
    double x_r = 1.0; // right border of compute domain
    double dx = 0.025; // grid size
    int nx = (int)((x_r - x_l)/dx); // the number of grid which equals to the number of points plus one

    double dt = 0.0025; // time step
    double t = 1.0; // final time
    int nt = (int) (t/dt); // the number of time steps between initial time and final time

    double a = 1.0/(PI * PI); // coefficient alpha

    double x[nx+1];
    double un_e[nx+1]; // matrix for exact value
    double un[nt+1][nx+1]; // two-dimensional matrix for u

    for(int i =0; i<nx+1;i++) // initial condition
    {
        x[i]=x_l + dx * i;
        un[0][i] = - sin(PI * x[i]);
        un_e[i]=-exp(-t) * sin(PI*x[i]);
    }
    cout << endl;
    un[0][0] = 0; // initial condition, when t equals to zero
    un[0][nx] = 0; // initial condition when t equals to zero


    double beta;
    beta = (a * dt)/(dx * dx);

    for(int i=1;i<nt+1;i++)
    {
        for(int j =1;j < nx;j++)
        {
            un[i][j] = un[i-1][j] + beta * (un[i-1][j+1] - 2.0 * un[i-1][j] + un[i-1][j-1]);
        }
        un[i][0]= 0; // boundary condition. this parameter is very important
        un[i][nx] = 0; // boundary condition. this parameter is very important as well
    }

    // print all the calculating results
    for(int i =0; i<nt+1;i++)
    {
        for(int j =0;j<nx+1;j++)
        {
            // cout << un[i][j] << "\t";
        }
        // cout << endl;
    }

    double un_final[nx+1]; // the last row of result
    for(int i = 0; i < nx+1; i++)
    {
        un_final[i] = un[nt][i];
    }
    cout << endl;

    // calculate the error between un and u_e (exact solution and numerical solution)
    double un_err[nx+1];
    for(int i = 0; i < nx+1; i++)
    {
        un_err[i] = un[nt][i] - un_e[i];
    }
    cout << endl;

    ofstream outfile("outcome.txt");

    double point[nx+1];

    for(int i = 0; i < nx+1;i++)
    {
        point [i] = -1 + i*dx;
        outfile << point[i] << "\t" << "\t";
        outfile << un_final[i] << "\t" << "\t";
        outfile << un_e[i] << "\t" << "\t";
        outfile << un_err[i] << "\n";
    }

    return 0;
}