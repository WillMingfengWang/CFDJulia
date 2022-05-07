#include <iostream>
#include <cmath>
using namespace std;
#include <fstream>
#define PI 3.1415926

int main() {

    double x_l = -1.0; // left border of compute domain
    double x_r = 1.0; // right border of compute domain
    double dx = 0.025; // grid size
    int nx = (int)((x_r-x_l)/dx); // the number of grid which equals to the number of points plus one

    double dt = 0.0025; // time step
    double t = 1.0; // final time
    int nt = (int) (t/dt); // the number of time steps between initial time and final time

    double a = 1.0/(PI * PI); // coefficient alpha

    double x[nx+1];
    double un_e[nx+1]; // matrix for exact value
    double un[nt+1][nx+1]; // two-dimensional matrix for u
    double ut[nt+1][nx+1]; // temporary array

    for(int i =0; i<nx+1;i++) // initial condition
    {
        x[i]=x_l + dx * i;
        un[0][i] = - sin(PI * x[i]);
        un_e[i]=-exp(-t) * sin(PI*x[i]);
    }

    for (int i = 0; i < nt+1; i++)
    {
        un[i][0] = un[i][nx] = ut[i][0] = ut[i][nx] = 0; // dirichlet boundary condition
    }


    double beta;
    beta = (a * dt)/(dx * dx);

    for(int i = 1; i< nt+1;i++) {

        for (int j = 1; j < nx; j++) // 1st step
        {
            ut[i][j] = un[i - 1][j] + beta * (un[i - 1][j + 1] - 2 * un[i - 1][j] + un[i - 1][j - 1]);
        }

        for (int k = 1; k < nx; k++) // 2ed step
        {
            ut[i][k] = 0.75 * un[i-1][k] + 0.25 * ut[i][k] +
                       0.25 * beta * (ut[i][k+1] - 2*ut[i][k] + ut[i][k-1]);
        }

        for (int l = 1; l < nx; l++) // 3rd step
        {
            un[i][l] = 0.33333333 * un[i-1][l] + 0.66666667*ut[i][l] +
                       0.66666667 * beta * (ut[i][l+1] - 2*ut[i][l] + ut[i][l-1]);
        }

    }


    /*
    for(int i = 0; i < nt +1; i ++)
    {
        for (int j =0; j < nx+1; j++)
        {
            cout << j << "\t"<< un[i][j] << "\t";
        }
        cout << endl;
    }
    */

    double un_final[nx+1];
    double un_err[nx+1];

    for(int i = 0; i < nx+1; i++)
    {
        un_final[i] = un[nt][i];
        un_err[i] = un_final[i] - un_e[i];

        cout << x_l + i * dx << "\t";
        cout << un_final[i] << "\t";
        cout << un_e[i] << "\t";
        cout << un_err[i] << "\n";
    }

    return 0;
}