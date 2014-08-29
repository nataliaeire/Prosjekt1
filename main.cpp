#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;


vec tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n);
void create_vectors(int n, vec &a, vec &b, vec &c, vec &f, vec &u_real);

int main()
{
    int n_exponent, n;
    double start, finish, operation_time, relative_error;

    // Creating a file to write the results to

    ofstream myfile;
    myfile.open("tekstfil_prosjekt1.dat");
    myfile.setf(ios::scientific);               // Forcing scientific notation

    // Creating for-loop to try the solving mechanism for different size matrices
    for(n_exponent=1; n_exponent<8; n_exponent++){

        n = pow(10,n_exponent);                 // Dimension of matrix
        vec a, b, c, f, u_real;                 // Please see create_vectors for explanation
        vec u = zeros(n+2);                     // Vector to store the solution of the problem

        create_vectors(n, a, b, c, f, u_real);  // Creating vectors for the tridiagonal matrix system

        start = clock();                        // Starting clock
        u = tridiagonal_matrix(a, b, c, f, n);  // Solving the tridiagonal matrix problem
        finish = clock();                       // Ending clock

        operation_time = (finish - start)/(double) CLOCKS_PER_SEC;  // Calculating time in seconds
        relative_error = log10(max(abs((u_real(span(1,n))-u(span(1,n)))/u_real(span(1,n)))));

        // Writing results to file
        if(n_exponent == 1){                    // Setting information about the columns in the first row
            myfile << setiosflags(ios::showpoint | ios:: uppercase);
            myfile << setw(10) << setprecision(8) << "n-exponent " << "log10 rel.error" << "time" << endl;
        }
        myfile << setiosflags(ios::showpoint | ios:: uppercase);
        myfile << setw(10) << setprecision(8) << n_exponent << " "
               << setw(10) << setprecision(8) << relative_error << " "
               << setw(10) << setprecision(8) << operation_time << endl;

    } // End for loop over n values

    myfile.close();
    return 0;

} // End of program

void create_vectors(int n, vec &a, vec &b, vec &c, vec &f, vec &u_real){

    // Resizing the vector to size needed for the current for-loop
    a.resize(n+2);              // Vector storing elements a_ij (j = i-1) for tridiagonal matrix A
    b.resize(n+2);              // Vector storing elements a_ii           for tridiagonal matrix A
    c.resize(n+2);              // Vector storing elements a_ij (j = i+1) for tridiagonal matrix A
    f.resize(n+2);              // Vector storing elements in vector f in the problem Ax = f
    u_real.resize(n+2);         // Vector storing the values of the actual solution to the problem

    // Filling all elements of the vectors with zeros
    a.zeros();
    b.zeros();
    c.zeros();
    f.zeros();
    u_real.zeros();

    double h = (1.-0.)/(n+1);   // Step size of variable

    // Assign values to the elements of a, b, c
    for(int i=0; i<=n+1; i++){

        double x = i*h;         // Variable of functions

        if(i == 1){
            a[i] = 0;
            b[i] = 2;
            c[i] = -1;
            f[i] = h*h*100*exp(-10*x);
            u_real[i] = 0;
        }if(i == n+1){
            a[i] = -1;
            b[i] = 2;
            c[i] = 0;
            f[i] = h*h*100*exp(-10*x);
            u_real[i] = 0;
        }else{
            a[i] = -1;
            b[i] = 2;
            c[i] = -1;
            f[i] = h*h*100*exp(-10*x);
            u_real[i] = 1-(1-exp(-10))*x - exp(-10*x);
        } // End if statements
    } // End vector value assignment
} // End the create_vectors-function

vec tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n){

    int i;
    double factor;
    vec b_new=zeros(n+2);       // Vector storing modified values of vector b after Gaussian elimination
    vec f_new=zeros(n+2);       // Vector storing modified values of vector f after Gaussian elimination
    vec u=zeros(n+2);           // Vector storing solution to the problem

    // Initialise
    b_new = b;
    f_new = f;

    // Perform forward substitution
    factor = a[2]/b[1];

    for(i=2; i<=n; i++){
        b_new[i] = b[i] - factor*c[i-1];
        f_new[i] = f[i] - factor*f_new[i-1];
        factor = a[i+1]/b_new[i];
    }

    // Perform backward substitution
    u[n] = f_new[n]/b_new[n];

    for(i=n-1; i>0; i--){
        u[i] = (f_new[i]-c[i]*u[i+1])/b_new[i];
    }

    // cout << u << endl;

    /* Note that in the special case of a tridiagonal matrix of second dervatives, I could've just set
     * [c-1] = -1, so that b_new[i] = b[i} + factor, as well as u[i] = (f_new[i]+u[i+1])/b_new[i],
     * which would've reduced the number of operations by 2n */

    return u;
}

