#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;

vec tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n);
void create_vectors(int n, vec &a, vec &b, vec &c, vec &f, vec &u_real);
void armadillo_solve(int n, vec &a, vec &b, vec &c, vec &f, vec &u_arma);


/* Indices in create_vectors.. they aren't correct, are they?
 * What are the correct values for i = 0, i = 1, i = n and i = n+1. */


int main()
{
    int n_exponent, n;
    double start, finish, operation_time, relative_error, start_armadillo, finish_armadillo, operation_time_armadillo;

    // Creating a file to write the time and error to
    ofstream timeanderror;
    timeanderror.open("timeanderror.dat");
    timeanderror.setf(ios::scientific);         // Forcing scientific notation

    // Creating a file to write the found function values to
    ofstream uvalues;
    uvalues.open("uvalues.dat");
    uvalues.setf(ios::scientific);              // Forcing scientific notation

    // Creating for-loop to try the solving mechanism for different size matrices
    for(n_exponent=1; n_exponent < 8; n_exponent++){

        n = pow(10,n_exponent);                 // Dimension of matrix
        vec a, b, c, f, u_real, u_arma;         // Please see sub functions for explanation
        vec u = zeros(n+2);                     // Vector to store the solution of the problem

        create_vectors(n, a, b, c, f, u_real);  // Creating vectors for the tridiagonal matrix system

        start = clock();                        // Starting clock
        u = tridiagonal_matrix(a, b, c, f, n);  // Solving the tridiagonal matrix problem
        finish = clock();                       // Ending clock

        operation_time = (finish - start)/(double) CLOCKS_PER_SEC;  // Calculating time in seconds
        relative_error = log10(max(abs((u_real(span(1,n))-u(span(1,n)))/u_real(span(1,n)))));

        // Solving problem using LU-decomposition of the matrix is small
        if(n_exponent < 2){
            start_armadillo = clock();
            armadillo_solve(n, a, b, c, f, u_arma);
            finish_armadillo = clock();
            operation_time_armadillo = (finish_armadillo - start_armadillo)/(double) CLOCKS_PER_SEC;
        }

        // Writing error and time to file
        if(n_exponent == 1){                    // Setting information about the columns in the first row
            timeanderror << setiosflags(ios::showpoint | ios:: uppercase);
            timeanderror << setw(10) << setprecision(8) << "exponent"
                         << setw(10) << setprecision(8) << "log10 rel.error"
                         << setw(10) << setprecision(8) << "time" << endl;
        }
        timeanderror << setiosflags(ios::showpoint | ios:: uppercase);
        timeanderror << setw(10) << setprecision(8) << n_exponent << " "
                     << setw(10) << setprecision(8) << relative_error << " "
                     << setw(10) << setprecision(8) << operation_time << endl;


        // Writing found function values to file
        uvalues << setiosflags(ios::showpoint | ios:: uppercase);
        uvalues << setw(10) << setprecision(8) << u << endl;

    } // End for loop over n values

    // Closing files
    timeanderror.close();
    uvalues.close();

    return 0;

} // End of program

void create_vectors(int n, vec &a, vec &b, vec &c, vec &f, vec &u_real){

    // Resizing the vector to size needed for the current for-loop
    a.resize(n+2);              // Vector storing elements a_ij (j = i-1) for tridiagonal matrix A
    b.resize(n+2);              // Vector storing elements a_ii           for tridiagonal matrix A
    c.resize(n+2);              // Vector storing elements a_ij (j = i+1) for tridiagonal matrix A
    f.resize(n+2);              // Vector storing elements in vector f in the problem Au = f
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

    /* Note that in the special case of a tridiagonal matrix of second dervatives, I could've just set
     * [c-1] = -1, so that b_new[i] = b[i} + factor, as well as u[i] = (f_new[i]+u[i+1])/b_new[i],
     * which would've reduced the number of operations by 2n */

    return u;
}


void armadillo_solve(int n, vec &a, vec &b, vec &c, vec &f, vec &u_arma){

    mat A = mat(n,n);           // Create full matris from tridiagonal matrix vectors
    A.zeros();                  // Fill matrix with zeros

    u_arma.resize(n);           // Create solution vector
    u_arma.zeros();             // Fill elements of solution vector with zeros

    vec f_arma;                 // Create vector f in the problem Au = f
    f_arma.resize(n);           // Size vector f to the right number of elements
    f_arma = f(span(1,n));      // Fill elements of f_new with elements of f

    // Assign values to elements of matrix
    for (int i = 0; i<n; i++){
        if (i == 0){
            A(i,i)   = b[i+1];
            A(i,i+1) = c[i+1];
        }if (i==n-1){
            A(i,i)   = b[i+1];
            A(i,i-1) = a[i+1];
        }else{
            A(i,i)   = b[i+1];
            A(i,i+1) = c[i+1];
            A(i,i-1) = a[i+1];
        } // End if statements
    } // End creating matrix

    cout << A << endl;
    //u_arma = solve(A,f_arma);

}// End of armadillo_solve-function
