#include <iostream>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include "time.h"

using namespace std;
using namespace arma;


vec tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n);

int main()
{
    int n_exponent, i, n;
    double h, start, finish, operation_time;

    ofstream myfile;
    myfile.open("tekstfil_prosjekt1.dat");

    for(n_exponent=1; n_exponent<8; n_exponent++){
        n = pow(10,n_exponent);
        //vec a, b, c, f;
        vec a = zeros(n+2);
        vec b = zeros(n+2);
        vec c = zeros(n+2);
        vec f = zeros(n+2);
        vec u_real = zeros(n+2);
        vec u = zeros(n+2);
        h = (1.-0.)/(n+1);

        // Assign values to the elements of a, b, c
        for(i=0; i<=n+1; i++){
            if(i == 1){
                a[i] = 0;
                b[i] = 2;
                c[i] = -1;
                f[i] = h*h*100*exp(-10*i*h);
                u_real[i] = 0;
            }if(i == n+1){
                a[i] = -1;
                b[i] = 2;
                c[i] = 0;
                f[i] = h*h*100*exp(-10*i*h);
                u_real[i] = 0;
            }else{
                a[i] = -1;
                b[i] = 2;
                c[i] = -1;
                f[i] = h*h*100*exp(-10*i*h);
                u_real[i] = 1-(1-exp(-10))*i*h - exp(-10*i*h);
            } // End if statements
        } // End vector value assignment


        start = clock();

        // cout << u_real << endl;
        u = tridiagonal_matrix(a, b, c, f, n);

        finish = clock();
        operation_time = finish - start;

        if(n_exponent == 1){
            myfile << setiosflags(ios::showpoint | ios:: uppercase);
            myfile << setw(20) << setprecision(15) << "n-exponent " << "relative error (log10)" << endl;
        }

        myfile << setiosflags(ios::showpoint | ios:: uppercase);
        myfile << setw(20) << setprecision(15) << n_exponent << " " <<
                 log10(max(abs((u_real(span(1,n))-u(span(1,n)))/u_real(span(1,n))))) << endl;

    } // End for loop over n values

    return 0;
}

vec tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n){

    int i;
    double factor;
    vec b_new=zeros(n+2);
    vec u=zeros(n+2);
    vec f_new=zeros(n+2);
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

