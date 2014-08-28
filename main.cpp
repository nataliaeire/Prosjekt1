#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;


int tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n);

int main()
{
    int i, n;
    double h;

    n = 10;
    //vec a, b, c, f;
    vec a = zeros(n+2);
    vec b = zeros(n+2);
    vec c = zeros(n+2);
    vec f = zeros(n+2);
    vec u_real = zeros(n+2);
    h = (1.-0.)/(n+1);




    // Assign values to the elements of a, b, c
    for(i=0; i<=n+1; i++){
        if(i == 1){
            a[i] = 0;
            b[i] = 2;
            c[i] = -1;
            f[i] = 100*exp(-10*i*h);
            u_real[i] = 1-(1-exp(-10))*i*h - exp(-10*i*h);
        }if(i == n+1){
            a[i] = -1;
            b[i] = 2;
            c[i] = 0;
            f[i] = 100*exp(-10*i*h);
            u_real[i] = 1-(1-exp(-10))*i*h - exp(-10*i*h);
        }else{
            a[i] = -1;
            b[i] = 2;
            c[i] = -1;
            f[i] = 100*exp(-10*i*h);
            u_real[i] = 1-(1-exp(-10))*i*h - exp(-10*i*h);
        }
    }

    cout << u_real << endl;
    tridiagonal_matrix(a, b, c, f, n);

    return 0;
}

int tridiagonal_matrix(vec &a, vec &b, vec &c, vec &f, int n){

    int i;
    double factor;
    vec b_new=zeros(n+2);
    vec u=zeros(n+2);
    vec f_new=zeros(n+2);
    b_new = b;
    f_new = f;



    // Perform forward substitution
    factor = a[1]/b[0];

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

    cout << u << endl;

    return 0;
}

