
#include <cstdlib>
#include<iostream>
#include<math.h>
#include <complex>

//#include "complex.h"
using namespace std;

/*
 *
 */
const int n = 20;
const double lymda = 0.05;
const double a = 0;
const double b = 1;
const complex<double> icomp(0, 1);
complex<double> K(double x, double y) {

    return(log(abs(x - y)));



}
complex<double> U0(double xi) {
    complex<double> ed(1, 0);
    return (1 - lymda * (xi * log(abs(xi)) - log(abs(xi - 1)) * xi + log(abs(xi - 1)) - 1));

}

complex<double> middlepryam(double a1, double b1) {
    double nn = 10000, h, x, i;
    complex<double>in(0, 0);
    h = (b1 - a1) / nn;
    x = a1 + (h / 2);
    while (x <= b1 - (h / 2)) {
        in = in + (U0(x));
        x = x + h;


    }

    return in * h;
}

complex<double>Green(double p) {
    return(icomp*(_j0(p)+icomp*_y0(p))/4.0);
    //return(1.0 / (4.0 * icomp) * exp(icomp * p));

}
complex<double> middlepryam2(double a2, double b2, double a1, double b1) {
    double nn = 1000, h, h1, x, x1, i;
    complex<double> in(0, 0);
    h = (b2 - a2) / nn;
    h1 = (b1 - a1) / nn;
    x = a2 + (h / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    while (x1 < b1) {
        while (x < b2) {

            if ((x1 - x) != 0) {
                in = in + (K(x1, x));
                x = x + h;
            }
            else {
                in += 0;
                x = x + h;
            }
        }
        x1 = x1 + h1;
    }

    return in * h*h1;

}


double del(int i, int j) {

    if (i == j) {
        return(1);
    }
    else
        return(0);


}

double phi(double xi, int i) {
    double x[n + 1], h, s;
    int j;

    h = (b - a) / n;


    for (j = 0; j < n + 1; j++) {

        x[j] = a + j * h;

    }
    if ((xi >= x[i]) && (xi <= x[i + 1]))
        s = 1;

    else
        s = 0;

    return(s);

}

complex<double> u(double xi, complex<double> c[n]) {
    int i;
    complex<double> s(0, 0);
    for (i = 0; i < n; i++) {
        s = s + c[i] * phi(xi, i);

    }
    return(s);
}

void Gauss(int k, complex<double> Matrix[n][n + 1]) {

    if (Matrix[k][k] != (1.0, 0.0)) {
        complex<double> T = Matrix[k][k];
        for (int j = k; j < n + 1; j++) {
            Matrix[k][j] = Matrix[k][j] / T;
        }
    }
    for (int i = 0; i < n; i++) {
        if ((Matrix[i][k] != (0.0, 0.0)) && (i != k)) {
            complex<double> T = Matrix[i][k];
            Matrix[i][k] = (0, 0);
            for (int j = k + 1; j < n + 1; j++) {
                Matrix[i][j] -= Matrix[k][j] * T;
            }
        }
    }
    if (k < n - 1) {
        Gauss(k + 1, Matrix);
    }
}

int main(int argc, char** argv) {

    double h, x[n + 1], xi[n];
    int i, j, k;

    complex<double> A[n][n + 1], c[n];
    complex<double> ed(1, 0);
    h = (b - a) / n;

    for (i = 0; i < n + 1; i++) {

        x[i] = a + i * h;
        cout << x[i] << endl;
    }
    for (i = 0; i < n; i++) {

        xi[i] = x[i] + (h / 2);
        cout << xi[i] << endl;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            // cout<<" j= "<<j<<" x (j)= "<<x[j]<<" x (j+1)= "<<x[j+1]<<endl;
            A[i][j] = del(i, j)*h - lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);
        }

        // cout << "yawol= " << xi[i] << endl;
        A[i][n] = A[i][n] = middlepryam(x[i], x[i + 1]);
           // 1 - (xi[i] * log(abs(xi[i])) - log(abs(xi[i] - 1)) * xi[i] + log(abs(xi[i] - 1)) - 1);





        // exp(icomp * xi[i]) - (icomp * exp(icomp) - exp(icomp) * icomp * xi[i] - exp(icomp) + icomp * xi[i] + ed);
     // 
     //cos(xi[i])+icomp*sin(xi[i])- lymda *(icomp*exp(icomp)-exp(icomp)*icomp*xi[i]-exp(icomp)+icomp*xi[i]+(1.0,0.0)); 
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    Gauss(0, A);

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {

            cout << A[i][j] << " ";


        }
        c[i] = A[i][n];
        cout << endl;
    }
    //double bbc = sqrt((xi[i] * xi[i]) + (x[i] * x[i]));
    for (i = 0; i < n; i++) {
        cout << u(xi[i], c) << "  " << cos(xi[i]) + icomp * sin(xi[i]) << "  " << Green(xi[i] - x[i])<< 1.0 / (4.0 * icomp) * exp(icomp *(xi[i]-x[i])) << endl;
    }

    return 0;
}
