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
const double lymda = 1;
const double a =0;
//double cc = 0;
//double d =1;
//const int pi = acos(-1);

const double b = 1;
const complex<double> icomp(0, 1);
double phi(double a, double b, double xi, int i);
complex<double>middlepryam1(double a1, double b1, double y);

double x(double t) {
    return(cos(t));
}

double xii(double t) {
    return(sin(t));
}

complex<double>Green(double p) {
    return(icomp * (_j0(p) + icomp * _y0(p)) / 4.0);
    //return(1.0 / (4.0 * icomp) * exp(icomp * p));

}


complex<double> K(double x, double y) {
    //return(log(abs(x - y)));
    return Green(abs(x - y));
}

complex<double> U0(double xi) {
    //complex<double> ed(1, 0);
    return (1);

}

complex<double> middlepryam(double a1, double b1) {
    double nn = 1000, h1, x1, i,c;

    if (b1 < a1) {
        c =b1;
        b1 = a1;
        a1 = c;
    }

    complex<double>in(0, 0);
    h1 = (b1 - a1) / nn;
    x1 = a1 + (h1 / 2);
    while (x1 <= b1 - (h1 / 2)) {
      //  cout<<"ya ne zahozhu"<< endl;
        in = in + (U0(x1));
        x1 = x1 + h1;
    }
   // cout << " in= " << in << endl;
    return in * h1;
}
complex<double> middlepryam1(double a1, double b1,double xi) {
    double nn = 100, h1, x1, i=0;
    complex<double>in(0, 0);
    h1 = (b1 - a1) / nn;
    x1 = a1 + (h1 / 2);
    while (x1 <= b1 - (h1 / 2)) {
        in = in + K(xi,x1);
        x1 = x1 + h1;
        i++;
    }
    //cout << " in= " << in << endl;
    return in * h1;
}

// в этой подпрограмме ошибка 
// интеграл вычисляля неправильно
// провербьте на К=1. Я написал свое интегрирование. Теперь решается и с логарифмическим ядром.
complex<double> middlepryam2(double a2, double b2, double a1, double b1) {
    double nn = 100, h2, h1, x2, x1, i,c;
    complex<double> in(0, 0);
    if (b1 < a1) {
        c = b1;
        b1 = a1;
        a1 = c;
    }
    if (b2 < a2) {
        c = b2;
        b2 = a2;
        a2 = c;
    }
    h2 = (b2 - a2) / nn;
    h1 = (b1 - a1) / nn;
    x2 = a2 + (h2 / 2);
    x1 = a1 + (h1 / 2);
    //cout<<" x= "<<x<<endl;
    /*while (x1 < b1) {
        while (x2 < b2) {

            if ((x1 - x2) != 0) {
                in = in + (K(x1, x2));
                x2 = x2 + h2;
            }
            else {
                in += 0;
                x2 = x2 + h2;
            }
        }
        x1 = x1 + h1;
    }
        */

    for (int i1 = 0; i1 < nn; i1++){
        for (int i2 = 0; i2 < nn; i2++) {
            x1 = a1 + (i1 + 0.5) * h1;
            x2 = a2 + (i2 + 0.5) * h2;
            if (abs(x1 - x2) > 1e-10) in += K(x1, x2);
        }
}
    return in * h2 * h1;

}


double del(int i, int j) {

    if (i == j) {
        return(1);
    }
    else
        return(0);


}

// Зачем пересчитывать x[j] при каждом обращении к phi? 
//  массив x[j] заполнен в начале программы... (?)
double phi(double a, double b, double xi, int i) {
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
        s = s + c[i] * phi(a,b,xi, i);

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

    double h, t[n + 1], xi[n];
    int i, j, k;

    complex<double> A[n][n + 1], c[n];
    complex<double> ed(1, 0);
    h = (b - a) / n;

    for (i = 0; i < n + 1; i++) {

        t[i] = a + i * h;
        cout << t[i] << endl;
    }

   //cout << middlepryam2(x[0], x[1], x[4], x[5]) << endl;
    //system("pause");


    for (i = 0; i < n; i++) {

        //cout<<"xi = " << xii(t[i]) << endl;
    }
    for (i = 0; i < n; i++) { 
       // cout<<"x= " << x(t[i]) << endl;
        for (j = 0; j < n; j++) {

            // cout << " j= " << j << " x (j)= " << x[j] << " x (j+1)= " << x[j + 1] << endl;
            //if ((x[i]<cc)||(x[i]>d)||(x[j] < cc) || (x[j] > d)) {
            // A[i][j] = del(i,j)*h-lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);
            //    //cout << "voshlo" << endl;
            //}
            /*else {  */
           
                A[i][j] = lymda * middlepryam2(x(t[j]), x(t[j + 1]), x(t[i]), x(t[i + 1]));
             
           /* }*/
            
        }

        //cout << "yawol= " << x(t[i])<< " loway= "<< x(t[i + 1]) << endl;
        A[i][n] = middlepryam(x(t[i]), x(t[i + 1]));
        // 1 - (xi[i] * log(abs(xi[i])) - log(abs(xi[i] - 1)) * xi[i] + log(abs(xi[i] - 1)) - 1);





     // exp(icomp * xi[i]) - (icomp * exp(icomp) - exp(icomp) * icomp * xi[i] - exp(icomp) + icomp * xi[i] + ed);
  // 
  //cos(xi[i])+icomp*sin(xi[i])- lymda *(icomp*exp(icomp)-exp(icomp)*icomp*xi[i]-exp(icomp)+icomp*xi[i]+(1.0,0.0)); 
    }

    for (i = 0; i < n; i++) {
        for (j = 0; j < n + 1; j++) {
            if (i == j) {
                A[i][j] -= 0.000;
            }

        }
        cout << endl;
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
    cout << endl;
    //double bbc = sqrt((xi[i] * xi[i]) + (x[i] * x[i]));
    for (i = 0; i < n; i++) {
        cout << u(xii(t[i]), c) << "  " << endl;
    }

    return 0;
}