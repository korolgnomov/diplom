#include <cstdlib>
#include<iostream>
#include<math.h>
#include <complex>

//#include "complex.h"
using namespace std;

/*
*
*/
const double pi = acos(-1);
const int n = 10;
const double lymda = 1;
const double a = 0;
const double b = pi/2;
const int n2 = 2 * n;
double cc = pi;
double d =1.5*pi;
const complex<double> icomp(0, 1);
double phi(double a, double b, double xi, int i);
complex<double>middlepryam1(double a1, double b1, double y);

double x(double t) {
	return(cos(t));
}
double xii(double t) {
	return(sin(t));
}
double prx(double t) {
	return ((x(t + 0.000001) - x(t - 0.000001)) / (2 * 0.000001));
}
double pry(double t) {
	return((xii(t + 0.000001) - xii(t - 0.000001)) / (2 * 0.000001));
}
complex<double>Green(double p) {
	return(icomp * (_j0(p) + icomp * _y0(p)) / 4.0);
	//return(1.0 / (4.0 * icomp) * exp(icomp * p));

}


complex<double> K(double x1, double y1, double x2, double y2) {
	//return(log(abs(x - y)));
	double p = sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
	return Green(abs(p));
}

complex<double> U0(double x1, double x2) {
	//complex<double> ed(1, 0);
	return (1);

}

complex<double> middlepryam(double a1, double b1) {
	double nn = 100, h1, x1,y1,t, i, c;

	if (b1 < a1) {
		c = b1;
		b1 = a1;
		a1 = c;
	}

	complex<double>in(0, 0);
	h1 = (b1 - a1) / nn;
	for (i = 0; i < nn; i++) {
		t = a1 + (i + 0.5)* h1;
		x1 = x(t);
		y1 = xii(t);
		//cout<<"ya ne zahozhu"<< endl;
		in += (U0(x1,y1) * (sqrt(prx(t) * prx(t) + pry(t) * pry(t))));
		// cout << " = " << in << endl;

	}
	//cout << " in= " << in <<" h1= "<<h1 << endl;
	return in * h1;
}
//complex<double> middlepryam1(double a1, double b1, double xi) {
//	double nn = 100, h1, x1, i = 0;
//	complex<double>in(0, 0);
//	h1 = (b1 - a1) / nn;
//	x1 = a1 + (h1 / 2);
//	while (x1 <= b1 - (h1 / 2)) {
//		in = in + K(xi, x1);
//		x1 = x1 + h1;
//		i++;
//	}
//	//cout << " in= " << in << endl;
//	return in * h1;
//}

// в этой подпрограмме ошибка
// интеграл вычисляля неправильно
// провербьте на К=1. Я написал свое интегрирование. Теперь решается и с логарифмическим ядром.
complex<double> middlepryam2(double a2, double b2, double a1, double b1) {
	double nn = 100, h2, h1, x2, x1, i=0, c,t1,t2;
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
	
	t2 = a2 + (i + 0.5) * h2;
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

	for (int i1 = 0; i1 < nn; i1++) {
		t1 = a1 + (i1 + 0.5) * h1;
		for (int i2 = 0; i2 < nn; i2++) {
			t2 = a2 + (i2 + 0.5) * h2;
			if (abs(t1 - t2) > 1e-10) in += K(x(t1), xii(t1),x(t2),xii(t2)) * (sqrt(prx(t2) * prx(t2) + pry(t2) * pry(t2)));
		}
		in *= (sqrt(prx(t1) * prx(t1) + pry(t1) * pry(t1)));
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
// массив x[j] заполнен в начале программы... (?)
double phi(double t[n2], double xi, int i) {
	double h, s;
	//int j;

	//h = (b - a) / n;


	// for (j = 0; j < n + 1; j++) {

	// x[j] = a + j * h;

	// }
	if ((xi >=(t[i])) && (xi <= t[i + 1]))
		s = 1;

	else
		s = 0;

	return(s);

}

complex<double> u(double xi, complex<double> c[n2], double t[n2]) {
	int i;
	complex<double> s(0, 0);
	for (i = 0; i < n2; i++) {
		s = s + c[i] * phi(t, xi, i);

	}
	return(s);
}
void Gauss(int k, complex<double> Matrix[n2][n2+1]) {
	if (Matrix[k][k] != (1.0, 0.0)) {
		complex<double> T = Matrix[k][k];
		for (int j = k; j < n2 + 1; j++) {
			Matrix[k][j] = Matrix[k][j] / T;
		}
	}
	for (int i = 0; i < n2; i++) {
		if ((Matrix[i][k] != complex<double>(0.0, 0.0)) && (i != k)) {
			complex<double> T = Matrix[i][k];
			Matrix[i][k] = complex<double>(0.0, 0.0);
			for (int j = k + 1; j < n2 + 1; j++) {
				Matrix[i][j] -= Matrix[k][j] * T;
			}
		}
	}
	if (k < n2 - 1) {
		Gauss(k + 1, Matrix);
	}
}
int main(int argc, char** argv) {

	double h,h1, t[2*n+2], xi[n];
	int i, j, k;
	complex<double> A[n2][n2 + 1], c[n2];
	complex<double> ed(1, 0);
	h = (b - a) / n;
	h1 = (d -cc ) / n;
	// cout << "integral = " << middlepryam(0.5*pi , pi)<<" "<< 0.5 * pi<<" "<< pi << endl;
	for (i = 0; i < n+1; i++) {

		t[i] = a + i * h;
		
	}
	for (i = 0; i < n+1; i++) {

		t[i+n] = cc + i * h1;
		
	}
	//cout << middlepryam2(x[0], x[1], x[4], x[5]) << endl;
	//system("pause");
	for (i = 0; i < 2*n; i++) {

		cout << t[i] << endl;

	}

	for (i = 0; i < n; i++) {

		//cout<<"xi = " << xii(t[i]) << endl;
	}
	for (i = 0; i <2*n; i++) {
		// cout<<"x= " << x(t[i]) << endl;
		for (j = 0; j < 2*n; j++) {

			// cout << " j= " << j << " x (j)= " << x[j] << " x (j+1)= " << x[j + 1] << endl;
			//if ((x[i]<cc)||(x[i]>d)||(x[j] < cc) || (x[j] > d)) {
			// A[i][j] = del(i,j)*h-lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);
			// //cout << "voshlo" << endl;
			//}
			/*else { */
			if( (i < n) && (j < n)) {
				A[i][j] = lymda * middlepryam2(t[j], t[j + 1], t[i], t[i + 1]);
			}
			if ((i < n) && (j > n)) {
				A[i][j] = lymda * middlepryam2(t[i], t[i + 1],t[j], t[j + 1]);
			}
			if ((i > n) && (j < n)) {
				A[i][j] = lymda * middlepryam2(t[i], t[i + 1], t[j], t[j + 1]);
			}
			if ((i > n) && (j > n)) {
				A[i][j] = lymda * middlepryam2(t[j], t[j + 1],t[i], t[i + 1]);
			}

			/* }*/

		}
		//cout << "yawol= " << x(t[i])<< " loway= "<< x(t[i + 1]) << endl;
		A[i][n2] = middlepryam(t[i], t[i + 1]);
		// 1 - (xi[i] * log(abs(xi[i])) - log(abs(xi[i] - 1)) * xi[i] + log(abs(xi[i] - 1)) - 1);
		// exp(icomp * xi[i]) - (icomp * exp(icomp) - exp(icomp) * icomp * xi[i] - exp(icomp) + icomp * xi[i] + ed);
		//
		//cos(xi[i])+icomp*sin(xi[i])- lymda *(icomp*exp(icomp)-exp(icomp)*icomp*xi[i]-exp(icomp)+icomp*xi[i]+(1.0,0.0));
	}

	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n + 1; j++) {
	//		if (i == j) {
	//			A[i][j] -= 0.000;
	//		}

	//	}
	//	cout << endl;
	//}

	for (i = 0; i < 2*n; i++) {
		for (j = 0; j < 2*n + 1; j++) {
			cout << A[i][j] << " ";
		}
		cout << endl;
	}


	Gauss(0, A);

	for (i = 0; i <2* n; i++) {
		for (j = 0; j <2*n + 1; j++) {

			cout << A[i][j] << " ";


		}
		c[i] = A[i][n2];
		cout << endl;

	}
	cout << endl;
	//double bbc = sqrt((xi[i] * xi[i]) + (x[i] * x[i]));
	for (i = 0; i < 2*n; i++) {
		cout << u(t[i], c, t) << " " << endl;
	}

	return 0;
}