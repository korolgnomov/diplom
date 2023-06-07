#include <cstdlib>
#include<iostream>
#include<math.h>
#include <complex>
#include<fstream>
//#include "complex.h"
using namespace std;

/*
*
*/
const double pi = acos(-1);
const int n = 10;
const double lymda = 1;
const double a = 0;
const double b = 0.5 * pi;
const int n2 = 2 * n;
double cc =pi ;
double d =1.5*pi;
const complex<double> icomp(0, 1);
double phi(double a, double b, double xi, int i);
//complex<double>middlepryam1(double a1, double b1, double y);

double x(double t, int num) {
	switch (num)
	{
	case 2: return cos(t);
	case 1: return t;
	}
}
double xii(double t, int num) {
	switch (num)
	{
	case 2: return sin(t);
	case 1: return t;
	}
}
double prx(double t, int num) {
	return ((x(t + 0.000001, num) - x(t - 0.000001, num)) / (2.0 * 0.000001));
}
double pry(double t, int num) {
	return((xii(t + 0.000001, num) - xii(t - 0.000001, num)) / (2.0 * 0.000001));
}
complex<double>Green(double p) {
	return (- icomp * 0.25)* (_j0(p) + icomp * _y0(p));
	//return(1.0 / (4.0 * icomp) * exp(icomp * p));

}


complex<double> K(double x1, double y1, double x2, double y2) {
	//return(log(abs(x - y)));
	double p = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	return Green(abs(p));
}

complex<double> U0(double x1, double x2) {
	//complex<double> ed(1, 0);
	return (1);

}

complex<double> middlepryam(double a1, double b1, int num) {
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
		x1 = x(t,num);
		y1 = xii(t,num);
		//cout<<"ya ne zahozhu"<< endl;
		in += (U0(x1,y1) * (sqrt(prx(t,num) * prx(t, num) + pry(t, num ) * pry(t, num))));
		// cout << " = " << in << endl;

	}
	//cout << " in= " << in <<" h1= "<<h1 << endl;
	return in * h1;
}
complex<double> middlepryam1(double x1, double y1, double a1, double b1, complex<double> phi[n],int ires, int num) { // num убрать должна быть сумма двух интегралов
	double nn = 100, h1, t,t2, i, c,h2;


	complex<double>in(0, 0);
	
	h1 = (b1 - a1) / nn;
	

	for (i = 0; i < nn; i++) {
		t = a1 + (i + 0.5) * h1;
		//x1 = x(t, num);
		//y1 = xii(t, num);
		//cout<<"ya ne zahozhu"<< endl;
		in += (K(x1, y1, x(t,num), xii(t, num)) * phi[ires] * (sqrt(prx(t, num) * prx(t, num) + pry(t, num) * pry(t, num))));
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
complex<double> middlepryam2(double a1, double b1, double a2, double b2,  int num1 , int num2) {
	double nn = 100, h2, h1, x2, x1, i=0, c,t1,t2;
	complex<double> in(0, 0);
	//if (b1 < a1) {
	//	c = b1;
	//	b1 = a1;
	//	a1 = c;
	//}
	//if (b2 < a2) {
	//	c = b2;
	//	b2 = a2;
	//	a2 = c;
	//}
	h2 = (b2 - a2) / nn;
	h1 = (b1 - a1) / nn;
	//cout << " h1= " << h1 << " h2= " << h2 << endl;
	//t2 = a2 + (i + 0.5) * h2;
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
			//if (abs(t1 - t2) > 0) cout << " x1 = " << x(t1, num1)<<" y1 = " << xii(t1, num1)<< " x2= " <<x(t2, num2)<< " y2 =" << xii(t2, num2)<< " ker=  "<<  K(x(t1, num1), xii(t1, num1), x(t2, num2), xii(t2, num2))<< " kriv1 = " <<sqrt(prx(t2, num2) * prx(t2, num2) + pry(t2, num2) * pry(t2, num2))<< endl;
			if (abs(t1 - t2) > 0) in += K(x(t1, num1), xii(t1, num1),x(t2, num2),xii(t2, num2)) * (sqrt(prx(t2, num2) * prx(t2, num2) + pry(t2, num2) * pry(t2, num2)))*(sqrt(prx(t1, num1) * prx(t1, num1) + pry(t1, num1) * pry(t1, num1)));
		}
		//cout << "kriv 2 = " << (sqrt(prx(t1, num1) * prx(t1, num1) + pry(t1, num1) * pry(t1, num1))) << endl;
	}
	//cout <<" int = " << in <<" h2= "<<h2<<" h1= "<< h1 << " return= "<< in * h2 * h1 << endl;
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
double phi(double t[n+1], double xi, int i) {
	double h, s;
	//int j;

	//h = (b - a) / n;


	// for (j = 0; j < n + 1; j++) {

	// x[j] = a + j * h;

	// }
	if ((xi >=(t[i])) && (xi < t[i + 1]))
		s = 1;

	else
		s = 0;

	return(s);

}

complex<double> u(double xi, complex<double> c[n2], double t1[n+1],double t2[n+1]) {
	int i;
	complex<double> s(0, 0);
	for (i = 0; i < 2*n; i++) {
		if (i < n) {
			s = s + c[i] * phi(t1, xi, i);
		}
		if (i >= n) {
			s = s + c[i] * phi(t2, xi, i-n);
		}
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
	
	double h,h1, t1[n+1], t2[n + 1], xi[n];
	int i, j, k;
	complex<double> A[n2][n2 + 1], c[n2],res1[n+1], res2[n+1];
	complex<double> ed(1, 0);
	ofstream  out1("1ecr.txt");
	ofstream  out2("2ecr.txt");
	ofstream  out3("3ecr.txt");
	ofstream  out4("4ecr.txt");
	h = (b - a) / n;
	h1 = (d -cc ) / n;
	// cout << "integral = " << middlepryam(0.5*pi , pi)<<" "<< 0.5 * pi<<" "<< pi << endl;
	for (i = 0; i < n+1; i++) {

		t1[i] = a + i * h;
		
	}
	for (i = 0; i < n + 1; i++) {

		t2[i] = cc + i * h1;

	}

	/*cout << middlepryam2(x[0], x[1], x[4], x[5]) << endl;
	system("pause");*/
	for (i = 0; i < n+1; i++) {

		cout << t1[i] << endl;

	}
	for (i = 0; i < n + 1; i++) {

		cout << t2[i] << endl;

	}
	for (i = 0; i < n; i++) {

		//cout<<"xi = " << xii(t[i]) << endl;
	}
	for (i = 0; i <n2; i++) {
		// cout<<"x= " << x(t[i]) << endl;
		for (j = 0; j < n2; j++) {

			// cout << " j= " << j << " x (j)= " << x[j] << " x (j+1)= " << x[j + 1] << endl;
			//if ((x[i]<cc)||(x[i]>d)||(x[j] < cc) || (x[j] > d)) {
			// A[i][j] = del(i,j)*h-lymda * middlepryam2(x[j], x[j + 1], x[i], x[i + 1]);
			// //cout << "voshlo" << endl;
			//}
			/*else { */
			if( (i < n) && (j < n)) {
				A[i][j] = lymda * middlepryam2(t1[i], t1[i + 1],t1[j], t1[j + 1] ,1,1);
			}
			if ((i < n) && (j >= n)) {
				A[i][j] = lymda * middlepryam2(t1[i], t1[i + 1],t2[j-n], t2[j-n+ 1],1,2);
			}
			if ((i >= n) && (j < n)) {
				A[i][j] = lymda * middlepryam2(t2[i-n], t2[i-n + 1], t1[j], t1[j + 1],2,1);
			}
			if ((i >= n) && (j >= n)) {
				A[i][j] = lymda * middlepryam2(t2[i-n], t2[(i-n)+ 1],t2[j-n], t2[(j-n)+1],2,2);
			}

			/* }*/

		}
		//cout << "yawol= " << x(t[i])<< " loway= "<< x(t[i + 1]) << endl;
		if ((i < n)) {
			A[i][n2] = middlepryam(t1[i], t1[i + 1], 1);
		}
		if ((i >= n)) {
			A[i][n2] = middlepryam(t2[i-n], t2[(i + 1)-n], 2);
		}
		
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
		//for (j = 0; j <2*n + 1; j++) {

		//	/*cout << A[i][j] << " ";*/


		//}
		c[i] = A[i][n2];
		/*cout << endl;*/

	}
	/*cout << endl;*/
	//double bbc = sqrt((xi[i] * xi[i]) + (x[i] * x[i]));
	
	for (i = 0; i < 2*n; i++) {
		if (i < n) {
			cout <<abs(u(t1[i], c, t1,t2)) << " " << endl;
			res1[i] = u(t1[i], c, t1, t2);
			out1 << x(t1[i], 1) << " " << xii(t1[i], 1) << " " << abs(u(t1[i], c, t1, t2)) << endl;
		}
		if (i == n)cout << endl << "na 2-om ecrane" << endl;
		if (i >= n) {
			cout << abs(u(t2[i-n], c, t1, t2)) << " " << endl;
			res2[i-n] = u(t2[i - n], c, t1, t2);
			out1 << x(t2[i-n], 2) << " " << xii(t2[i-n], 2) << " " << abs(u(t2[i-n], c, t1, t2)) << endl;
		}
		
} 
	
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n; j++) {
	//		if ((t1[i] != t1[j])&&(t1[i]!=t1[j+1])) {
	//			V[i][j] = middlepryam1(t1[j],t1[j + 1] ,t1[i] , t1[i+1],res1, i,1);
	//		}
	//		if (t1[i] == t1[j]) {
	//			V[i][j] = middlepryam1(t1[j]+0.001, t1[j + 1], t1[i], t1[i + 1], res1, i, 1);
	//		}
	//		if (t1[i] == t1[j+1]) {
	//			V[i][j] = middlepryam1(t1[j] , t1[j + 1]+ 0.001, t1[i], t1[i + 1], res1, i, 1);
	//		}
	//		

	//	}
	//}
	//for (i = 0; i < n; i++) {
	//	for (j = 0; j < n; j++) {
	//		if ((t2[i] != t2[j]) && (t2[i] != t2[j + 1])) {
	//			V[i+n][j+n] = middlepryam1(t2[j], t2[j + 1], t2[i], t2[i + 1], res2, i, 2);
	//		}
	//		if (t2[i] == t2[j]) {
	//			V[i+n][j+n] = middlepryam1(t2[j] + 0.001, t2[j + 1], t1[i], t1[i + 1], res2, i, 2);
	//		}
	//		if (t2[i] == t2[j + 1]) {
	//			V[i+n][j+n] = middlepryam1(t2[j], t2[j + 1] + 0.001, t1[i], t1[i + 1], res2, i, 2);
	//		}
	//	}
	//}
	cout << endl << "vne" << endl;
	
	//double hvne2 = (0 -(-1)) / 100;
	complex<double> V=0;
	double tperv1, tperv2, tvtr1, tvtr2,nach=-2,kon=2;
	double hvne1 = (kon-nach)/100;
	for (k =0; k < 100; k++) {
	tperv1 = nach + k * hvne1;
	for (j = 0; j < 100; j++) {
		tperv2 = nach + j* hvne1;
		//cout << tperv2 << endl;
		for (i = 0; i < n2; i++) { 

			if (i < n) {
						//if (tperv2 == t1[i]) tperv2 + 0.00001;
						V += abs(middlepryam1(tperv1, tperv2, t1[i], t1[i+1], res1, i, 1));
						//cout << abs(middlepryam1(tperv1, tperv2, t1[i], t1[i+1], res1, i, 1)) << endl;;

			}
			if (i >= n) {
				//cout << "pishu zaebal";
						V += abs(middlepryam1(tperv1, tperv2, t2[i-n], t2[i-n+1], res2, i-n, 2));
						//cout << abs(middlepryam1(x(tvtr1, 2), x(tvtr2, 2), t2[i], t2[i+ 1], res2, i, 2)) << endl;;
			
			}	
		}

		out3 << tperv1 << " " << tperv2<< " " << abs(V) << endl;
		V = 0;
		}
}
		////if(i==n)cout << endl << "na 2-om ecrane" << endl;
		//for (i = 0; i < n; i++) {
		//	if(i<n-1){
		//	for (j = 0; j < 1000; j++) {
		//		
		//		tvtr1 = cc + j * hvne1;
		//		tvtr2 = d + (j + 1) * hvne1;
		//		if (tvtr1 == t2[i]) tvtr1 + 0.00001;
		//		if (tvtr2 == t2[i]) tvtr2 + 0.00001;
		//		V[i+n] = abs(middlepryam1(x(tvtr1, 2), x(tvtr2, 2), t2[i], t2[i + 1], res2, i, 2));
		//		cout << abs(middlepryam1(x(tvtr1, 2), x(tvtr2, 2), t2[i], t2[i+ 1], res2, i, 2)) << endl;;
		//		out4 << x(tvtr1, 2) << " " << xii(t2[i], 2) << " " << abs(V[i + n]) << endl;
		//	}
		//	}
		//}
		//
		//cout << endl;
	

	//cout << "doshel" << endl;
	system("pause");
	return 0;
}