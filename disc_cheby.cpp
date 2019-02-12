#include <bits/stdc++.h>

using namespace std;
double r = 1.0;
double func(double u, double v){
	double s = r/2.0;
	double pi = 3.14159265358979323846;
	double beta = u*r + s*(1-u)/abs(sin(pi/4.0+v*pi/2.0));
	double a1 = (r - 0.5/abs(sin(pi/4.0+ v*pi/2.0)))*sin(v*pi/2.0);
	double a2 = (r - 0.5/abs(sin(pi/4.0 + v*pi/2.0)))*cos(v*pi/2.0);
	double a3 = 0.0;
	double t=-1*s*(1-u)*(pi/2.0)/(sin(pi/4.0+v*pi/2.0)*tan(pi/4.0+v*pi/2.0));
	if (v>0.5){
		t*=-1;
	}
	double b1 = t*sin(v*pi/2.0) + pi*beta*cos(v*pi/2.0)/2.0;
	double b2 = t*cos(v*pi/2.0) - pi*beta*sin(v*pi/2.0)/2.0;
	double test = 1*s*(1-u)*(pi/2.0)*cos(v*pi/2.0)/(sin(pi/4.0+v*pi/2.0)*tan(pi/4.0+v*pi/2.0));
	//cout << test << endl;
	double b3 = 0.0;
	double res = sqrt(pow(a2*b3-a3*b2, 2) + pow(a3*b1-a1*b3, 2) + pow(a1*b2-a2*b1, 2));
	return res;
	//return 1.0/sqrt(1-u*u-v*v);
}

double nu_k(int k){
    double res = (k%2==0)? -2.0/(k*k - 1) : 0;
    return res;
}

double chebyshev_semi(int n, double a, double b, double y){
    double x_points[n];
    double pi = 3.14159265358979323846;
    for(int i=0;i<n;i++){
        x_points[i] = a + (b-a)*(1+cos((i+1/2.0)*pi/n))/2.0;
    }
    double w_cc[n];
    for(int j=0;j<n;j++){
        w_cc[j] = 1;
        for(int k=1;k<=n;k++){
            w_cc[j] += nu_k(k)*cos(k*pi*(j+0.5)/n);
        }
        w_cc[j] *= (b-a)/n;
    }
    double res = 0.0;
    for(int j=0;j<n;j++){
        res += w_cc[j]*func(x_points[j], y);
    }
    return res;
}

int main()
{
    int n;
    n = 16;
    double a, b, c, d;
    a =0;
    b = 1;
    c = 0;
    d = 1;

    double res = 0;
    double y_points[n];
    double pi = 3.14159265358979323846;
    for(int i=0;i<n;i++){
        y_points[i] = c + (d-c)*(1+cos((i+1/2.0)*pi/n))/2.0;
    }
    double w_cc[n];
    for(int j=0;j<n;j++){
        w_cc[j] = 1;
        for(int k=1;k<=n;k++){
            w_cc[j] += nu_k(k)*cos(k*pi*(j+0.5)/n);
        }
        w_cc[j] *= (d-c)/n;
    }
    for(int j=0;j<n;j++){
        res += w_cc[j]*chebyshev_semi(n, a, b, y_points[j]);
    }
    cout<< setprecision(16)<< abs(4*res + 1  - pi)<<endl;
    return 0;
}
