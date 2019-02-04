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
	//return u*v*v;
}
double simpson_s(int n, double a, double b, double x){
	double del_x = (b-a)/(n);
	double res = func(a, x) + func(b, x);
	for(int i=1;i<n;i+=2){
		res += 4*func(a+i*del_x, x);
	}
	for(int i=2;i<n;i+=2){
		res += 2*func(a+i*del_x, x);
	}
	 	res = res*del_x/3;
	return res;
}
int main() {
	int n;
	double a, b, c, d;
	cin>>n;
	//n = 4;
	//cin>>n>>a>>b>>c>>d;
	a = 0;
	b = 1;
	c = a;
	d = b;
	double del_x = (b-a)/(n);
	double res = simpson_s(n, c, d, a) + simpson_s(n, c, d, b);
	for(int i=1;i<n;i+=2){
		//cout<<i<<endl;
		res += 4*simpson_s(n, c, d, a+i*del_x);
	}
	for(int i=2;i<n;i+=2){
		//cout<<i<<endl;
		res += 2*simpson_s(n, c, d, a+i*del_x);
	}
	res = res*del_x/3;
    double pi = 3.14159265358979323846;
	cout<<"For n =  "<< n<<" error is "<< setprecision(16)<< abs(4*res + 1 - pi) <<endl;
	return 0;
}
