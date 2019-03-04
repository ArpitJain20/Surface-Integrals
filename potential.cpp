#include <bits/stdc++.h>
using namespace std;
class Sphere{
public:
    vector<vector<double>> cords(double u, double v){
        vector<vector<double>> T;
        double C[6][3];
        double t = pow(u,2)+pow(v,2)+1.0/3.0;
         C[0][0] = u/pow(t,0.5);
         C[0][1] = v/pow(t,0.5);
         C[0][2] = pow(3.0,-0.5)/pow(t,0.5);
        // 90 around x
         C[1][0] = C[0][0];
         C[1][1] = -C[0][2];
         C[1][2] = C[0][1];
        // -90 around x
         C[2][0] = C[0][0];
         C[2][1] = C[0][2];
         C[2][2] = -C[0][1];
        //180 around x
         C[3][0] = C[0][0];
         C[3][1] = -C[0][1];
         C[3][2] = -C[0][2];
        //90 around y
         C[4][0] = C[0][2];
         C[4][1] = C[0][1];
         C[4][2] = -C[0][0];
        //-90 around y
         C[5][0] = -C[0][2];
         C[5][1] = C[0][1];
         C[5][2] = C[0][0];

        vector<double> temp;
        for(int i=0;i<6;i++){
            for (int j=0;j<3;j++){
                temp.push_back(C[i][j]);
            }
            T.push_back(temp);
            temp.clear();
        }
        return T;

    }
    double jcb(double u, double v){
        double res;
        double t = pow(u,2)+pow(v,2)+1.0/3.0;
        double a1= (pow(v,2)+1/3.0)/(pow(t,1.5));
        double a2 = -u*v/pow(t,1.5);
        double a3 = -u/(pow(t,1.5)*sqrt(3));
        double b1 = -u*v/pow(t,1.5);
        double b2 = (pow(u,2)+1/3.0)/(pow(t,1.5));
        double b3 = -v/(pow(t,1.5)*sqrt(3));
        res = sqrt(pow(a2*b3-a3*b2, 2) + pow(a3*b1-a1*b3, 2) + pow(a1*b2-a2*b1, 2));
        return res;

    }
};


double pot(double u, double v, double x1, double x2, double x3){
    Sphere S;
    vector<vector<double>> Cs = S.cords(u,v);
    double J = S.jcb(u,v);
    double y1,y2,y3;
    double res=0;
    for(int i=0;i<6;i++){
        y1=Cs[i][0];
        y2=Cs[i][1];
        y3=Cs[i][2];
        double t = sqrt(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3-y3,2));
        res+= ((x1-y1)*y1 + (x2-y2)*y2 + (x3-y3)*y3)*J/pow(t,3);
        //res+=J;
    }
    return res;
}

double nu_k(int k){
    double res = (k%2==0)? -2.0/(k*k - 1) : 0;
    return res;
}

double chebyshev_semi(int n, double a, double b, double y, double x1, double x2, double x3){
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
        res += w_cc[j]*pot(x_points[j], y, x1, x2, x3);
    }
    return res;
}

int main()
{
    int n;
    n = 32;
    double a, b, c, d;
    a =-1/sqrt(3);
    b = 1/sqrt(3);
    c = -1/sqrt(3);
    d = 1/sqrt(3);

    double x1=0;
    double x2=0;
    double x3 =0;

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
        res += w_cc[j]*chebyshev_semi(n, a, b, y_points[j], x1, x2, x3);
    }
    cout<< setprecision(16)<<(res)*0.25/pi<<endl;
    return 0;
}



