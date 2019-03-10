#include <bits/stdc++.h>
using namespace std;
class Ellipse{
public:
    double a,b,c;
    Ellipse(double p1, double p2, double p3){
        a=p1;
        b=p2;
        c=p3;

    };
    vector<vector<double>> cords(double u, double v, int S){
        vector<vector<double>> T;
        double C[2][3];
        double t ;
        int surf = S;
        // z constant
        if(surf==1){
             t= pow(u/a,2)+pow(v/b,2)+1.0/3.0;
             C[0][0] = u/pow(t,0.5);
             C[0][1] = v/pow(t,0.5);
             C[0][2] = c*pow(3.0,-0.5)/pow(t,0.5);
            // 180 rot
             C[1][0] = C[0][0];
            C[1][1] = -C[0][1];
            C[1][2] = -C[0][2];

        }
        // x const
        if(surf==2){
             t= pow(u/c,2)+pow(v/b,2)+1.0/3.0;
             C[0][0] = a*pow(3.0,-0.5)/pow(t,0.5);
             C[0][1] = v/pow(t,0.5);
             C[0][2] = u/pow(t,0.5);
            // 180 rot
             C[1][0] = -C[0][0];
            C[1][1] = C[0][1];
            C[1][2] = -C[0][2];

        }
        // y const
        if(surf==3){
             t= pow(u/a,2)+pow(v/c,2)+1.0/3.0;
             C[0][0] = u/pow(t,0.5);
             C[0][1] = b*pow(3.0,-0.5)/pow(t,0.5);
             C[0][2] = v/pow(t,0.5);
            // 180 rot
             C[1][0] = C[0][0];
            C[1][1] = -C[0][1];
            C[1][2] = -C[0][2];

        }

        // -90 around x


        vector<double> temp;
        for(int i=0;i<2;i++){
            for (int j=0;j<3;j++){
                temp.push_back(C[i][j]);
            }
            T.push_back(temp);
            temp.clear();
        }
        return T;

    }
    double jcb(double u, double v, int S){
        int surf = S;
        double res,t,a1,a2,a3,b1,b2,b3;

        if(surf==1){
            t = pow(u/a,2)+pow(v/b,2)+1.0/3.0;
            a1= (pow(v/b,2)+1/3.0)/(pow(t,1.5));
            a2 = -u*v*pow(a,-2)/pow(t,1.5);
            a3 = -u*c*pow(a,-2)/(pow(t,1.5)*sqrt(3));
            b1 = -u*v*pow(b,-2)/pow(t,1.5);
            b2 = (pow(u/a,2)+1/3.0)/(pow(t,1.5));
            b3 = -v*c*pow(b,-2)/(pow(t,1.5)*sqrt(3));
            res = sqrt(pow(a2*b3-a3*b2, 2) + pow(a3*b1-a1*b3, 2) + pow(a1*b2-a2*b1, 2));

        }
        if(surf==2){
            t = pow(u/c,2)+pow(v/b,2)+1.0/3.0;
            a1= (pow(v/b,2)+1/3.0)/(pow(t,1.5));
            a2 = -u*v*pow(c,-2)/pow(t,1.5);
            a3 = -u*a*pow(c,-2)/(pow(t,1.5)*sqrt(3));
            b1 = -u*v*pow(b,-2)/pow(t,1.5);
            b2 = (pow(u/c,2)+1/3.0)/(pow(t,1.5));
            b3 = -v*a*pow(b,-2)/(pow(t,1.5)*sqrt(3));
            res = sqrt(pow(a2*b3-a3*b2, 2) + pow(a3*b1-a1*b3, 2) + pow(a1*b2-a2*b1, 2));

        }
        if(surf==3){
            t = pow(u/a,2)+pow(v/c,2)+1.0/3.0;
            a1= (pow(v/c,2)+1/3.0)/(pow(t,1.5));
            a2 = -u*v*pow(a,-2)/pow(t,1.5);
            a3 = -u*b*pow(a,-2)/(pow(t,1.5)*sqrt(3));
            b1 = -u*v*pow(c,-2)/pow(t,1.5);
            b2 = (pow(u/a,2)+1/3.0)/(pow(t,1.5));
            b3 = -v*b*pow(c,-2)/(pow(t,1.5)*sqrt(3));
            res = sqrt(pow(a2*b3-a3*b2, 2) + pow(a3*b1-a1*b3, 2) + pow(a1*b2-a2*b1, 2));

        }

        return res;

    }
};


double pot(double u, double v, double x1, double x2, double x3, int S, double s1, double s2, double s3){
    Ellipse E(s1,s2,s3);
    vector<vector<double>> Cs = E.cords(u,v,S);
    double J = E.jcb(u,v,S);
    double y1,y2,y3;
    double res=0;
    for(int i=0;i<2;i++){
        y1=Cs[i][0];
        y2=Cs[i][1];
        y3=Cs[i][2];
        double t = sqrt(pow(x1-y1,2)+pow(x2-y2,2)+pow(x3-y3,2));
        double den = sqrt(pow(2*y1*pow(s1,-2),2)+pow(2*y2*pow(s2,-2),2)+pow(2*y3*pow(s3,-2),2));
        res+= ((x1-y1)*y1*pow(s1,-2) + (x2-y2)*y2*pow(s2,-2) + (x3-y3)*y3*pow(s3,-2))*J*2/(pow(t,3)*den);
        //res+=J;
    }
    return res;
}

double nu_k(int k){
    double res = (k%2==0)? -2.0/(k*k - 1) : 0;
    return res;
}

double chebyshev_semi(int n, double a, double b, double y, double x1, double x2, double x3, int S, double s1, double s2, double s3){
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
        res += w_cc[j]*pot(x_points[j], y, x1, x2, x3,S,s1,s2,s3);
    }
    return res;
}

int main()
{
    double x1= 1;
    double x2=0;
    double x3 =0;
    int rat[3]={2,6,3};
    int N = 64;
    double sides[3]={1,2,3};
    double pi = 3.14159265358979323846;
    double a, b, c, d;
    double res = 0;
    int n;
    for(int itr=0;itr<3;itr++){
        if(itr==0){
            a =-sides[0]/sqrt(3);
            b = sides[0]/sqrt(3);
            c = -sides[1]/sqrt(3);
            d = sides[1]/sqrt(3);
        }
        else if(itr==1){
            a =-sides[2]/sqrt(3);
            b = sides[2]/sqrt(3);
            c = -sides[1]/sqrt(3);
            d = sides[1]/sqrt(3);
        }
        else {
            a =-sides[0]/sqrt(3);
            b = sides[0]/sqrt(3);
            c = -sides[2]/sqrt(3);
            d = sides[2]/sqrt(3);
        }
        n=N*rat[itr];
        double y_points[n];

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
            res += w_cc[j]*chebyshev_semi(n, a, b, y_points[j], x1, x2, x3,itr+1,sides[0],sides[1],sides[2]);
        }

    }

    cout<< setprecision(16)<<abs(res*0.25/pi + .5)<<endl;
    return 0;
}






