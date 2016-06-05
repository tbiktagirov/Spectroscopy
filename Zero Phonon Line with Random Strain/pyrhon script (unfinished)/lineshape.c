//lineshape.c
//accounts for the split by random defoemation
//the effect of energy shift not included
//adapted for a particular zpl on 1190.1917 meV
#include <stdio.h>      
#include <math.h>      
double h(int n, double args[n]){
    double x1 = (1-args[0])/args[0];
    double x2 = (1-args[1])/args[1];
    double x3 = (1-args[2])/args[2];
    double x4 = (1-args[3])/args[3];
    double freq = args[4];
    double G = args[5];
    double gamma = args[6];
    double v1 = args[7];
    double v2 = args[8];
	
    //double J = 1/((1-args[0])*(1-args[1])*(1-args[2])*(1-args[3])*(1-args[0])*(1-args[1])*(1-args[2])*(1-args[3]));
    double J = 1/(args[0]*args[1]*args[2]*args[3]*args[0]*args[1]*args[2]*args[3]);
    double delta = sqrt((x1+x2) * (x1+x2) + (x3+x4) * (x3+x4));
    double L1 = exp(-(freq - 1190.1917 + delta)*(freq - 1190.1917 + delta)/2/(G*G));
    double L2 = exp(-(freq - 1190.1917 - delta)*(freq - 1190.1917 - delta)/2/(G*G));
    double distr = pow(1/(x1*x1/v1/v1+x2*x2/v2/v2+x3*x3/v1/v1+x4*x4/v2/v2+gamma*gamma),2.5);
    double h = (L1+L2)*distr*J;
    return h;
}
