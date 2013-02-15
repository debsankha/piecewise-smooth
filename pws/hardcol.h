#ifndef HARDCOL_H_INCLUDED
#define HARDCOL_H_INCLUDED

#include <rk4.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string.h>
#include <limits>
#define TMAX 10000	//Time to evolve for bifurc diagrams
#define ORB  50		//set it to 2+max periodicity that we'll attempt to detect.(extra 1 for redundancy)
#define N 2		//dimensionality of system
#define randdouble(min,max) min+(max-(min))*(rand()*1.0/RAND_MAX)	


extern float Sigma,G,F,W,m,K1;

int plotpoincare(vec, double, double);
int detect_period(double *, double *, double *);
double time_to_stabilize(vec x, double tmax);
int plotmap3d(int npts, int n, float newF, float newG);
int plotmap(int npts, int n, float newF, float newG);
int plotbasin(int npts,double tmax);

#endif
