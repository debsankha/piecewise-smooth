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

int plottraj(vec, float);	//returns period :P
int plotpoincare(vec, double, double, double);
int plotbifurc_F(float minF, float maxF, int npts);
int detect_period(double *, double *, double *);
double time_to_stabilize(vec x, double tmax);

int plotmap(int npts, int n, float newF, float newG);

#endif
