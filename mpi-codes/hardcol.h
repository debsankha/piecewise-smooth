#ifndef HARDCOL_H_INCLUDED
#define HARDCOL_H_INCLUDED

#include <rk4.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string.h>
#include <limits>
#define TMAX 1000	//Time to evolve for bifurc diagrams
#define ORB  50		//set it to 2+max periodicity that we'll attempt to detect.(extra 1 for redundancy)
#define N 2		//dimensionality of system
#define randdouble(min,max) min+(rand()*(max-(min))*1.0)/RAND_MAX


extern float Sigma,G,F,W,m,K1;

int plottraj(vec, float);	//returns period :P
int plotpoincare(vec, double, double, double);
int plotbifurc_F(float minF, float maxF, int npts);
int detect_period(double *, double *, double *);
int plotbasin(int npts,double tmax, int rank);
void ischaos_n(float nmin,float nmax);

#endif
