#include<iostream>
#include<stdlib.h>
#include<cmath>
#define N 3
#define h 0.01
class vec
{
	public:
		double *arr;
		vec();
		vec(double arg[]);
		vec operator +(vec other);
		vec operator %(double other);
		vec operator /(double other);
		void show();
};

void f(double, vec, vec*);	
void rk4(double t, vec *x);
