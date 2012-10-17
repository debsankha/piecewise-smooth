#include<iostream>
#include<stdlib.h>
#include<cmath>
#define h 0.01
class vec
{
	public:
		double *arr;
		int len;
		vec(int N);
		vec(int N, double arg[]);
		vec(const vec& other);	//COPY CONST

		vec& operator =(const vec& other);
		vec operator +(const vec other);
		vec operator %(const double other);
		vec operator /(const double other);
		void show();
		~vec();
};

void f(double, vec, vec*);	
void rk4(double t, vec *x);
