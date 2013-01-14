#ifndef RK4_H_included
#define RK4_H_included

#include<iostream>
#include<stdlib.h>
#define h 0.001

class vec			//a vector class: contains just an array
{
	public:
		double *arr;				//R1. DYNAMIC ALLOCATION
		int len;	//the length
		vec(int N);
		vec(int N, double arg[]);	//deepcopies arg to arrr 
		vec(const vec& other);			//R2. COPY CONST
		
		vec& operator =(const vec& other);	//R3. ASSIGNMENT OVERLOADED
		vec operator +(const vec other);
		vec operator -(const vec other);
		vec operator %(const double other);//scalar multiplication
		vec operator /(const double other);
		void show();
		double norm();
		~vec();					//R4. DESTRUCTOR
};

void f(double, vec, vec*);	
void rk4(double t, vec *x);

#endif

