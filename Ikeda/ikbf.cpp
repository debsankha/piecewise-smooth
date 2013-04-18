#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<complex>
#define NUM_ITER 1000
#define NUM_DISCARD 10000
#define NUM_KEEP 1024
#define xmin -30
#define xmax 30

using namespace std;
complex<double> J(0,1);

complex<double> next(complex<double> z,double u)
{

	extern complex<double> J;
	complex<double> zn;
	return z*exp(J*(complex<double> (0.4-6/(pow(abs(z),2)+1))))*u+(complex<double>)1;
}

void plot_attractor(double u, int ninit)
{

	srand(time(NULL));
	for (int nit=0;nit<ninit;nit++)				//do for 20000 diff init. point.
	{
	complex<double> z(xmin+((float)rand()/RAND_MAX)*(xmax-xmin),xmin+((float)rand()/RAND_MAX)*(xmax-xmin));
	complex<double> zold;
	
	complex<double> zorig=z;

	for (int i=0;i<NUM_DISCARD;i++){z=next(z,u);}	//discard 1st few iter
	for (int i=0;i<NUM_KEEP;i++)
	{
		z=next(z,u);
		cout<<u<<'\t'<<real(z)<<'\t'<<imag(z)<<endl;
	}
	}
}


main(int argc, char **argv)
{
	double umin=atof(argv[1]);
	double umax=atof(argv[2]);

	double du=(umax-umin)/1000;

	for (double u=umin;u<umax;u+=du)
	{
		plot_attractor(u,100);
	}
}
