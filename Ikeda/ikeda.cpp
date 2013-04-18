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
		zold=z;
		z=next(z,u);
		
//		if(z==zold)		//Detecting if fixed pt conv. is occuring or not
//		{
//			cout<<real(z)<<'\t'<<imag(z)<<endl;

//			break;
//		}
//		else
//		{
//			1==1;
//			cerr<<real(zorig)<<'\t'<<imag(zorig)<<endl; //The pts which does NOT lead to f.p. the basin for chao attr
//		}
		cout<<real(z)<<'\t'<<imag(z)<<endl;
	}
	}
}


main(int argc, char **argv)
{
	float u=atof(argv[1]);
	int ninit=atoi(argv[2]);
//	complex<double> az(3,5);
//	az=next(az,u);
//	cout<<az<<endl;
//	cout<<real(az)<<'\t'<<imag(az)<<'\t'<<abs(az)<<endl;
	plot_attractor(u,ninit);
//
//	cout<<pow(abs(exp(J*3.14)),2)<<endl;
}
