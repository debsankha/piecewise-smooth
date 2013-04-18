#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<complex>
#define NUM_ITER 1000
#define NUM_DISCARD 4096
#define NUM_KEEP 18
#define ORB 17
#include <limits>

using namespace std;
complex<double> J(0,1);
float xmin,xmax,ymin,ymax;

complex<double> next(complex<double> z,double u)
{
	extern complex<double>J;
	complex<double> zn;
	return z*exp(J*(complex<double> (0.4-6/(pow(abs(z),2)+1))))*u+(complex<double>)1;
}

void plot_attractor(double u, int ninit)
{

	extern float xmin,xmax,ymin,ymax;
	srand(time(NULL));
	for (int nit=0;nit<ninit;nit++)				//do for 20000 diff init. point.
	{

//	int rn=rand(); complex<double> z(xmin+rn*(xmax-xmin)/RAND_MAX,ymin+rn*(ymax-ymin)/RAND_MAX);	//diagonal random pts, a la SoB
	complex<double> z(xmin+((float)rand()/RAND_MAX)*(xmax-xmin),xmin+((float)rand()/RAND_MAX)*(xmax-xmin));
	
	complex<double> zorig=z;	
	complex<double> orb[ORB];

	for (int ind=0;ind<ORB;ind++){orb[ind]=ind;}

	for (int i=0;i<NUM_DISCARD;i++){z=next(z,u);}	//discard 1st few iter

	bool not_periodic=1;
	for (int i=0;(i<NUM_KEEP)&&(not_periodic);i++)
	{
		z=next(z,u);
		if (i%(ORB)==0)
		{
			if (abs(orb[0]-orb[16])<std::numeric_limits<double>::epsilon()*1000)
			{
				cerr<<real(zorig)<<'\t'<<imag(zorig)<<endl;
				not_periodic=0;
				break;
			}

			for (int ind=0;ind<ORB;ind++){orb[ind]=rand();}
		}

		orb[i%(ORB)]=z;
	}
	if (not_periodic){cout<<real(zorig)<<'\t'<<imag(zorig)<<endl;}
	}
}


main(int argc, char **argv)
{
	extern float xmin,xmax,ymin,ymax;
	float u=atof(argv[1]);
	xmin=atof(argv[2]);
	xmax=atof(argv[3]);
	ymin=atof(argv[4]);
	ymax=atof(argv[5]);
	int ninit=atoi(argv[6]);
	plot_attractor(u,ninit);
}
