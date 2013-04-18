#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<complex>
#define NUM_ITER 1000
#define NUM_DISCARD 10000
#define NUM_KEEP 16384
#define ORB 49
#include <limits>

using namespace std;
complex<double> J(0,1);

complex<double> next(complex<double> z,double u)
{
	complex<double> zn;
	return z*exp(J*(complex<double> (0.4-6/(pow(abs(z),2)+1))))*u+(complex<double>)1;
}

void plot_attractor(double u, complex<double> z)
{

	complex<double> orb[ORB];

	for (int ind=0;ind<ORB;ind++){orb[ind]=ind;}

	for (int i=0;i<NUM_DISCARD;i++){z=next(z,u);}	//discard 1st few iter
	
	bool not_periodic=1;
	for (int i=0;(i<NUM_KEEP)&&(not_periodic);i++)
	{
		z=next(z,u);
		if (i%(ORB)==0)
		{
			for (int ind=1;ind<ORB;ind++)
			{
				if (abs(orb[0]-orb[ind])<std::numeric_limits<double>::epsilon()*1000)
				{
					cerr<<"Z:"<<orb[0]<<'\t'<<ind<<endl;
					not_periodic=0;
					break;
				}
			}

			for (int ind=0;ind<ORB;ind++){orb[ind]=rand();}
		}

		orb[i%(ORB)]=z;
		cout<<real(z)<<'\t'<<imag(z)<<endl;
	}
}


main(int argc, char **argv)
{
	extern float xmin,xmax,ymin,ymax;
	float u=atof(argv[1]);
	float x=atof(argv[2]);
	float y=atof(argv[3]);
	plot_attractor(u,complex<double>(x,y));
}
