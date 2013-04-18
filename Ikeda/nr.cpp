#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<complex>
#include<vector>
#include<algorithm>
#define NTRY 100

using namespace std;
complex<double> J(0,1);
double  u;

complex<double> nextx(complex<double> z,int n);
complex<double> nraph(complex<double> X, int n);
double f(complex<double> X,int n);	//|x(m+n)-x(m)|**2;
complex<double> numgrad(complex<double> X, int n);
void eig(complex<double> X, int n);
double myrand(double xmin, double xmax);

bool comp(complex<double> x, complex<double> y)
{
	if (real(x)>real(y)){return 1;}
	else if (real(x)<real(y)){return 0;}
	else
	{
		if (imag(x)>imag(y)){return 1;}
		else {return 0;}
	}
}



main(int argc,char **argv)
{
	extern double u;
	u=atof(argv[1]);
	float xmin=atof(argv[2]);
	float xmax=atof(argv[3]);
	float ymin=atof(argv[4]);
	float ymax=atof(argv[5]);
	float n=atoi(argv[6]);

	double tol=10e-8;
	double tolsa=10e-5;
//	found=[];
	
	srand(time(NULL));	
	double dx;
	complex<double> X, Xn;
	
	bool notin=1;
	int howmany=0;
	double df;
	complex<double> fp[NTRY];

	for (int ntry=0; ntry<NTRY;ntry++)
	{
		double dx=1;
		complex<double> X(myrand(xmin,xmax),myrand(ymin,ymax));
		complex<double> Xn;
		
		while (dx>tol)
		{
			X=nraph(X,n);
			Xn=nextx(X,n);
			dx=abs(Xn-X);
		}

		notin=1;
		for (int i=0;i<howmany;i++)
		{
			df=abs(fp[i]-X);
			
			if (df<tolsa)
			{
				notin=0;
				break;
			}
		}
			
		if (notin)
		{
			fp[howmany]=X;
			howmany++;
		}
	}
	
	
	vector<complex<double> > vecarr(fp,fp+howmany);

	sort(vecarr.begin(),vecarr.begin()+howmany,comp);

	for (int i=0;i<howmany;i++)
	{
		X=vecarr[i];
		cout<<X<<'\t'<<nextx(X,n)<<'\t';
		eig(X,n);
	}
}


complex<double> nextx(complex<double> z,int n)
{
	extern complex<double> J;
	extern double u;
	for (int i=0;i<n;i++)
	{
		z=z*exp(J*(complex<double> (0.4-6/(pow(abs(z),2)+1))))*u+(complex<double>)1;
	}
	return z;
}

complex<double> nraph(complex<double> X, int n)
{
	complex<double> gr=numgrad(X,n);
	complex<double> alphgr=(-f(X,n)/(pow(abs(gr),2)))*gr;
	return X+alphgr;
}

double f(complex<double> X,int n)	//|x(m+n)-x(m)|**2
{
	complex<double> Y=nextx(X,n);
	complex<double> diff=Y-X;
	return pow(abs(diff),2);
}

complex<double> numgrad(complex<double> X, int n)
{
	double gx,gy;
	double h=10e-4;
	gx=(-f(X+complex<double>(2*h,0),n)+f(X+complex<double>(h,0),n)*8.0-f(X-complex<double>(h,0),n)*8.0+f(X-complex<double>(2*h,0),n))/(12*h);
	gy=(-f(X+complex<double>(0,2*h),n)+8*f(X+complex<double>(0,h),n)-8*f(X-complex<double>(0,h),n)+f(X-complex<double>(0,2*h),n))/(12*h);
	
	return complex<double>(gx,gy);
}

void eig(complex<double> X, int n)
{
	complex<double> gx,gy;
	double h=10e-4;
	gx=(-nextx(X+complex<double>(2*h,0),n)+nextx(X+complex<double>(h,0),n)*8.0-nextx(X-complex<double>(h,0),n)*8.0+nextx(X-complex<double>(2*h,0),n))/(12*h);
	gy=(-nextx(X+complex<double>(0,2*h),n)+nextx(X+complex<double>(0,h),n)*8.0-nextx(X-complex<double>(0,h),n)*8.0+nextx(X-complex<double>(0,2*h),n))/(12*h);
	double a=real(gx);
	double c=real(gy);
	double b=imag(gx);
	double d=imag(gy);
	
	complex<double> D=pow(complex<double>(pow(a-d,2)+b*c*4.0),0.5);
	
	if (imag(D)==0) {cout<<real(D+a+d)/2.0<<'\t'<<real(-D+a+d)/2.0<<'\n';}
	else {cout<<(D+a+d)/2.0<<'\t'<<(-D+a+d)/2.0<<'\n';}
}


double myrand(double xmin, double xmax)
{
	return xmin+((double)rand()/RAND_MAX)*(xmax-xmin);
}


