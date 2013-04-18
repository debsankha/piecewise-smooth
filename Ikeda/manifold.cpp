#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<complex>
#include<vector>
#include<algorithm>
#define NTRY 100
#define NUM_PTS 2000
#define NUM_KEEP 50
#define MAXDIST 30

using namespace std;
complex<double> J(0,1);
complex<double> v1;
complex<double> v2;
double  u;

complex<double> nextx(complex<double> z,int n);
complex<double> prevx(complex<double> z,int n);
complex<double> nraph(complex<double> X, int n);
double f(complex<double> X,int n);	//|x(m+n)-x(m)|**2;
complex<double> numgrad(complex<double> X, int n);
bool issaddle(complex<double> X, int n);
double myrand(double xmin, double xmax);

void manifold(complex<double>x);

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
	int n=1;
	double tol=10e-8;
	double tolsa=10e-4;
//	found=[];
	
	srand(time(NULL));	
	double dx;
	complex<double> X, Xn;
	
	bool notin=1;
	int howmany=0;
	double df;
	complex<double> fp[NTRY];
	
	int nstep,maxstep;
	maxstep=4000;

	for (int ntry=0; ntry<NTRY;ntry++)
	{
		double dx=1;
//		complex<double> X(myrand(xmin,xmax),myrand(ymin,ymax));
		int rn=rand(); complex<double> X(xmin+rn*(xmax-xmin)/RAND_MAX,ymin+rn*(ymax-ymin)/RAND_MAX);	//diagonal random pts, a la SoB
		complex<double> Xn;
		
		nstep=0;
		while ((dx>tol) && (nstep<maxstep))
		{
			X=nraph(X,n);
			Xn=nextx(X,n);
			dx=abs(Xn-X);
			nstep++;
		}
	
		if (nstep<maxstep)
		{	
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
	}
	
	
	vector<complex<double> > vecarr(fp,fp+howmany);

	sort(vecarr.begin(),vecarr.begin()+howmany,comp);

	for (int i=0;i<howmany;i++)
	{
		X=vecarr[i];
		if (issaddle(X,1))
		{
			cout<<real(X)<<'\t'<<imag(X)<<endl;
			manifold(X);
			break;
		}
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

complex<double> prevx(complex<double> z,int n)
{
	extern complex<double> J;
	extern double u;
	complex<double> zmd;
	for (int i=0;i<n;i++)
	{
		zmd=(z-1.0)/u;
		z=zmd*exp(-J*(complex<double> (0.4-6/(pow(abs(zmd),2)+1))));
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


bool issaddle(complex<double> X, int n)
{
	extern complex<double> v1,v2;
	complex<double> gx,gy;
	double h=10e-4;
	gx=(-nextx(X+complex<double>(2*h,0),n)+nextx(X+complex<double>(h,0),n)*8.0-nextx(X-complex<double>(h,0),n)*8.0+nextx(X-complex<double>(2*h,0),n))/(12*h);
	gy=(-nextx(X+complex<double>(0,2*h),n)+nextx(X+complex<double>(0,h),n)*8.0-nextx(X-complex<double>(0,h),n)*8.0+nextx(X-complex<double>(0,2*h),n))/(12*h);
	double a=real(gx);
	double b=real(gy);
	double c=imag(gx);
	double d=imag(gy);
	
	complex<double> D=pow(complex<double>(pow(a-d,2)+b*c*4.0),0.5);
	
	if (imag(D)==0) 
	{
		double l1=real(D+a+d)/2;
		double l2=real(-D+a+d)/2;

				
		if (abs(l1)>abs(l2)){double temp=l1;l1=l2;l2=temp;}
		if (abs(l1)<1 && abs(l2)>1)
		{

			v1=(b,(l1-a));		//stable manif. starts here
			v2=(b,(l2-a));

			v1=v1/abs(v1);
			v2=v2/abs(v2);
			return 1;
		}
		else return 0;
	}
	else {return 0;}
}


void manifold(complex<double>x)
{
	extern complex<double> v1,v2;
	complex<double> xn;
	int num=1;
	int ind=0;
	double dx=10e-6/NUM_PTS;
	for (num=-NUM_PTS;num<NUM_PTS;num++)
	{
		xn=x+v1*(num*dx);
		
		for (ind=0;(ind<NUM_KEEP) && (abs(xn)<MAXDIST) ;ind++)
		{
			xn=prevx(xn,1);
			cerr<<real(xn)<<'\t'<<imag(xn)<<endl;
		}

	}

	for (num=-NUM_PTS;num<NUM_PTS;num++)
	{
		xn=x+v2*(num*dx);	//unstable manif
		
		for (ind=0;(ind<NUM_KEEP) && (abs(xn)<MAXDIST);ind++)
		{
			xn=nextx(xn,1);
			cout<<real(xn)<<'\t'<<imag(xn)<<endl;
		}
	}
}



double myrand(double xmin, double xmax)
{
	return xmin+((double)rand()/RAND_MAX)*(xmax-xmin);
}

