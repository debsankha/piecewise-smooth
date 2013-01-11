#include<iostream>
#include<cmath>
#include<stdlib.h>
#define sigma 5

using namespace std;

main(int argc, char **argv)
{
	double x0,v0,x,v;
	x0=0;
	v0=12;
	x=atof(argv[1]);
	v=atof(argv[2]);

	double T=0.56;
	double tcol=0.493;
	double t=0;
	double w1=2;
	double w2=pow(8,0.5);
	double vcol=6.61768;
	
	double dx=x-x0;
	double dv=v-v0;
	
	dx=cos(w1*tcol)*dx+sin(w1*tcol)*dv/w1;	
	dv=-w1*sin(w1*tcol)*dx+cos(w1*tcol)*dv;
	
	dv+=-sigma*dx*w1*w1/6.61768;

	dx=cos(w2*(T-tcol))*dx+sin(w2*(T-tcol))*dv/w2;	
	dv=-w2*sin(w2*(T-tcol))*dx+cos(w2*(T-tcol))*dv;
	
	cout<<dx<<'\t'<<dv<<endl;	
}

