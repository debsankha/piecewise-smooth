#include <hardcol.h>
#define NPTS 10		//# of pts to take for each param velue in bifurc diagram
#define EPSILON 0.001

#include <stdlib.h>
#include <stdio.h>
#include <complex>
using namespace std;

float Sigma=1;		//the boundary: x=Sigma
float G=0.08;	//damping
float F=0.393094; //1.4881;		//forcing amplitude
float W=1;	//forcing freq: sin(W*t)	NOTE: w_0=1
float m=2.3556;
float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;


void f(double t, vec x, vec *out) 
{
	//SHM hitting with a hard wall , x=(x,v)
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0]-G*x.arr[1]+F*sin(W*t);
}

int plottraj(vec x, float tmax)
{
	double t;
	double time_to_stable=0;
	double oldvel=0;
	double poinc_x[ORB];							//array to store poincare x vals in to detect period
	double poinc_t[ORB];
	int i=0;
	int period=0;
			
	t=0;
	cerr<<"F="<<F<<endl;	
	while (t<tmax)
	{
		oldvel=x.arr[1];
		cout<<t<<'\t';x.show();						//responsible for the output to stdout
		rk4(t,&x);
		if (x.arr[0]>Sigma)						//The reset map on hard collision
		{
			x.arr[0]=Sigma;
			x.arr[1]*=-1;
		}

		if ((oldvel<0) && (x.arr[1]>0))
		{
			cerr<<t<<'\t'<<x.arr[0]<<endl;				//outputs the poincare mappings to stderr
			poinc_x[i]=x.arr[0];
			poinc_t[i]=t;
			i++;
		 
			if (i==ORB)
			{
				i=0;
				period=detect_period(&poinc_x[0],&poinc_t[0],&time_to_stable);
				
				if (period>0)
				{
					cerr<<"# Period: "<<period<<endl;
					cerr<<"time to stabilize: "<<time_to_stable<<endl;
					return period;
				}
			}
		}	
		t+=h;
	}
	return 0;
}


int plotmap(int npts, int n, float newF, float newG)
{
	F=newF;
	G=newG;

	double t;
	double T=2*M_PI/W;
	double oldvel=0;
	vec x(2);
	
	double xini;
	//Suposing interesting things happen only in th interval (-8,1)
	double xinincr=9*1.0/npts;
	
	int numpoint;
	double initphase;
	double Amp=-F/pow(pow(W*W-K1,2)+W*W*G*G,0.5);

	cerr<<"#G_graz\t"<<pow(F*F/(Sigma*Sigma)-(K1-W*W)*(K1-W*W),0.5)/W<<'\t'<<"F_graz\t"<<Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5)<<'\t'<<"Amp\t"<< Amp <<endl;

	for (xini=-8;xini<1;xini+=xinincr)
	{
		numpoint=0;
		initphase=(-M_PI/2-atan(G*W/(W*W-K1)))/W;//Doesn't really matter in the new way of getting the map: rely on suucessive stroboscopic points and map strob_pt[0] to strobe_pt[1]
		t=initphase;
		x.arr[0]=xini;
		x.arr[1]=0;

		cout<<x.arr[0]<<'\t';

		while (t<n*T+initphase)
		{
			oldvel=x.arr[1];
			rk4(t,&x);
			if (x.arr[0]>Sigma)		//The reset map on hard collision
			{
				x.arr[0]=Sigma;
				x.arr[1]*=-1;
			}
	
			t+=h;

		}
		cout<<x.arr[0]<<endl;
		if (abs(xini-Amp)<0.25) 
		{
			xinincr=9*0.5/npts;		//look more closely near the f.p.
		}
		else xinincr=9*1.0/npts;	
	}
	return 0;
}


int plotmap3d(int npts, int n, float newF, float newG)
{
	F=newF;
	G=newG;

	double t;
	double T=2*M_PI/W;
	double oldvel=0;
	vec x(2);
	
	double xini,vini;
	//Suposing interesting things happen only in th interval (-8,1)
	double xinincr=9*1.0/npts;
	
	int numpoint;
	double initphase;
	double Amp=-F/pow(pow(W*W-K1,2)+W*W*G*G,0.5);

	cerr<<"#G_graz\t"<<pow(F*F/(Sigma*Sigma)-(K1-W*W)*(K1-W*W),0.5)/W<<'\t'<<"F_graz\t"<<Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5)<<'\t'<<"Amp\t"<< Amp <<endl;
	

	for (vini=-8;vini<8;vini+=xinincr)
	{	
	for (xini=-8;xini<1;xini+=xinincr)
	{
		numpoint=0;
		initphase=(-M_PI/2-atan(G*W/(W*W-K1)))/W;//Doesn't really matter in the new way of getting the map: rely on suucessive stroboscopic points and map strob_pt[0] to strobe_pt[1]
		t=initphase;
		x.arr[0]=xini;
		x.arr[1]=vini;

		cout<<x.arr[0]<<'\t'<<x.arr[1]<<'\t';

		while (t<n*T+initphase)
		{
			oldvel=x.arr[1];
			rk4(t,&x);
			if (x.arr[0]>Sigma)		//The reset map on hard collision
			{
				x.arr[0]=Sigma;
				x.arr[1]*=-1;
			}
	
			t+=h;

		}
		cout<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
//		if (abs(xini-Amp)<0.25) 
//		{
//			xinincr=9*0.5/npts;		//look more closely near the f.p.
//		}
//		else xinincr=9*1.0/npts;	
	}
	}
	return 0;
}



int plotpoincare(vec x, double tmin, double tmax)
{
	double t=tmin;
	double T=2*M_PI/W;
	double oldvel=0;
	int i=0;	
	int period=0;
	
	i=0;
	
	double oldx=x.arr[0];

	double poinc_x[ORB], poinc_t[ORB], time_to_stable;
	while (t<tmax)
	{
		cerr<<t<<'\t'<<x.arr[0]<<'\t'<<x.arr[1]<<endl;

		oldvel=x.arr[1];
		rk4(t,&x);
		if (x.arr[0]>Sigma)		//The reset map on hard collision
		{
			x.arr[0]=Sigma;
			x.arr[1]*=-1;
		}
		if (fmod(t-tmin,T)<h) //((oldvel<0) && (x.arr[1]>0)) 
		{
			poinc_x[i]=x.arr[0];
			i++;
			if (i==ORB)
			{
				i=0;
				period=detect_period(&poinc_x[0],&poinc_t[0],&time_to_stable);
				
				if (period>0)
				{
					cerr<<"#time to stabilize: "<<time_to_stable<<endl;
					cerr<<"#period: "<<period<<endl;
					return period;
				}
			}

			cout<<oldx<<'\t'<<x.arr[0]<<endl;
			cout<<x.arr[0]<<'\t'<<x.arr[0]<<endl;		//for cobweb
			oldx=x.arr[0];

		}
		t+=h;
	}
	return -1;
}



int plotbasin(int npts,double tmax,double newf,double newg)
{
	F=newf;
	G=newg;

	float minx=-8;
	float maxx=-1;
	float miny=-8;
	float maxy=8;
	
	double T=2*M_PI/W;
	double initphase=(-M_PI/2-atan(G*W/(W*W-K1)))/W;
	double t,time_to_stable;

	int i,period,boxxind,boxyind;

	vec x(2);

	int numxdiv=npts;
	int numydiv=npts;

	short int basin[numydiv][numxdiv];
	bool tempbasin[numxdiv][numydiv];

	for (int xind=0;xind<numxdiv;xind++)
	{
		for (int yind=0;yind<numydiv;yind++)
		{
			basin[yind][xind]=2;
		}
	}
	
	double poinc_x[ORB], poinc_t[ORB];
	for (int xptco=0;xptco<npts+1;xptco++)
	{
	for (int yptco=0;yptco<npts+1;yptco++)
	{
	
		i=0;
		t=0;

		x.arr[0]=minx+(maxx-minx)*xptco*1.0/npts;
		x.arr[1]=miny+(maxy-miny)*yptco*1.0/npts;

//		cerr<<"starting from: "<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
		for (int xind=0;xind<numxdiv;xind++)
		{
			for (int yind=0;yind<numydiv;yind++)
			{
				tempbasin[yind][xind]=0;
			}
		}

		
		while (t<tmax)
		{
			boxxind=(x.arr[0]-minx)*numxdiv/(maxx-minx);
			boxyind=(x.arr[1]-miny)*numydiv/(maxy-miny);

			rk4(t,&x);
			if (x.arr[0]>Sigma)		//The reset map on hard collision
			{
				x.arr[0]=Sigma;
				x.arr[1]*=-1;
			}
			
			

			
			if ((boxxind>-1) && (boxxind<numxdiv) && (boxyind>-1) && (boxyind<numydiv))
			{
				tempbasin[boxyind][boxxind]=1;
				switch (basin[boxyind][boxxind])
				{
					case 1:
					{
						period=1;
						break;
					}
					case 0:
					{
						period=3;
						break;
					}
						
						
					case -1:
					{
						//point not on any basin drawn so far	
					}
				}
			}


			if (fmod(t-initphase,T)<h) //((oldvel<0) && (x.arr[1]>0)) 
			{
				poinc_x[i]=x.arr[0];
				i++;
				if (i==ORB)
				{
					i=0;
					period=detect_period(&poinc_x[0],&poinc_t[0],&time_to_stable);
					
					if (period>0)
					{
						break;
					}
				}
	
			}
			t+=h;
		}
		
		cerr<<"period: "<<period<<endl;

		switch (period)
		{
			case 1:
			{
				for (int xind=0;xind<numxdiv;xind++)
				{
					for (int yind=0;yind<numydiv;yind++)
					{
						if (tempbasin[yind][xind]==1) basin[yind][xind]=1;
					}
				}
				break;
			}

			default:
			{
				for (int xind=0;xind<numxdiv;xind++)
				{
					for (int yind=0;yind<numydiv;yind++)
					{
						if (tempbasin[yind][xind]==1) basin[yind][xind]=0;
					}
				}

				break;

			}

		}
	
	}
	}

	for (int yind=0;yind<numydiv;yind++)
	{
		for (int xind=0;xind<numxdiv;xind++)
		{
			if (basin[yind][xind]==1) cout<<"\\filled;& ";
			else if (basin[yind][xind]==0)  cout<<"\\unfilled;& ";
			else cout<<"\\undef;& ";
		}
		cout<<"\\\\\n";
	}
}


int detect_period(double *arr, double *t_arr, double *time_to_stable)
{
	double inc;
	int stretch;
	bool stillthere;

	for (int per=1;per<ORB/2;per++)	//per==periodicity to detect
	{
		stretch=0;
		stillthere=0;
		for (int startpt=0;startpt<ORB-per;startpt++)
		{
			inc=abs(arr[startpt]-arr[startpt+per]);
			if (inc<EPSILON)
			{
				stillthere=1;
				stretch++;
			}
			else
			{
				stillthere=0;
				stretch=0;		
			}
		}

		if (stillthere && (stretch>4))
		{
			*time_to_stable=t_arr[ORB-stretch-1];
			cerr<<"#stretch: "<<stretch<<endl;
//			for (int i=0; i<per; i++) cerr<<F<<'\t'<<G<<'\t'<<arr[ORB-i-1]<<endl;//get a bifurcation giadram for free
			return per;
		}
	}
	*time_to_stable=-1;	//Just some absurd value: should never be accepted by the calling function
	return 0;
}



double time_to_stabilize(vec x, double tmax)
{
	double t=0;
	double oldvel=0;
	double poinc_x[ORB];							//array to store poincare x vals in to detect period
	double poinc_t[ORB];							//array to store poincare x vals in to detect period
	int i=0;
	int period=0;
			
	t=0;
	double time_to_stable=0;
	while (t<tmax)
	{
		oldvel=x.arr[1];
		rk4(t,&x);
		if (x.arr[0]>Sigma)						//The reset map on hard collision
		{
			x.arr[0]=Sigma;
			x.arr[1]*=-1;
		}

		if ((oldvel<0) && (x.arr[1]>0))
		{
			poinc_x[i]=x.arr[0];
			poinc_t[i]=t;
			i++;
		
			if (i==ORB)
			{
				period=detect_period(&poinc_x[0],&poinc_t[0], &time_to_stable);
				
				if (period>0)
				{
					cerr<<"#F: "<<F<<"\tPeriod: "<<period<<endl;
					return time_to_stable;
				}
				i=0;
			}
		}
		t+=h;
	}
	return tmax;		//Remember: you'd get t=0 if chaotic orb or periodicity >48
}
