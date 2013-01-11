#include <rk4.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <string.h>
#include <limits>
#include <stdlib.h>
#define TMAX 500	//Time to evolve for bifurc diagrams
#define NPTS 1		//# of pts to take for each param velue in bifurc diagram
#define ORB  50		//set it to 2+max periodicity that we'll attempt to detect.(extra 1 for redundancy)
#define N 2		//dimensionality of system
#define randdouble(min,max) min+rand()*(max-min)*1.0/RAND_MAX


float Sigma=1;		//the boundary: x=Sigma
float G=0.062;	//damping
float F=4.43741; //1.4881;		//forcing amplitude
float W=1;	//forcing freq: sin(W*t)	NOTE: w_0=1
float m=2.3556;
float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;

using namespace std;
int plottraj(vec, float);
int plotpoincare(vec, double, double);
int plotbifurc_F(float minF, float maxF, int npts);
int detect_period(double *);

void f(double t, vec x, vec *out) 
{
	//SHM hitting with a hard wall , x=(x,v)
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0]-G*x.arr[1]+F*sin(W*t);
}


int main(int argc , char *argv[])
{
	srand(0);
	if (strcmp(argv[1],"traj")==0)
	{
		int i;	
		double tmax,tmp[N];
	
		for (i=1;i<N+1;i++) { tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
		vec x(N,tmp);	
		tmax=atof(argv[1+i]);
		F=atof(argv[2+i]);
		plottraj(x,tmax);
	}
	if (strcmp(argv[1],"poincare")==0)
	{
		int i;	
		double tmax,tmp[N];
	
		for (i=1;i<N+1;i++) { tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
		vec x(N,tmp);
	
		tmax=atof(argv[1+i]);
		F=atof(argv[2+i]);
		plotpoincare(x,0,tmax);
	}

	if (strcmp(argv[1],"plotbifurc")==0)
	{
		float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
		cerr<<"F_graz="<<F_graz<<endl;
		plotbifurc_F(F_graz-atof(argv[2]),F_graz+atof(argv[3]),atoi(argv[4]));
	}

}

int plottraj(vec x, float tmax)
{
	double t;
	double oldvel=0;
	double poinc_x[ORB];							//array to store poincare x vals in to detect period
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
			i++;
		} 
		
		if (i==ORB)
		{
			period=detect_period(&poinc_x[0]);
			
			if (period>0)
			{
				cerr<<"# Period: "<<period<<endl;
				cerr<<"time to stabilize: "<<t<<endl;
				return 0;
			}
			i=0;
		}
		
		t+=h;
	}
	return 0;
}

int plotpoincare(vec x, double tmin, double tmax)
{
	double t=tmin;
	double T=2*M_PI/W;
	double oldvel=0;
	while (t<tmax)
	{
		oldvel=x.arr[1];
		rk4(t,&x);
		if (x.arr[0]>Sigma)		//The reset map on hard collision
		{
			x.arr[0]=Sigma;
			x.arr[1]*=-1;
		}
		t+=h;
		if ((oldvel<0) && (x.arr[1]>0)) //(fmod(t,T)<h)  DOES NOT WORK well due to non-zero time step
		{
			cout<<F<<'\t'<<t<<'\t'<<x.arr[0]<<endl;			//responsible for the output
		}
	}
	return 0;

}

int plotbifurc_F(float minF, float maxF, int npts)
{
	double dF=(maxF-minF)/npts;
	double  t=0;
	vec x(N);
	for (F=minF;F<maxF;F+=dF)
	{
		for (int cnt=0;cnt<NPTS;cnt++)
		{
			x.arr[0]=-Sigma+randdouble(-0.01,Sigma);	
			x.arr[1]=1.1*pow(K1*(Sigma+x.arr[0])*(Sigma-x.arr[0]),0.5);	
			cout<<"#Starting from: "<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
			t=0;
			while (t<TMAX*0.9)
			{
				rk4(t,&x);
				if (x.arr[0]>Sigma)		//The reset map on hard collision
				{
					x.arr[0]=Sigma;
					x.arr[1]*=-1;
				}
				t+=h;
			}//discarding the transient
		
			plotpoincare(x,t,TMAX);
			cout<<endl;				//Crucial for block detection by gnuplot
		}
	}
}


int detect_period(double *arr)
{
	for (int i=1;i<ORB-1;i++)
	{
		if (abs(arr[i]-arr[0])<0.001)
		{
			if (abs(arr[i+1]-arr[1])<0.001) return i;
		}
	}
	return 0;
}
