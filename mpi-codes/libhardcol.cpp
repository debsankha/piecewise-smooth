#include <hardcol.h>
#define NPTS 10		//# of pts to take for each param velue in bifurc diagram
#define EPSILON 0.001
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

float Sigma=1;		//the boundary: x=Sigma
float G=0.08;	//damping
float F=0.393094; //1.4881;		//forcing amplitude
float W=1;	//forcing freq: sin(W*t)	NOTE: w_0=1
float m=2.3556;


void f(double t, vec x, vec *out) 
{
	//SHM hitting with a hard wall , x=(x,v)
	float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
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

int plotpoincare(vec x, double tmin, double tmax, double t_startprint)
{
	double t=tmin;
	double T=2*M_PI/W;
	double oldvel=0;
	int i=0;	
	double poinc_x[ORB], poinc_t[ORB], time_to_stable;
	int period=0;
		
	t=0;
	while (t<t_startprint)
	{
		cerr<<t<<'\t'<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
		oldvel=x.arr[1];
		rk4(t,&x);
		if (x.arr[0]>Sigma)		//The reset map on hard collision
		{
			x.arr[0]=Sigma;
			x.arr[1]*=-1;
		}


		if ((oldvel<0) && (x.arr[1]>0)) //(fmod(t,T)<h)  DOES NOT WORK well due to non-zero time step
		{
			poinc_x[i]=x.arr[0];
			i++;
			if (i==ORB)
			{
				i=0;
				period=detect_period(&poinc_x[0],&poinc_t[0],&time_to_stable);
				
				if (period>0)
				{

					cout<<"#period: "<<period<<endl;
					for (int j=ORB;j>0;j--) cout<<F<<'\t'<<poinc_x[ORB-j]<<endl;
					return period;
				}
			}
		}
	
		t+=h;
	}//discarding the transient
	
	i=0;
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
		if ((oldvel<0) && (x.arr[1]>0)) //(fmod(t,T)<h)  DOES NOT WORK well due to non-zero time step
		{
			poinc_x[i]=x.arr[0];
			i++;
			if (i==ORB)
			{
				i=0;
				period=detect_period(&poinc_x[0],&poinc_t[0],&time_to_stable);
				
				if (period>0)
				{
					cerr<<"#period: "<<period<<endl;
					return period;
				}
			}
		}
		t+=h;
		cout<<F<<'\t'<<x.arr[0]<<endl;			//responsible for the output
	}
	return 0;
}

int plotbifurc_F(float minF, float maxF, int npts)
{
	double dF=(maxF-minF)/npts;
	double  t=0;
	vec x(N);
	float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
	for (F=minF;F<maxF;F+=dF)
	{
		for (int cnt=0;cnt<NPTS;cnt++)
		{
			x.arr[0]=-Sigma+randdouble(-0.01,Sigma);	
			x.arr[1]=1.1*pow(K1*(Sigma+x.arr[0])*(Sigma-x.arr[0]),0.5);	
			cout<<"#Starting from: "<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
			plotpoincare(x,t,TMAX,TMAX*0.9);
			cout<<endl;				//Crucial for block detection by gnuplot
		}
	}
}

int plotbasin(int npts,double tmax, int rank)
{
	int i,period;
	char fname[20],command[30];
	sprintf(fname, "basindat%d.dat",rank);
	vec x(2);

	float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
	double initphase=(-M_PI/2-atan(G*W/(W*W-K1)))/W;

	FILE *outf;

	for (i=0;i<npts;i++)
	{
		x.arr[0]=randdouble(-8,1);
		x.arr[1]=randdouble(-8,8);
		
		//dup cerr to some file and cout to /dev/null
		outf=freopen(fname, "w", stderr );
		period=plotpoincare(x,initphase,tmax, tmax+1);
		
		sprintf(command,"cat %s >>period%d.dat",fname,period);
		system(command);

		fclose(outf);
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
			cerr<<"#Period:"<<per<<endl;
//			for (int i=0; i<per; i++) cerr<<F<<'\t'<<G<<'\t'<<arr[ORB-i-1]<<endl;//get a bifurcation giadram for free
			return per;
		}
	}
	*time_to_stable=-1;	//Just some absurd value: should never be accepted by the calling function
	return 0;
}
