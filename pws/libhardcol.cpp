#include <hardcol.h>
#define NPTS 10		//# of pts to take for each param velue in bifurc diagram
#define EPSILON 0.00000001
using namespace std;

float Sigma=1;		//the boundary: x=Sigma
float G=0.062;	//damping
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
	cerr<<"G="<<G<<endl;	
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



int plotpoincare(vec x, double tmin, double tmax)
{
	double t=tmin;
	double T=2*M_PI/W;
	double oldvel=0;
	int i=0;	
	double poinc_x[ORB], poinc_t[ORB], time_to_stable;
	int period=0;
	
	i=0;
	
	double oldx=x.arr[0];

	while (t<tmax)
	{
		cerr<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
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
					cerr<<"#period: "<<period<<endl;
					break;
				}
			}

			cout<<oldx<<'\t'<<x.arr[0]<<endl;
			cout<<x.arr[0]<<'\t'<<x.arr[0]<<endl;		//for cobweb
			oldx=x.arr[0];

		}
		t+=h;
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
			plotpoincare(x,t,TMAX);
			cout<<endl;				//Crucial for block detection by gnuplot
		}
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
			cout<<"#stretch: "<<stretch<<endl;
			for (int i=0; i<per; i++) cerr<<F<<'\t'<<G<<'\t'<<arr[ORB-i-1]<<endl;//get a bifurcation giadram for free
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
