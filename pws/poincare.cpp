#define _USE_MATH_DEFINES
#include <cmath>
#include <complex>
#include <rk4.h>
#define NPTS 200	//number of x values to take for bifurc diag
#define NPTSEACHX 10
#define NMAX 300
#define randdouble(min,max) min+rand()*(max-min)*1.0/RAND_MAX
#define ORB 100
#define energy(x,y) y*y+K1*x*x

using namespace std;

float Sigma=1;		//the boundary: x=Sigma
float G=0.062;	//damping
float F=4.43741; //1.4881;		//forcing amplitude
float W=1;	//forcing freq: sin(W*t)	NOTE: w_0=1
float m=2.336;
float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
double Tau=M_PI/W;	//time when stable orbit grazes
double Amp;

int detect_period(complex<double> *arr);
void mtintox(double t,vec, vec *);
void poinc_x(vec *);

void f(double, vec, vec*){;}

int main(int argc, char **argv)
{
	cerr<<"K1: "<<K1<<endl;
	cerr<<"#delta_t\txstart\txhcol\tvhcol\tvpcol\n";
	srand(0);
	float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
	cout<<"#F_graz: "<<F_graz<<endl;
	float F_range=atof(argv[1]);
	float dF=F_range/NPTS;
	float Fmax=F_graz+F_range/2;
	float Fmin=F_graz-F_range/2;
	int period,startprint;
	
	complex<double> arr[ORB];
		
	vec x(2);
	for (F=Fmin;F<Fmax;F+=dF)
	{
		for (int n=0;n<NPTSEACHX;n++)
		{
			//poincare all the way
			x.arr[0]=randdouble(-0.01,0.01);
			x.arr[1]=randdouble(-0.01,0.01);

			int n_iter=0;

			while (n_iter<NMAX)
			{
				arr[n_iter%ORB]=complex<double>(x.arr[0],x.arr[1]);
				poinc_x(&x);
				n_iter++;

				if ((n_iter%ORB)==0)
				{
					period=detect_period(&arr[0]);
					if (period!=0) break;
				}
			}
			
			cout<<"#period: "<<period<<endl;	
			if (period==0) startprint=0;
			else startprint=ORB-period;
			for (int i=startprint;i<ORB;i++) cout<<F<<'\t'<<-Amp+real(arr[i])<<'\t'<<imag(arr[i])<<endl;
		}
		cout<<endl;
	}
	
}

void poinc_x(vec *x)
{
	vec rhcol(2);
	mtintox(Tau,*x, &rhcol);
	double xh=rhcol.arr[0];
	double vh=rhcol.arr[1];

//	cout<<"xh: "<<xh<<" vh: "<<vh<<endl;

	Amp=F/pow((K1-W*W)*(K1-W*W)+W*W*G*G,0.5);	//amplitude of stable orb

	double delta_t;
	
	double D=vh*vh+2*Amp*W*W*(xh+Amp-Sigma);

	if (D>0)
	{
		if (vh>0) delta_t=(-vh+pow(D,0.5))/(-Amp*W*W);
		else delta_t=(-vh-pow(D,0.5))/(-Amp*W*W);

		double vhcol=vh*(1-G*delta_t)-K1*delta_t*xh;
		double vpcol=-Amp*W*sin(W*delta_t);
		
//		cout<<'#'<<delta_t<<'\t'<<x->arr[0]<<'\t'<<xh<<'\t'<<vhcol<<'\t'<<vpcol<<endl;	
		
		vec extra_x(2);
		extra_x.arr[0]=0;
		extra_x.arr[1]=-2*(vhcol+vpcol);
		
		vec firstterm(2),secondterm(2);
	
		mtintox(Tau-delta_t,extra_x,&secondterm);	//because in our case, T=2*Tau
		mtintox(2*M_PI/W,*x,&firstterm);
		*x=firstterm+secondterm;
	}
	else
	{
		vec firstterm(2);
		mtintox(2*M_PI/W,*x,&firstterm);
		*x=firstterm;
	}
	
//	cout<<delta_t<<endl;	

}


void mtintox(double t,vec x, vec *y)
{
	double wg=pow(4*K1-G*G,0.5)/2;
	double damping=exp(-G*t/2)/wg;
	double coswgt=cos(wg*t);
	double sinwgt=sin(wg*t);

	y->arr[0]=damping*((wg*coswgt+G*sinwgt/2)*x.arr[0]+sinwgt*x.arr[1]);
	y->arr[1]=damping*(-K1*sinwgt*x.arr[0]+(wg*coswgt-G*sinwgt/2)*x.arr[1]);
}





int detect_period(complex<double> *arr)
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
			if (inc<0.001)
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
			return per;
		}
	}
	return 0;
}
