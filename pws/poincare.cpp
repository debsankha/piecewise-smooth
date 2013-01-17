#define _USE_MATH_DEFINES
#include <cmath>
#include <rk4.h>
#define NPTS 100	//number of x values to take for bifurc diag
#define NPTSEACHX 1
#define NMAX 200
#define randdouble(min,max) min+rand()*(max-min)*1.0/RAND_MAX
#define ORB 50

using namespace std;

float Sigma=1;		//the boundary: x=Sigma
float G=0.062;	//damping
float F=4.43741; //1.4881;		//forcing amplitude
float W=1;	//forcing freq: sin(W*t)	NOTE: w_0=1
float m=2.3556;
float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
double Tau=M_PI/W;	//time when stable orbit grazes

int detect_period(double *arr, int *n_iter);
void mtintox(double t,vec, vec *);
void poinc_x(vec *);

void f(double, vec, vec*){;}

int main(int argc, char **argv)
{
	srand(0);
	float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
	float F_range=atof(argv[1]);
	float dF=F_range/NPTS;
	float Fmax=F_graz+F_range/2;
	float Fmin=F_graz-F_range/2;
	
	cerr<<"dt"<<"xhcol"<<"vhcol"<<"vpcol"<<endl;

	for (F=Fmin;F<Fmax;F+=dF)
	{
		for (int n=0;n<NPTSEACHX;n++)
		{
			//poincare all the way
			double tmp[]={randdouble(-0.001,0.001),randdouble(-0.001,0.001)};
			vec x(2,tmp);
			int n_iter=0;

			while (n_iter<NMAX)
			{
				x.show();
				poinc_x(&x);
				n_iter++;
			}

		}
	}


}


void poinc_x(vec *x)
{
	vec xhcol(2);
	mtintox(Tau,*x, &xhcol);
	
	double delta_t=-1*xhcol.arr[0]/xhcol.arr[1];
	
	double vhcol=xhcol.arr[1]*(1-G*delta_t)-K1*delta_t*xhcol.arr[0];
	double vpcol=-F*W/pow((K1-W*W)*(K1-W*W)+W*W*G*G,0.5)*sin(delta_t);
	
	vec extra_x(2);
	extra_x.arr[0]=0;
	extra_x.arr[1]=-2*(vhcol+vhcol);
	
	vec firstterm(2),secondterm(2);

	mtintox(Tau-delta_t,extra_x,&secondterm);	//because in our case, T=2*Tau

	mtintox(2*M_PI/W,*x,&firstterm);
	*x=firstterm+secondterm;
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





int detect_period(double *arr, int *n_iter)
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
			*n_iter-=stretch;
			cerr<<"stretch: "<<stretch<<endl;
			return per;
		}
	}
	return 0;
}
