#include <rk4.h>
#include <stdlib.h>
#include <cmath>
#define N 4		//dimensionality of system
#define lmin 0		//The boundary of the box from which init points will be chosen randomly
#define lmax 10
#define	phimin 0 
#define	phimax 6.28

//x=(l,l',phi,phi')

using namespace std;
double mu;

void f(double t, vec x, vec *out) 
{
	//Pulley and swinging ball. See Tel & Gruiz Sec 7.4.3
	out->arr[0]=x.arr[1];
	out->arr[2]=x.arr[3];
	out->arr[1]=mu*(cos(x.arr[2])+x.arr[0]*pow(x.arr[3],2)+1)-1;
	out->arr[3]=-(sin(x.arr[2])+2*x.arr[1]*x.arr[3])/x.arr[0];
}


void plot_asymptote(vec *x)
{
	double t, tmp[N],tmax;
	tmax=5000;
	t=0;

	while ((t<tmax*0.8) && (x->norm()<100))
	{
		
		rk4(t,x);
		t+=h;
	}

	while ((t<tmax) && (x->norm()<100))
	{
		
		rk4(t,x);
		t+=h;
		cout<<mu<<'\t'<<t<<'\t'<<x->arr[0]<<'\t'<<x->arr[1]<<'\t'<<x->arr[2]<<'\t'<<x->arr[3]<<endl;
	}

}



int main(int argc , char *argv[])
{
	double mua=atof(argv[1]);
	double mub=atof(argv[2]);
	double dmu=atof(argv[3]);

	extern double mu;
	srand(202);
	double tmp[N];
	for (mu=mua;mu<mub;mu+=dmu)
	{
		for (int i=0;i<100;i++)
		{
			tmp[0]=lmin+(lmax-lmin)*((double)rand()/RAND_MAX);
			tmp[1]= -4+8*((double)rand()/RAND_MAX);
			tmp[2]=phimin+(phimax-phimin)*((double)rand()/RAND_MAX);
			tmp[3]=-4+8*((double)rand()/RAND_MAX);
			vec x(N,tmp);
			cout<<"#";x.show();

			plot_asymptote(&x);
		}
	}
}

