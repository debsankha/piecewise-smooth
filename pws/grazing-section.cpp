#include <rk4.h>
#include <cmath>
#define N 2		//dimensionality of system
#define K1 1		
#define K2 4
#define G 0.001		//damping
#define G2 6.2			//extra damping
#define W 2.176		//forcing freq: sin(W*t)
#define ORB 49
#include <limits>
#include<stdlib.h>

using namespace std;

double F=1.4881;		//forcing amplitude
double Sigma=1.23;
void f(double t, vec x, vec *out) 
{
	//SHM hitting with a bouncy wall , x=(x,v,t)
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0]-G*x.arr[1]+F*sin(W*t);
	if (x.arr[0]>Sigma) out->arr[1]-=K2*(x.arr[0]-Sigma);//out->arr[1]-=G2*x.arr[1];}
}




void next(vec x)
{
	double a,b;
	b=pow(4*K1-pow(G,2),0.5)/2;
	a=G/2;
	

}


void printpoinc(double varf)
{
	F=varf;
	Sigma=F/pow(pow((K1-pow(W,2)),2)+pow(W*G,2),0.5);
	cerr<<"SIgma= "<<Sigma<<endl;
	double t, tmp[N],tmax;
	tmp[0]=-Sigma*1.01;tmp[1]=0;	
	tmax=10000;
	t=0;
	
	vec x(N,tmp);	
	
}



int main(int argc, char **argv)
{
	double f=atof(argv[1]);
	double fmax=atof(argv[2]);
	double df=atof(argv[3]);


	while (f<fmax) {printpoinc(f);f+=df;}
}

