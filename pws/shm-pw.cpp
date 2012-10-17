#include <rk4.h>
#define R 28
#define B 2.6667
#define SIGMA 10
#define K2 1.5
#define N 3

using namespace std;

void f(double t, vec x, vec *out) 
{
	double x0,x1,x2;
	out->arr[0]=(x.arr[1]-x.arr[0])*SIGMA;

	out->arr[1]=x.arr[0]*(R-x.arr[2])-x.arr[1];

	out->arr[2]=(x.arr[0])*(x.arr[1])-(x.arr[2])*B;
}


int main(int argc , char *argv[])
{
	double t, tmp[N],tmax;
	int i;	
	for (i=0;i<N;i++) { tmp[i]=atof(argv[1+i]);}

	tmax=atof(argv[1+i]);
	t=0;
	
	vec x(N,tmp);	

	while (t<tmax)
	{
		x.show();
		rk4(t,&x);
		t+=h;
	}
}
