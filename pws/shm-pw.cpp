#include <rk4.h>
#define K1 1
#define SIGMA 2
#define K2 1.5


using namespace std;

void f(double t, vec x, vec *out) 
{
	double args[3],x0,x1,x2;
	args[0]=(x.arr[1]-x.arr[0])*SIGMA;

	args[1]=x.arr[0]*(R-x.arr[2])-x.arr[1];

	args[2]=(x.arr[0])*(x.arr[1])-(x.arr[2])*B;

	out->arr=args;
}


int main(int argc , char *argv[])
{
	double t, tmp[N],tmax;
	int i;	
	for (i=0;i<N;i++) { tmp[i]=atof(argv[1+i]);}

	tmax=atof(argv[1+i]);
	t=0;
	
	vec x;	
	x=vec(tmp);

	x.show();
	while (t<tmax)
	{
		x.show();
		rk4(t,&x);
		t+=h;
	}
}
