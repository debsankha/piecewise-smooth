#include <rk4.h>
#include <cmath>
#define N 2
#define Sigma 5
#define K1 4
#define K2 8

using namespace std;

void f(double t, vec x, vec *out) 
{
	//SHM hitting with a bouncy wall 
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0];//-G*x.arr[1];
	if (x.arr[0]>Sigma) out->arr[1]-=K2*x.arr[0];
}


int main(int argc , char *argv[])
{
	double t, tmp[N],tmax;
	double tmp2[]={0,12};
	int i;	
	for (i=0;i<N;i++) { tmp[i]=atof(argv[1+i]);}

	tmax=0.56;
	t=0;
	
	vec x(N,tmp);	
	vec x0(N,tmp2);
//	(x-x0).show();
	while (t<tmax)
	{
		rk4(t,&x);
		rk4(t,&x0);
		t+=h;
		cout<<t<<'\t'
	}
//	(x-x0).show();
}
