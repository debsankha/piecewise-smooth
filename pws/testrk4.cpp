#include <rk4.h>
#include <cmath>
#define N 2
#define K1 4

using namespace std;

void f(double t, vec x, vec *out) 
{
	//SHM
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0];//-G*x.arr[1];
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
		rk4(t,&x);
		t+=h;
		cout<<t<<'\t'<<x.arr[0]-cos(2*t)*tmp[0]-sin(2*t)*tmp[1]/2<<endl;
			
	}
}