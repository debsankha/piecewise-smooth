#include <rk4.h>
#include <cmath>
#define N 2
#define Sigma 5
#define K1 4
#define K2 3
#define G 0
#define F 0
#define W 2

using namespace std;

void f(double t, vec x, vec *out) 
{
	//SHM hitting with a bouncy wall 
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0]-G*x.arr[1];
	if (x.arr[0]>Sigma) out->arr[1]-=K2*x.arr[0];
//	out->arr[2]=1;	
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
		cout<<t<<'\t';
		x.show();
		rk4(t,&x);
		t+=h;
	}
}
