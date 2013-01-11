#include <rk4.h>
#include <cmath>
#define N 3		//dimensionality of system
#define Sigma 1		//the boundary: x=5
#define K1 1		
#define K2 1000
#define G 0.00062		//damping
#define G2 6.2			//extra damping
#define F 1.4881		//forcing amplitude
#define W 2.176
	//forcing freq: sin(W*t)

using namespace std;

void f(double t, vec x, vec *out) 
{
	//SHM hitting with a bouncy wall , x=(x,v,t)
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0]-G*x.arr[1]+F*sin(W*x.arr[2]);
	if (x.arr[0]>Sigma) {out->arr[1]-=K2*(x.arr[0]-Sigma);}//out->arr[1]-=G2*x.arr[1];}
	out->arr[2]=1;	
}


int main(int argc , char *argv[])
{
	double t, tmp[N],tmax;
	int i;	
	for (i=0;i<N;i++) { tmp[i]=atof(argv[1+i]);}	//take x from stdin

	tmax=atof(argv[1+i]);
	t=0;
	
	vec x(N,tmp);	

	while (t<tmax)
	{
		rk4(t,&x);
		t+=h;
		cout<<t<<'\t';x.show();			//responsible for the output
	}
}
