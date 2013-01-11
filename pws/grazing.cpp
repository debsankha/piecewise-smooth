#include <rk4.h>
#include <cmath>
#define N 2		//dimensionality of system
#define K1 1		
#define K2 4
#define G 0.001		//damping
#define G2 6.2			//extra damping
#define F 1.4881		//forcing amplitude
#define W 2.176		//forcing freq: sin(W*t)

using namespace std;

double Sigma=F/pow(pow((K1-pow(W,2)),2)+pow(W*G,2),0.5);

void f(double t, vec x, vec *out) 
{
	//SHM hitting with a bouncy wall , x=(x,v,t)
	out->arr[0]=x.arr[1];
	out->arr[1]=-K1*x.arr[0]-G*x.arr[1]+F*sin(W*t);
	if (x.arr[0]>Sigma) out->arr[1]-=K2*(x.arr[0]-Sigma);//out->arr[1]-=G2*x.arr[1];}
}


int main(int argc , char *argv[])
{
	cerr<<"Sigma= "<<Sigma<<endl;
	double t, tmp[N],tmax;
	int i=0;	
//	for (i=0;i<N;i++) { tmp[i]=atof(argv[1+i]);}	//take x from stdin
	double mult=atof(argv[1+i]);i++;	
	tmp[0]=-Sigma*mult;tmp[1]=0;	
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
