#include <rk4.h>
#include <cmath>
#define N 4		//dimensionality of system
#define mu 0.15

//x=(l,l',phi,phi')

using namespace std;

void f(double t, vec x, vec *out) 
{
	//Pulley and swinging ball. See Tel & Gruiz Sec 7.4.3
	out->arr[0]=x.arr[1];
	out->arr[2]=x.arr[3];
	out->arr[1]=mu*(cos(x.arr[2])+x.arr[0]*pow(x.arr[3],2)+1)-1;
	out->arr[3]=-(sin(x.arr[2])+2*x.arr[1]*x.arr[3])/x.arr[0];
}


int main(int argc , char *argv[])
{
	double t, tmp[N],tmax;
	int i;	
	for (i=0;i<N;i++) {tmp[i]=atof(argv[1+i]);}	//take x from stdin

	tmax=atof(argv[1+i]);
	t=0;
	
	vec x(N,tmp);
	
	double oldsinphi=sin(tmp[2]);
	while (t<tmax)
	{
		rk4(t,&x);
		t+=h;
		if ((oldsinphi*sin(x.arr[2])<0)&&(x.arr[3]>0)) cerr<<t<<'\t'<<x.arr[0]<<'\t'<<x.arr[1]<<endl;
		cout<<t<<'\t'<<x.arr[0]<<'\t'<<x.arr[1]<<'\t'<<x.arr[2]<<'\t'<<x.arr[3]<<endl;
		oldsinphi=sin(x.arr[2]);
	}
}
