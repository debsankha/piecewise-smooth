#include <rk4.h>

using namespace std;
vec::vec(){double newarr[3];arr=&newarr[0];};
vec::vec(double arg[])
{
	arr=&arg[0];
}

vec vec::operator +(vec other)
{
	vec temp=vec();
	for (int i=0;i<N;i++) {temp.arr[i]=arr[i]+other.arr[i];}
	return temp;
}

vec vec::operator %(double other)
{
	vec temp=vec();
	for (int i=0;i<N;i++) {temp.arr[i]=arr[i]*other;}
	return temp;
}
vec vec::operator /(double other)
{
	vec temp=vec();
	for (int i=0;i<N;i++) {temp.arr[i]=arr[i]/other;}
	return temp;
}

void vec::show()
{
	for (int i=0;i<N;i++) {cout<<arr[i]<<"\t";}
	cout<<endl;
}

void rk4(double t, vec *x)
{
	vec k1,k2,k3,k4;
	f(t,*x,&k1);k1=k1%h;
	f(t+h/2,*x+k1/2,&k2);k2=k2%h;
	f(t+h/2,*x+k2/2,&k3);k3=k3%h;
	f(t+h,*x+k3,&k4);k4=k4%h;
	*x=*x+(k1+k2%2+k3%2+k4)/6;	
}
