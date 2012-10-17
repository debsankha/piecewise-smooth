#include <rk4.h>

using namespace std;
vec::vec(int N)
{
	len=N;
	arr=(double*)calloc(len,sizeof(double));
//	double newarr[3];
//	arr=&newarr[0];};
}
vec::vec(int N, double arg[])
{
	len=N;
	arr=(double*)calloc(len,sizeof(double));
	for (int i=0;i<len;i++) {*(arr+i)=arg[i];}
}

vec vec::operator +(const vec other)
{
	vec temp(len);
	for (int i=0;i<len;i++) {temp.arr[i]=arr[i]+other.arr[i];}
	return temp;
}

vec::vec(const vec& other)	//COPY CONST
{
	len=other.len;
	arr=(double*)calloc(len,sizeof(double));
	for (int i=0;i<len;i++) {*(arr+i)=other.arr[i];}
}

vec& vec::operator =(const vec& other)
{
	for (int i=0;i<len;i++) {arr[i]=other.arr[i];}
	return *this;
}


vec vec::operator %(const double other)
{
	vec temp(len);
	for (int i=0;i<len;i++) {temp.arr[i]=arr[i]*other;}
	return temp;
}
vec vec::operator /(const double other)
{
	vec temp(len);
	for (int i=0;i<len;i++) {temp.arr[i]=arr[i]/other;}
	return temp;
}

void vec::show()
{
	for (int i=0;i<len;i++) {cout<<*(arr+i)<<"\t";}
	cout<<endl;
}

vec::~vec()
{
	free(arr);
}

void rk4(double t, vec *x)
{
	int N=(*x).len;
	vec k1(N);
	vec k2(N);
	vec k3(N);
	vec k4(N);
	f(t,*x,&k1);
	k1=k1%h;
	f(t+h/2,*x+k1/2,&k2);k2=k2%h;
	f(t+h/2,*x+k2/2,&k3);k3=k3%h;
	f(t+h,*x+k3,&k4);k4=k4%h;
	*x=*x+(k1+k2%2+k3%2+k4)/6;	
}
