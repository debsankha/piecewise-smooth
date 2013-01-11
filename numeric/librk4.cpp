#include <rk4.h>

using namespace std;
vec::vec(int N)		
{
	len=N;
	arr=(double*)calloc(len,sizeof(double));
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

vec vec::operator -(const vec other)
{
	vec temp(len);
	for (int i=0;i<len;i++) {temp.arr[i]=arr[i]-other.arr[i];}
	return temp;
}


vec::vec(const vec& other)	//COPY CONST
{
	len=other.len;
	arr=(double*)calloc(len,sizeof(double));
	for (int i=0;i<len;i++) {*(arr+i)=other.arr[i];}
}

vec& vec::operator =(const vec& other)		//overloaded assignment
{	//A reasonably fail-safe deepcopy, can meddle with no-throw swap for better security
	if (this!=&other)
	{
	double *temp;
	temp=(double*)calloc(len,sizeof(double));
	
	for (int i=0;i<len;i++) {*(temp+i)=other.arr[i];}
	free(arr);

	arr=temp;
	}
	return *this;
}


vec vec::operator %(const double other)		//Scalar multiplication. 
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

double vec::norm()
{
	double res=0;
	for (int i=0;i<len;i++) {res+=(*(arr+i))*(*(arr+i));}
	return res;
}


vec::~vec()
{
	free(arr);
}

void rk4(double t, vec *x)
{	//Normal RK4 stuff
	int N=(*x).len;
	vec k1(N);
	vec k2(N);
	vec k3(N);
	vec k4(N);
	f(t,*x,&k1);k1=k1%h;
	f(t+h/2,*x+k1/2,&k2);k2=k2%h;
	f(t+h/2,*x+k2/2,&k3);k3=k3%h;
	f(t+h,*x+k3,&k4);k4=k4%h;
	*x=*x+(k1+k2%2+k3%2+k4)/6;
}
