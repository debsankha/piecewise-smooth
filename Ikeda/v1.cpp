#include<iostream>
#include<stdlib.h>
#include<cmath>
#define h 0.01
#define N 3
#define B 8.0/3
#define SIGMA 10
#define R 26

using namespace std;

class Vector	// A simple column Vectortor with N elements, +,*,/ are overloaded
{
	public:
		double *arr;
		Vector()
		{
			double temp[N];
			for (int i=0;i<N;i++) {temp[i]=0;}
			arr=&temp[0];
		}

		Vector(double arg[])
		{
			(*this).arr=&arg[0];
		}
		
		Vector operator +(const Vector other)
		{
			double tarr[N];

			for (int i=0;i<N;i++)
			{
				tarr[i]=arr[i]+other.arr[i];
			}
			
			return Vector(tarr);
		}

		Vector operator %(const double other)
		{

			double tarr[N];
			for (int i=0;i<N;i++) {tarr[i]=(other)*(arr[i]);}
			return Vector(tarr);
		}
		Vector operator /(const double other)
		{

			double tarr[N];
			for (int i=0;i<N;i++) {tarr[i]=(1/other)*arr[i];}
			return Vector(tarr);
		}
		
//		Vector operator =(const Vector& other)
//		{
//			double *tarr;
//			tarr=new double[N];
//			for (int i=0;i<N;i++) {tarr[i]=arr[i];}
//			return Vector(tarr);
//		}


		void show()
		{
			for (int i=0;i<N;i++) {cout<<(this->arr)[i]<<"\t";}
			cout<<endl;
		}
};

void f(double, Vector, Vector*);	

void rk4(double ,Vector*);

int main(int argc , char *argv[])
{
	double t, tmp[N],tmax;
	int i;	
	for (i=0;i<N;i++) { tmp[i]=atof(argv[1+i]);}

	tmax=atof(argv[1+i]);
	t=0;
	
	Vector x=Vector(tmp);

	x.show();
	
//	Vector res;
//	res=x%0.01;
//
//	res.show();

	while (t<tmax)
	{
		x.show();
		rk4(t,&x);
		t+=h;
	}
}

void f(double t, Vector x, Vector *out) 
{
	double args[]={1,2,3};
	args[0]=(x.arr[1]-x.arr[0])*SIGMA*h;
	args[1]=(x.arr[0]*(R-x.arr[2])-x.arr[1])*h;
	args[2]=((x.arr[0])*(x.arr[1])-(x.arr[2])*B)*h;

//	cout<<args[0]<<"\t"<<args[1]<<"\t"<<args[2]<<endl;	

	for (int i=0;i<N;i++)
	{
//		cout<<"i=="<<i<<endl;
		(*out).arr[i]=args[i];	
//		cout<<args[i]<<endl;
	}
}

void rk4(double t, Vector *x)
{
	double a1[]={1,2,3};
	double a2[]={1,2,3};
	double a3[]={1,2,3};
	double a4[]={1,2,3};
	Vector k1=Vector(a1);
	Vector k2=Vector(a2);
	Vector k3=Vector(a3);
	Vector k4=Vector(a4);

	f(t,*x,&k1);
//	cout<<"k1: ";
//
//	k1.show();

	f(t+h/2,*x+k1/2,&k2);
//	cout<<"k2: ";
//	k2.show();

	f(t+h/2,*x+k2/2,&k3);
//	cout<<"k3: ";
//	k3.show();

	f(t+h,*x+k3,&k4);
	
//	cout<<"k4: ";
//	k4.show();
	Vector x1;
	*x=*x+(k1+k2%2+k3%2+k4)/6;
//	x1.show();
}
