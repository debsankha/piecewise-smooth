#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<cmath>
#define PI 3.141592653589793
#define NUM_ITER 1000
#define NUM_DISCARD 128
#define NUM_KEEP 64	//Better be 2**n so that upto nth compositions can be studied easily

using namespace std;
float bifurc(float f(float,float),float,float,float,float,int);

float f(float x,float mu)
{
//	return mu*x*(1-x);//	logistic
//	return (1+mu)*x-pow(x,3);
//	return mu-pow(x,2);
	
	float r=x+0.37-mu*sin(2*PI*x)/(2*PI);

	return r-floor(r);
}

main(int argc, char **argv)
{
	bifurc(f,atof(argv[1]),atof(argv[2]),atof(argv[3]),atof(argv[4]),atoi(argv[5]));	//5 args: xmin,xmax,pmin,pmax,num_same	
}


float bifurc(float f(float,float),float xmin,float xmax,float pmin,float pmax,int num_same)
{
	srand(time(NULL));
	float x;
	float dp=(pmax-pmin)/NUM_ITER;

	for (float p=pmin;p<pmax;p+=dp)
	{
		srand(20);
		for (int num=0;num<num_same;num++)
		{
			int j;
			x=xmin+((float)rand()/RAND_MAX)*(xmax-xmin);
			for (j=0;j<NUM_DISCARD;j++){x=f(x,p);}
			for (j=0;j<NUM_KEEP;j++)
			{
				x=f(x,p);
				cout<<p<<"\t"<<x<<endl;
			}
		}
	}
}

