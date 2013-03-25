#include<iostream>
#define _USE_MATH_DEFINES
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<cstring>
#define randdouble(min,max) min+(max-(min))*((float) rand())/RAND_MAX
using namespace std;


int NPTS=300;
int mcpts=1000000;
float F=0.29;	//only a default value
float Gamma=0.062;	//only a default value
float m=1;	//mass
float sigma=1;	//switching manifold
float w=1;	//force freq
float n_val=2.3556;	//2*w_g/w
float K1=(w*n_val/2)*(w*n_val/2.0)+Gamma*Gamma/4.0;
float w_0=sqrt(K1/m);
float w_g=n_val*w/2;	//sqrt(w_0**2-Gamma**2/4)
float Ast=0;
float bigbox_area=0;
float bigellipse_area=0;

inline float C(float x, float v);
inline float B(float x, float v);
float nextcol(float x,float v);
float  vx_area(float f,float g);
void scan_f(float fmin,float fmax,float g,int npts);
void scan_g(float fmin,float fmax,float g,int npts);


int main(int argc, char **argv)
{
	srand(time(NULL));

	if (strcmp(argv[1],"f")==0) scan_f(atof(argv[2]),atof(argv[3]),atof(argv[4]),atoi(argv[5]));
	else
	{
		if (strcmp(argv[1],"g")==0) scan_g(atof(argv[2]),atof(argv[3]),atof(argv[4]),atoi(argv[5]));
		else cerr<<"invalid option, 1st argument must be f or g"<<endl;
	}

}


inline float C(float x, float v)
{
	return -atan(((Gamma/2*x+v)/w_g)/x);
}

inline float B(float x, float v)
{
	return pow((pow(x,2)+pow(((Gamma/2*x+v)/w_g),2)),0.5);
}

float nextcol(float x,float v)
{
	float c=C(x,v);
	float b=B(x,v);

	float res=(w_g>w)? (M_PI-c)/(w_g-w):(M_PI+c)/(w-w_g);

	float nextpeak_height=(Ast+b*exp(-Gamma*res/2));
	float heightnow=-Ast+b*cos(c);		//(Ast+b)*sin(c/2)
	float tc;

	if (nextpeak_height>heightnow)	//max height reached at the next peak of the envelope
	{
		tc=res;
		return nextpeak_height;

	}
	else
	{
		cerr<<"oh no!\n";				//max height reached at next peak of the high freq. signal itself, problem with calc
		tc=(2*M_PI)/(w+w_g);		//very rough estimate
		return Ast+b*exp(-Gamma*tc/2)*abs(sin((c+(w_g-w)*tc)/2.0));
	}
}	

float  vx_area(float f,float g)
{
	F=f;
	Gamma=g;
	K1=(w*n_val/2)*(w*n_val/2.0)+Gamma*Gamma/4.0;
	w_0=sqrt(K1/m);
	Ast=F/(m*pow(pow((pow(w,2)-pow(w_0,2)),2)+pow(Gamma*w,2),0.5));

	if (Ast>sigma) return 0;

	float chi_far=(sigma-Ast)*exp(3*Gamma*M_PI/(2*(w_g-w)));	//assume next peak at next envelope
	float chi_near=(sigma-Ast)*exp(2*Gamma*M_PI/(w_g+w));	//assume next peak at same envelope
	float chi_avg=(sigma-Ast)*exp(Gamma*M_PI/(2*(w_g-w)));	//(cmax+cmin)/2
	
	float chi=chi_avg;	
	float xmin=-1*chi;
	float xmax=(chi<1)?chi:1;
	
	bigbox_area=(xmax-xmin)*2*chi*w_g;
//	bigellipse_area=M_PI*w_g*pow(chi,2);
	bigellipse_area=w_g*((xmax*pow(chi*chi-xmax*xmax,0.5)+chi*chi*atan(xmax/pow(chi*chi-xmax*xmax,0.5))) 
			+chi*chi*M_PI/2);

	int sofar=0;
	int istrue=0;
	
	float x,ymin,ymax,y,nexthei;	
	while (sofar<mcpts)
	{
		x=randdouble(xmin,xmax);
		ymin=-1*chi*w_g-Gamma*x/2;
		ymax=chi*w_g-Gamma*x/2;
		y=randdouble(ymin,ymax);
	
		nexthei=nextcol(x,y);
		if (nexthei<sigma)
		{	
			istrue+=1;
		}
		sofar+=1;

	}

	return (istrue*1.0/mcpts)*bigbox_area;
}


void scan_f(float fmin,float fmax,float g,int npts)
{
	float f=fmin;
	float dF=(fmax-fmin)/npts;
	
	float area;
	while (f<fmax)
	{
		area=vx_area(f,g);
		cout<<f<<'\t'<<area<<'\t'<<bigellipse_area<<endl;
		f+=dF;
		cerr<<"F: "<<f<<endl;
	}

}

void scan_g(float gmin,float gmax,float f,int npts)
{
	float g=gmin;
	float dg=(gmax-gmin)/npts;
	
	float area;
	while (g<gmax)
	{
		area=vx_area(f,g);
		cout<<g<<'\t'<<area<<'\t'<<bigellipse_area<<endl;
		g+=dg;
		cerr<<"G: "<<g<<endl;
	}

}
