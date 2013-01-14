#include <hardcol.h>

using namespace std;

extern float Sigma,G,F,W,m,K1;
double time_to_stabilize(vec, double);

int main(int argc, char **argv)
{

	//Look at time to stabilize, while F is changed gradually to grazing

	cout<<"#F"<<'\t'<<"time_to_stabilize"<<endl;
	float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
	cerr<<"#F_graz: "<<F_graz<<endl;

	double F_range=atof(argv[1]);
	double tmax=atof(argv[2]);

	double dF=F_range/100;

	for (F=F_graz-F_range;F<F_graz;F+=dF)
	{
		for (int cnt=0;cnt<10;cnt++)
		{
			double tmp[2];
			tmp[0]=randdouble(-Sigma*1.1,-Sigma*0.8);
			tmp[1]=pow(K1*(Sigma*1.1+tmp[0])*(Sigma*1.1-tmp[0]),0.5);

			vec x(2,tmp);
			cout<<F<<'\t'<<time_to_stabilize(x, tmax)<<endl;
			cerr<<"\tstarting from: "<<tmp[0]<<'\t'<<tmp[1]<<endl;
		}
	}
}



double time_to_stabilize(vec x, double tmax)
{
	double t=0;
	double time_to_stable=0;
	double oldvel=0;
	double poinc_x[ORB];							//array to store poincare x vals in to detect period
	double poinc_t[ORB];							//array to store poincare x vals in to detect period
	int i=0;
	int period=0;
			
	t=0;
	while (t<tmax)
	{
		oldvel=x.arr[1];
		rk4(t,&x);
		if (x.arr[0]>Sigma)						//The reset map on hard collision
		{
			x.arr[0]=Sigma;
			x.arr[1]*=-1;
		}

		if ((oldvel<0) && (x.arr[1]>0))
		{
			poinc_x[i]=x.arr[0];
			poinc_t[i]=t;
			i++;
		
			if (i==ORB)
			{
				period=detect_period(&poinc_x[0],&poinc_t[0], &time_to_stable);
				
				if (period>0)
				{
					cerr<<"#F: "<<F<<"\tPeriod: "<<period<<endl;
					return time_to_stable;
				}
				i=0;
			}
		}
		t+=h;
	}
	return 0;		//Remember: you'd get t=0 if chaotic orb or periodicity >48
}
