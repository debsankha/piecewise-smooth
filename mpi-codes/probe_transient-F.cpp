#include <hardcol.h>
#define NPTSEACHX 300	//Number of points for each F value
#define NPTS 50	//Number of x values
#include <mpi.h>
#include <ctime>
//time negative for some conditions. FIX ASAP
using namespace std;

extern float Sigma,G,F,W,m,K1;
double time_to_stabilize(vec, double);

int main(int argc, char **argv)
{
	srand(time(NULL));
	int rank,size;	//roughly speaking, #instance,#totalinstances
	MPI::Init(argc, argv);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();

	//Look at time to stabilize, while F is changed gradually to grazing
	

	cout<<"#F"<<'\t'<<"time_to_stabilize"<<endl;
	float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
	cout<<"#F_graz: "<<F_graz<<endl;

	double F_min=atof(argv[1]);
	double F_max=atof(argv[2]);
	double tmax=atof(argv[3]);
	
	double F_range=F_max-F_min;

	double dF=F_range/NPTS;
	double stopF=F_max-rank*F_range/size;
	double startF=F_max-(rank+1)*F_range/size;
	
	bool timeceilreached=0;
	double tau;


	for (F=startF;(F<stopF) && (timeceilreached==0);F+=dF)
	{
		for (int cnt=0;(cnt<NPTSEACHX) && (timeceilreached==0);cnt++)
		{
			double tmp[2];
			tmp[0]=randdouble(-Sigma*1.1,-Sigma*0.8);		//This choice of initial pts will have same impact 
			tmp[1]=pow(K1*(Sigma*1.1+tmp[0])*(Sigma*1.1-tmp[0]),0.5);//velocity, barring impacts

			vec x(2,tmp);
			tau=time_to_stabilize(x, tmax);

			if (tau>0)
			{
				cout<<F<<'\t'<<tau<<endl;
				cout<<"#starting from: "<<tmp[0]<<'\t'<<tmp[1]<<endl;
			}
			else timeceilreached=1;
		}
	}

	MPI::Finalize();
	return 0;
}



double time_to_stabilize(vec x, double tmax)
{
	double time_to_stable=0;
	double t=0;
	double oldvel=0;
	double poinc_x[ORB];							//array to store poincare x vals in to detect period
	double poinc_t[ORB];							//array to store poincare x vals in to detect period
	int i=0;
	int period=0;
			
	t=0;
	time_to_stable=0;
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
					cout<<"#F: "<<F<<"\tPeriod: "<<period<<endl;
					return time_to_stable;
				}
				i=0;
			}
		}
		t+=h;
	}
	return -1;		//Remember: you'd get t=0 if chaotic orb or periodicity >48
}
