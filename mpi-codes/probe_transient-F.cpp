#include <hardcol.h>
#define NPTSEACHX 300	//Number of points for each F value
#define NPTS 50	//Number of x values
#include <mpi.h>
#include <ctime>
#define ORB 4
#define EPSILON 0.0001
//time negative for some conditions. FIX ASAP
using namespace std;

<<<<<<< HEAD
extern float Sigma,G,F,W,m;
=======
extern float Sigma,G,F,W,m,K1;
>>>>>>> 891ac1e3d951fca45ce5ad384c6ed7798cd5f0f9
double time_to_stabilize(vec, double, int *);

int main(int argc, char **argv)
{
	srand(time(NULL));
	int rank,size;	//roughly speaking, #instance,#totalinstances
	MPI::Init(argc, argv);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();

	//Look at time to stabilize, while F is changed gradually to grazing
	

	float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
	cout<<"#F_graz: "<<F_graz<<endl;

	float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
	double F_min=atof(argv[1]);
	double F_max=atof(argv[2]);
	double tmax=atof(argv[3]);
		
	double F_range=F_max-F_min;

	double dF=F_range/NPTS;
	double stopF=F_max-rank*F_range/size;
	double startF=F_max-(rank+1)*F_range/size;
	
	bool timeceilreached=0;
	double tau;
	int period;

	for (F=startF;(F<stopF) && (timeceilreached==0);F+=dF)
	{
		for (int cnt=0;(cnt<NPTSEACHX) && (timeceilreached==0);cnt++)
		{
			double tmp[2];
			tmp[0]=randdouble(-8,1);//randdouble(-Sigma*1.1,-Sigma*0.8);		//This choice of initial pts will have same impact 
			tmp[1]=randdouble(-8,8);//pow(K1*(Sigma*1.1+tmp[0])*(Sigma*1.1-tmp[0]),0.5);//velocity, barring impacts

			vec x(2,tmp);
			tau=time_to_stabilize(x, tmax, &period);

			if ((tau>0) && (period==1))
			{
				cout<<F<<'\t'<<tau<<endl;
			}
		//	else timeceilreached=1;
		}
	}

	MPI::Finalize();
	return 0;
}



double time_to_stabilize(vec x, double tmax, int *per)
{
	double time_to_stable=0;
	double t=0;
	double oldvel=0;
	double poinc_x[ORB];							//array to store poincare x vals in to detect period
	int i=0;
	bool period;			
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
			i++;
		
			if (i==ORB)
			{
				period=pow((pow(poinc_x[0],2)+pow(poinc_x[1],2)+pow(poinc_x[2],2)+pow(poinc_x[3],2))/4.0-pow((poinc_x[0]+poinc_x[1]+poinc_x[2]+poinc_x[3])/4.0,2),0.5)<EPSILON;
				
				if (period==1)
				{
					*per=period;
					cerr<<"GOT PERIOD\n";
					return t;
				}
				i=0;
			}
		}
		t+=h;
	}
	return -1;		//Remember: you'd get t=0 if chaotic orb or periodicity >48
}
