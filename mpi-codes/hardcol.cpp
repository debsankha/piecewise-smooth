#include <hardcol.h>
#include <mpi.h>
using namespace std;

int main(int argc , char *argv[])
{
	int rank,size;	//roughly speaking, #instance,#totalinstances
	srand(time(NULL));
	MPI::Init(argc, argv);
	rank = MPI::COMM_WORLD.Get_rank();
	size = MPI::COMM_WORLD.Get_size();

	if (strcmp(argv[1],"plotbasin")==0)
	{
		int npts=atoi(argv[2]);
		double tmax=atof(argv[3]);
		plotbasin(npts,tmax,rank);
	}


	if (strcmp(argv[1],"traj")==0)
	{
		int i;	
		double tmax,tmp[N];
	
		for (i=1;i<N+1;i++) {tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
		vec x(N,tmp);	
		tmax=atof(argv[1+i]);
		F=atof(argv[2+i]);
		plottraj(x,tmax);
	}
	if (strcmp(argv[1],"poincare")==0)
	{
		int i;	
		double tmax,tmp[N];
	
		for (i=1;i<N+1;i++) {tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
		vec x(N,tmp);
	
		tmax=atof(argv[1+i]);
		F=atof(argv[2+i]);
		plotpoincare(x,0,tmax,tmax*0.9);
	}

	if (strcmp(argv[1],"plotbifurc")==0)
	{
		float K1=(W*m/2)*(W*m/2.0)+G*G/4.0;
		float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
		cerr<<"F_graz="<<F_graz<<endl;
		plotbifurc_F(F_graz-atof(argv[2]),F_graz+atof(argv[3]),atoi(argv[4]));
	}

	if (strcmp(argv[1],"ischaos_n")==0)
	{
		double n0=atof(argv[2]);
		double n1=atof(argv[3]);
		ischaos_n(n0+rank*(n1-n0)/((float)size), n0+(rank+1)*(n1-n0)/((float)size));
	}
}
