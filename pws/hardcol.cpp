#include <hardcol.h>
using namespace std;

int main(int argc , char *argv[])
{
	srand(0);
	if (strcmp(argv[1],"traj")==0)
	{
		int i;	
		double tmax,tmp[N];
	
		for (i=1;i<N+1;i++) { tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
		vec x(N,tmp);	
		tmax=atof(argv[1+i]);
		F=atof(argv[2+i]);
		plottraj(x,tmax);
	}
	if (strcmp(argv[1],"poincare")==0)
	{
		int i;	
		double tmax,tmp[N];
	
		for (i=1;i<N+1;i++) { tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
		vec x(N,tmp);
	
		tmax=atof(argv[1+i]);
		F=atof(argv[2+i]);
		plotpoincare(x,0,tmax,tmax*0.9);
	}

	if (strcmp(argv[1],"plotbifurc")==0)
	{
		float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
		cerr<<"F_graz="<<F_graz<<endl;
		plotbifurc_F(F_graz-atof(argv[2]),F_graz+atof(argv[3]),atoi(argv[4]));
	}

}
