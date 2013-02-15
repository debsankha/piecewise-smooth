#include <hardcol.h>
using namespace std;

int main(int argc , char *argv[])
{
	srand(time(NULL));

	if (strcmp(argv[1],"plotmap")==0)
	{
		plotmap(atoi(argv[2]),atoi(argv[3]),atof(argv[4]),atof(argv[5]));
	}

	if (strcmp(argv[1],"plotmap3d")==0)
	{
		plotmap3d(atoi(argv[2]),atoi(argv[3]),atof(argv[4]),atof(argv[5]));
	}

	if (strcmp(argv[1],"poincare")==0)
	{
		double tmax,tmp[N];
	
		tmp[0]=atof(argv[2]);	//take x from stdin
		tmp[1]=0;
		vec x(N,tmp);
		tmax=atof(argv[3]);
		F=atof(argv[4]);
		G=atof(argv[5]);
		
		cerr<<"F_graz="<<Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5)<<endl;
		double initphase=(-M_PI/2-atan(G*W/(W*W-K1)))/W;
		plotpoincare(x,initphase,tmax);
	}

	if (strcmp(argv[1],"plotbasin")==0)
	{
		int npts=atoi(argv[2]);
		double tmax=atof(argv[3]);
		plotbasin(npts,tmax);
	}
}
