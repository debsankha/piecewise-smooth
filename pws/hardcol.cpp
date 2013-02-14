#include <hardcol.h>
using namespace std;

int main(int argc , char *argv[])
{
	srand(time(NULL));

	if (strcmp(argv[1],"plotmap")==0)
	{
		plotmap(atoi(argv[2]),atoi(argv[3]),atof(argv[4]),atof(argv[5]));
	}
	

	if (strcmp(argv[1],"traj")==0)
	{
		int i;	
		double tmax,tmp[N];
		
		if (argc>5)
		{
			for (i=1;i<N+1;i++) { tmp[i-1]=atof(argv[1+i]);}	//take x from stdin
			vec x(N,tmp);	
			tmax=atof(argv[1+i]);
			F=atof(argv[2+i]);
			plottraj(x,tmax);
		}

		else
		{
			tmp[0]=randdouble(-Sigma*1.1,-Sigma*0.8);		//This choice of initial pts will have same impact 
			tmp[1]=pow(K1*(Sigma*1.1+tmp[0])*(Sigma*1.1-tmp[0]),0.5);//velocity, barring impacts
			vec x(N,tmp);
			tmax=atof(argv[2]);
			F=atof(argv[3]);
			plottraj(x,tmax);

		}

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

	if (strcmp(argv[1],"plotbifurc")==0)
	{
		float F_graz=Sigma*pow(pow(W*W-K1,2)+W*G*W*G,0.5);
		cerr<<"F_graz="<<F_graz<<endl;
		plotbifurc_F(F_graz-atof(argv[2]),F_graz+atof(argv[3]),atoi(argv[4]));
	}

}
