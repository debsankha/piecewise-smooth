#include <hardcol.h>
#define NPTSEACHF 10	//Number of points for each F value


using namespace std;

extern float Sigma,G,F,W,m,K1;

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
		for (int cnt=0;cnt<NPTSEACHF;cnt++)
		{
			double tmp[2];
			tmp[0]=randdouble(-Sigma*1.1,-Sigma*0.8);		//This choice of initial pts will have same impact 
			tmp[1]=pow(K1*(Sigma*1.1+tmp[0])*(Sigma*1.1-tmp[0]),0.5);//velocity, barring impacts

			vec x(2,tmp);
			cout<<F<<'\t'<<time_to_stabilize(x, tmax)<<endl;
			cerr<<"\tstarting from: "<<tmp[0]<<'\t'<<tmp[1]<<endl;
		}
	}
}



