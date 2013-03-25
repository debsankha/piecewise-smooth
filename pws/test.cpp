#include<iostream>
#include<stdlib.h>
#include<hardcol.h>
using namespace std;

main(int argc, char **argv)
{
	srand(0);
	int i=0;
	
	float min=atof(argv[2]);
	float max=atof(argv[3]);
	
	double xini;

	while (i<atoi(argv[1])) 
	{
		xini=randdouble(min,max);
		cout<<xini<<endl;
		i++;
	}
}

