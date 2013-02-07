#include<iostream>
#include<stdlib.h>
#define randdouble(min,max) min+rand()*(max-(min))*1.0/RAND_MAX

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

