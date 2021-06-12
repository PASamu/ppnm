#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main (int argc, char** argv){
	if (argc<2)
	printf("There were no arguments\n");
	else{
		for(int i=1;i<argc;i++){
			double x = atof(argv[i]);
			printf("x=%g \t sin(x)=%g \t cos(x)=%g\n",x,sin(x),cos(x));
		}
	}
return 0;
}
