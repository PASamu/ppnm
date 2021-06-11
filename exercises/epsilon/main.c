#include<stdio.h>
#include<limits.h>
#include<float.h>

int main(){
{
	//Checking maximum integer values for different loops
	int i=1; while(i+1>i) {i++;}
	printf("my max int (while) = %i\n",i);

	for(i=1;i+1>i;) {i++;}
	printf("my max int (for)  = %i\n",i);

	i=1; do {i++;} while(i+1>i);
	printf("my max int (do while) = %i\n",i);

	i=INT_MAX;
	printf("int_max = %i\n",i);
}
{
	//Checking minimum integer value for different loops
	int i=1; while(i-1<i) {i--;}
	printf("my min int (while) = %i\n",i);

	for(i=1;i-1<i;) {i--;}
	printf("my min int (for) = %i\n",i);

	i=1; do {i--;} while(i-1<i);
 	printf("my min int (do while) = %i\n",i);

	i=INT_MIN;
	printf("int_min = %i\n",i);
}
{
	printf("machine epsilon while loops\n");
	double x=1; while(1+x!=1){x/=2;} x*=2;
	printf("double = %g\n",x);

	float y=1; while(1+y!=1){y/=2;} y*=2;
	printf("float = %g\n",y);

	long double z=1; while(1+z!=1){z/=2;} z*=2;
	printf("long double = %Lg\n",z);
}
{
	printf("machine epsilon for loops\n");
	double x; for(x=1; 1+x!=1; x/=2){} x*=2;
	printf("double = %g\n",x);

	float y; for(y=1; 1+y!=1; y/=2){} y*=2;
	printf("float = %g\n",y);
	
	long double z; for(z=1; 1+z!=1; z/=2){} z*=2;
	printf("long double = %Lg\n",z);
}
{
	printf("machine epsilon do while loops\n");
	double x=1; do{x/=2;} while(1+x!=1); x*=2;
	printf("double = %g\n",x);

	float y=1; do{y/=2;} while(1+y!=1); y*=2;
	printf("float = %g\n",y);

	long double z=1; do{z/=2;} while(1+z!=1); z*=2;
	printf("long double = %Lg\n",z);

	printf("flt_epsilon=%g\n",FLT_EPSILON);
	printf("dbl_epsilon=%g\n",DBL_EPSILON);
	printf("ldbl_epsilon=%Lg\n",LDBL_EPSILON);
}
{
	int max=INT_MAX/4;
	int i;
	float sum_up_float=0;
	for(i=1; i<=max;i++)
	{sum_up_float+=1.0f/i;}
	printf("sum_up_float=%g\n",sum_up_float);

	float sum_down_float=0;
	for(i=max; i>0;i--)
	{sum_down_float+=1.0f/i;}
	printf("sum_down_float=%g\n",sum_down_float);
	printf("diff of sum_down and sum_up=%g\n",sum_up_float-sum_down_float);

	printf("When changing the value of max, the difference stays the same. So I do not see a convergence.\n");
}
{
	int max=INT_MAX/3;
	int i;
	double sum_up_double=0;
	for (i=1;i<=max;i++)
	{sum_up_double+=1.0/i;}
	printf("sum_up_double=%g\n",sum_up_double);

	double sum_down_double=0;
	for (i=max; i>0;i--)
	{sum_down_double+=1.0/i;} 
	printf("sum_down_double=%g\n",sum_down_double);

	printf("When using double instead of float, the sums are the same. This makes sense, since double has a higher precision (more decimals) than float\n");
}
}
