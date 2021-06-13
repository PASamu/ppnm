#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>

struct params {int n, counts; unsigned int seed;};

void* throw_points(void* arg){
	struct params * p=(struct params*)arg;
	p->counts=0;
	for(int i=0;i < p->n;i++){
		double x=(double)rand_r(&(p->seed))/RAND_MAX;
		double y=(double)rand_r(&(p->seed))/RAND_MAX;
		if(x*x+y*y<1)p->counts++;
	}
	return NULL;
}

int main(int argc, char** argv){
	int totalN=(int)1e6;
	if(argc>1) totalN=(int)atof(argv[1]);
	struct params n1 = {.n=totalN/3, .counts=0, .seed=1};
	struct params n2 = {.n=totalN/3, .counts=0, .seed=13};
	struct params n3 = {.n=totalN/3, .counts=0, .seed=42};
	pthread_t t1,t2,t3;
	pthread_create(&t1,NULL,throw_points,(void*)&n1);
	pthread_create(&t2,NULL,throw_points,(void*)&n2);
	pthread_create(&t3,NULL,throw_points,(void*)&n3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	int Nin =n1.counts+n2.counts+n3.counts;
	int N=n1.n+n2.n+n3.n;
	double pi=4*(double)Nin/N;
	printf("total number of points = %i\n",N);
	printf("pi from monte carlo sim = %g\n",pi);
	printf("simulated pi - actual pi = %g\n",fabs(pi-M_PI));
return 0;
}

