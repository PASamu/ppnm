#include"komplex.h"
#include"stdio.h"
#define TINY 1e-6

int main(){
	komplex a = {2,3}, b = {3,4};

	printf("testing komplex_add..\n");
	komplex r = komplex_add(a,b);
	komplex R = {5,7};
	komplex p = komplex_sub(a,b);
	komplex P = {-1,-1};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should   = ", R);
	komplex_print("a+b actually = ", r);
	komplex_print("a-b should   = ", P);
	komplex_print("a-b actually = ", p);
}
