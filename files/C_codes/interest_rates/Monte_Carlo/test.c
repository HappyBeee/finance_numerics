#include <stdio.h>
#include <stdlib.h>

void function1(const void * ptr)
{
	double * p = (double *)ptr;
	p[1] = 3.14;
	printf("%.5lf \n", ((double *)ptr)[1]);
}

int main()
{
	const double * array1 = malloc(5*8);
	for (int i=0; i<5; i++)
	{
		array1[i] = i;
	}
	double nums[] = {2.0, 4.0, 6.0, 8.0, 10.0};
	const double * array2 = nums;

	function1(array1);
	function1(array2);
}

