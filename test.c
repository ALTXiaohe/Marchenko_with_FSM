#include "xiaohe.h"
int main(void)
{
	float complex *a, b;

	int N=100, i;
	float *real, imag, *real1;

	real = (float *) calloc (N, sizeof(float));
	real1 = (float *) calloc (N, sizeof(float));
	a = (float complex *) calloc (N, sizeof(float complex));
	for (i=0; i<N; i++)
	{
		real[i] = 1.11101*i;
	}
	DFFT(a, real, N);
	IDFFT(a, real1, N);
	b = 3 + 4 *I;

//	memset(a, 0, N*2*sizeof(float));
	for (i=0; i<N; i++)
	{
		printf("%f\t %f\n", crealf(a[i]), cimagf(a[i]));
	//	printf("%p\n", &cimag(a[1]));
	}
}
