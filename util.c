#include <stdio.h>
#include <stdlib.h>

#define min(a,b) (((a)<(b))?(a):(b))

void fill_mat(double* A, int m, int n)
{
	int i, j;
	
	for(i = 0; i < m; i++)
		for(j = 0; j < n; j++)	
			A[i*n+j] = rand() % 10;
}

void print_mat(double* A, int m, int n)
{
	int i, j;

	for(i = 0; i < m; i++)
	{
		for(j = 0; j < n; j++)
			printf("%.2lf ", A[i*n+j]);
		printf("\n");
	}
	printf("\n");
}

double sum_diag(double *A, int m, int n)
{
	int i;
	double sum = 0.0;

	for(i = 0; i < min(m,n); i++)
		sum += A[i*n+i];

	return sum;
}

