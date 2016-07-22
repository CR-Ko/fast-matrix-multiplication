#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mxm.h"
#include "util.h"

#define N 1000
#define M 1000
#define P 1000
#define PRINT_MATS 0
#define PRINT_RESULT 0

void main()
{
	int m, n, p;	
	double *a, *b, *c;
	static struct timeval str, end;
	unsigned long long time;
	
	m = M;
	n = N;
	p = P;
	a = (double*)malloc((m * n) * sizeof(double));
	b = (double*)malloc((n * p) * sizeof(double));
	c = (double*)calloc(m * p, sizeof(double));
	
	fill_mat(a, m, n);
	fill_mat(b, n, p);
 	if(PRINT_MATS)
	{
		printf("A\n"); print_mat(a, m, n);
		printf("B\n"); print_mat(b, n, p);
	}

    	gettimeofday(&str, NULL);
	mxm_product(c, a, b, m, n, p);
    	gettimeofday(&end, NULL);
    
	time = 1000 * (end.tv_sec - str.tv_sec) + (end.tv_usec - str.tv_usec) / 1000;
	printf("Diagonal sum: %.2lf\n", sum_diag(c, m, p));
	printf("Time: %llu ms\n", time);
	if(PRINT_RESULT) print_mat(c, m, p);
	free(a); free(b); free(c);
}

