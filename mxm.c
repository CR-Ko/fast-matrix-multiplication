#include <stdio.h>
#include <stdlib.h>

#define BLK_SIZE 128
#define min(a,b) (((a)<(b))?(a):(b))

void mxm_naive(double *C, const double *A, const double *B, int m, int n, int p)
{
    int i, j, k;

    for(i = 0; i < m; i++)
        for(j = 0; j < p; j++)
			for(k = 0; k < n; k++)
                C[p*i+j] += A[n*i+k] * B[p*k+j];			      
}

void mxm_block(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, bs = BLK_SIZE;

	for(ii = 0; ii < m; ii += bs)
		for(jj = 0; jj < p; jj += bs)
			for(kk = 0; kk < n; kk += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(j = jj; j < min(p, jj+bs); j++)
						for(k = kk; k < min(n, kk+bs); k++)
							C[p*i+j] += A[n*i+k] * B[p*k+j];											
}

void mxm_block_reorder(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, bs = BLK_SIZE;

	for(ii = 0; ii < m; ii += bs)
		for(kk = 0; kk < n; kk += bs)
			for(jj = 0; jj < p; jj += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(k = kk; k < min(n, kk+bs); k++)
						for(j = jj; j < min(p, jj+bs); j++)
							C[p*i+j] += A[n*i+k] * B[p*k+j];						
				
}

void mxm_block_reorder_reuse(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, Aik, bs = BLK_SIZE;
	
	for(ii = 0; ii < m; ii += bs)
		for(kk = 0; kk < n; kk += bs)
			for(jj = 0; jj < p; jj += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(k = kk; k < min(n, kk+bs); k++)
					{
						Aik = A[n*i+k];
						for(j = jj; j < min(p, jj+bs); j++)
							C[p*i+j] += Aik * B[p*k+j];		
					}					
}

void mxm_block_reorder_reuse_unroll_2(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, Aik, bs = BLK_SIZE;
	
	for(ii = 0; ii < m; ii += bs)
		for(kk = 0; kk < n; kk += bs)
			for(jj = 0; jj < p; jj += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(k = kk; k < min(n, kk+bs); k++)
					{
						Aik = A[n*i+k];
						for(j = jj; j < min(p, jj+bs); j+=2)
						{
							C[p*i+j] += Aik * B[p*k+j];
							C[p*i+j+1] += Aik * B[p*k+j+1];		
						}
					}					
}

void mxm_block_reorder_reuse_unroll_4(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, Aik, bs = BLK_SIZE;

	for(ii = 0; ii < m; ii += bs)
		for(kk = 0; kk < n; kk += bs)
			for(jj = 0; jj < p; jj += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(k = kk; k < min(n, kk+bs); k++)
					{
						Aik = A[n*i+k];
						for(j = jj; j < min(p, jj+bs); j+=4)
						{
							C[p*i+j] += Aik * B[p*k+j];
							C[p*i+j+1] += Aik * B[p*k+j+1];		
							C[p*i+j+2] += Aik * B[p*k+j+2];	
							C[p*i+j+3] += Aik * B[p*k+j+3];		
						}
					}					
}

void mxm_block_reorder_reuse_unroll_8(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, Aik, bs = BLK_SIZE;

	for(ii = 0; ii < m; ii += bs)
		for(kk = 0; kk < n; kk += bs)
			for(jj = 0; jj < p; jj += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(k = kk; k < min(n, kk+bs); k++)
					{
						Aik = A[n*i+k];
						for(j = jj; j < min(p, jj+bs); j+=8)
						{
							C[p*i+j] += Aik * B[p*k+j];
							C[p*i+j+1] += Aik * B[p*k+j+1];		
							C[p*i+j+2] += Aik * B[p*k+j+2];	
							C[p*i+j+3] += Aik * B[p*k+j+3];	
							C[p*i+j+4] += Aik * B[p*k+j+4];		
							C[p*i+j+5] += Aik * B[p*k+j+5];	
							C[p*i+j+6] += Aik * B[p*k+j+6];	
							C[p*i+j+7] += Aik * B[p*k+j+7];	
						}
					}					
}

void mxm_block_reorder_reuse_unroll_16(double *C, const double *A, const double *B, int m, int n, int p)
{ 
	int i, j, k, ii, jj, kk, Aik, bs = BLK_SIZE;

	for(ii = 0; ii < m; ii += bs)
		for(kk = 0; kk < n; kk += bs)
			for(jj = 0; jj < p; jj += bs)
				for(i = ii; i < min(m, ii+bs); i++)
					for(k = kk; k < min(n, kk+bs); k++)
					{
						Aik = A[n*i+k];
						for(j = jj; j < min(p, jj+bs); j+=16)
						{
							C[p*i+j] += Aik * B[p*k+j];
							C[p*i+j+1] += Aik * B[p*k+j+1];		
							C[p*i+j+2] += Aik * B[p*k+j+2];	
							C[p*i+j+3] += Aik * B[p*k+j+3];	
							C[p*i+j+4] += Aik * B[p*k+j+4];		
							C[p*i+j+5] += Aik * B[p*k+j+5];	
							C[p*i+j+6] += Aik * B[p*k+j+6];	
							C[p*i+j+7] += Aik * B[p*k+j+7];
							C[p*i+j+8] += Aik * B[p*k+j+8];		
							C[p*i+j+9] += Aik * B[p*k+j+9];	
							C[p*i+j+10] += Aik * B[p*k+j+10];	
							C[p*i+j+11] += Aik * B[p*k+j+11];		
							C[p*i+j+12] += Aik * B[p*k+j+12];	
							C[p*i+j+13] += Aik * B[p*k+j+13];	
							C[p*i+j+14] += Aik * B[p*k+j+14];	
							C[p*i+j+15] += Aik * B[p*k+j+15];	
						}
					}					
}

void mxm(double *C, const double *A, const double *B, int m, int n, int p, int uf)
{
	//mxm_naive(C, A, B, m, n, p);
	//mxm_block(C, A, B, m, n, p);
	//mxm_block_reorder(C, A, B, m, n, p);
	switch(uf)
	{
		case 16:mxm_block_reorder_reuse_unroll_16(C, A, B, m, n, p);
				break;
		case 8:	mxm_block_reorder_reuse_unroll_8(C, A, B, m, n, p);
				break;
  		case 4: mxm_block_reorder_reuse_unroll_4(C, A, B, m, n, p);
				break;
  		case 2: mxm_block_reorder_reuse_unroll_2(C, A, B, m, n, p);
				break;
  		default: mxm_block_reorder_reuse(C, A, B, m, n, p);
				break;
	}
}

// Cmxp = Amxn * Bnxp
void mxm_product(double *dest, const double *a, const double *b, int m, int n, int p)
{
	int uf; 
	
	if(p % 16 == 0) 		uf = 16;
	else if(p % 8 == 0) 	uf = 8;
	else if(p % 4 == 0) 	uf = 4;
	else if(p % 2 == 0)		uf = 2;
	else					uf = 0;
	
	printf("Matrix dimensions %dx%d\n", m, n);
	printf("Block size: %d\n", BLK_SIZE);	
	printf("Unrolls: %d\n", uf);
	
	mxm(dest, a, b, m, n, p, uf); 				
}
