#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include<time.h>

void printmatrix (float **A, int N)
{
int i,j;
for(i=1;i<=N;i++){
	for(j=1;j<=N;j++){
		printf("%f\t",A[i][j]);
		}
		printf("\n");
	}
}

void uppertriangular(float **S, float **A, int N)
{
int i,j;
for(i=1;i<=N;i++){
	for(j=1;j<=N;j++){
		if(i<=j)
			A[i][j] = S[i][j];
		else
			A[i][j]= 0.0;
		}
	}
}

void lowertriangular(float **A, float **B, float p[], int N)
{
int i,j;
for(i=1;i<=N;i++){
	for(j=1;j<=N;j++){
		if(i>j)
			A[i][j] = B[i][j];
		else if (i==j)
			A[i][j] = p[i];
		else
			A[i][j] = 0;
		}
	}
}

void transpose(float **A, float **B, int N)
{
int i,j;
for(i=1;i<=N;i++)
	for(j=1;j<=N;j++)
		A[i][j] = B[j][i];
}

void inverse(float **inv, float **a, int N)
{
int i,j;
float d,*col;
int *indx;
indx = ivector(1,N);
col = vector(1,N);
ludcmp(a,N,indx,&d);

for(j=1;j<=N;j++) {
	for(i=1;i<=N;i++) col[i]=0.0;
		col[j]=1.0;
lubksb(a,N,indx,col);
for(i=1;i<=N;i++) inv[i][j]=col[i];
}
free_vector(col,1,N);
free_ivector(indx,1,N);
}

void matmult(float **A, float **B, float **C, int N)
{
int i,j,k;
float sum;
for(i=1;i<=N;i++){
	for(j=1;j<=N;j++){
		sum=0.0;
		for(k=1;k<=N;k++){
			sum = sum + A[i][k]*B[k][j];}
		C[i][j]=sum;		
		}
	}
}


int main()
{
clock_t start = clock();
int N,i,j,m,n;
printf("Enter the number of basis vectors:");
scanf("%d",&N);
float **H, **S, **L, **Lt, **U;
H = matrix(1,N,1,N);
S = matrix(1,N,1,N);
U = matrix(1,N,1,N);
L = matrix(1,N,1,N);
Lt = matrix(1,N,1,N);
for(m=0;m<N;m++){
	for(n=0;n<N;n++){
		if((m+n)%2==0)
			S[m+1][n+1] = 2.0/(n+m+5.0) - 4.0/(m+n+3.0) + 2.0/(m+n+1);
		else
			S[m+1][n+1]=0;
		}
	}
for(m=0;m<N;m++){
	for(n=0;n<N;n++){
		if((m+n)%2==0)
			H[m+1][n+1] = -8.0*(1-m-n-2*m*n)/((m+n+3)*(m+n+1)*(m+n-1));
		else
			H[m+1][n+1]=0;
		}
	}

uppertriangular(S,U,N);
float *p;
p = vector(1,N);
choldc(U,N,p);
lowertriangular(L,U,p,N);
transpose(Lt,L,N);

float **y,**z;
y = matrix(1,N,1,N);
z = matrix(1,N,1,N);
inverse(y,L,N);   // Inverse of L
inverse(z,Lt,N);  // Inverse of Transpose[L]
float **temp, **C;
temp = matrix(1,N,1,N); C = matrix(1,N,1,N);
matmult(y,H,temp,N);
matmult(temp,z,C,N);

float *d,*e;
d = vector(1,N);
e = vector(1,N);
tred2(C,N,d,e);
tqli(d,e,N,C);
eigsrt(d,C,N);
for(i=1;i<=5;i++)
	printf("Eigenvalue %d = \t %f\n",i,d[i]);
//for free....................................................................
free_matrix(H,1,N,1,N);
free_matrix(S,1,N,1,N);
free_matrix(L,1,N,1,N);
free_matrix(Lt,1,N,1,N);
free_matrix(U,1,N,1,N);
free_matrix(y,1,N,1,N);
free_matrix(z,1,N,1,N);
free_matrix(temp,1,N,1,N);
free_matrix(C,1,N,1,N);
free_vector(p,1,N);
free_vector(d,1,N);
free_vector(e,1,N);

clock_t end = clock();
float seconds = (float)(end - start) / CLOCKS_PER_SEC;
printf("Execution time = %f seconds\n",seconds);
return 0;
}
#undef NRANSI
