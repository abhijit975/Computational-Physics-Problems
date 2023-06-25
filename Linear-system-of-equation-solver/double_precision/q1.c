#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include<time.h>

void printmatrix (double **A, int N)
{
int i,j;
for(i=1;i<=N;i++){
	for(j=1;j<=N;j++){
		printf("%.4f\t",A[i][j]);
		}
		printf("\n");
	}
}

void uppertriangular(double **S, double **A, int N)
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

void lowertriangular(double **A, double **B, double p[], int N)
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

void transpose(double **A, double **B, int N)
{
int i,j;
for(i=1;i<=N;i++)
	for(j=1;j<=N;j++)
		A[i][j] = B[j][i];
}

void inverse(double **inv, double **a, int N)
{
int i,j;
double d,*col;
int *indx;
indx = ivector(1,N);
col = dvector(1,N);
ludcmp(a,N,indx,&d);

for(j=1;j<=N;j++) {
	for(i=1;i<=N;i++) col[i]=0.0;
		col[j]=1.0;
lubksb(a,N,indx,col);
for(i=1;i<=N;i++) inv[i][j]=col[i];
}
free_dvector(col,1,N);
free_ivector(indx,1,N);
}

void matmult(double **A, double **B, double **C, int N)
{
int i,j,k;
double sum;
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
double **H, **S, **L, **Lt, **U;
H = dmatrix(1,N,1,N);
S = dmatrix(1,N,1,N);
U = dmatrix(1,N,1,N);
L = dmatrix(1,N,1,N);
Lt = dmatrix(1,N,1,N);
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
double *p;
p = dvector(1,N);
choldc(U,N,p);
lowertriangular(L,U,p,N);
transpose(Lt,L,N);

double **y,**z;
y = dmatrix(1,N,1,N);
z = dmatrix(1,N,1,N);
inverse(y,L,N);  
inverse(z,Lt,N);  
double **temp, **C;
temp = dmatrix(1,N,1,N); C = dmatrix(1,N,1,N);
matmult(y,H,temp,N);
matmult(temp,z,C,N);

double *d,*e;
d = dvector(1,N);
e = dvector(1,N);
tred2(C,N,d,e);
tqli(d,e,N,C);
eigsrt(d,C,N);
for(i=1;i<=5;i++)
	printf("Eigenvalue %d = \t %.6f\n",i,d[i]);
//for free....................................................................
free_dmatrix(H,1,N,1,N);
free_dmatrix(S,1,N,1,N);
free_dmatrix(L,1,N,1,N);
free_dmatrix(Lt,1,N,1,N);
free_dmatrix(U,1,N,1,N);
free_dmatrix(y,1,N,1,N);
free_dmatrix(z,1,N,1,N);
free_dmatrix(temp,1,N,1,N);
free_dmatrix(C,1,N,1,N);
free_dvector(p,1,N);
free_dvector(d,1,N);
free_dvector(e,1,N);

clock_t end = clock();
double seconds = (double)(end - start) / CLOCKS_PER_SEC;
printf("Execution time = %lf seconds\n",seconds);
return 0;
}
#undef NRANSI
