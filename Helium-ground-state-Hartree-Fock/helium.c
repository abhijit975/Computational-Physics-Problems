#include <stdio.h>
#include<stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include <math.h>
#define pi acos(-1.0)

void normalize(double **a, double c[])
{
	double sum = 0.0,N;
	int i,j;
	for(i=1;i<=4;i++)
	{
		for(j=1;j<=4;j++)
		{
			sum += c[i]*a[i][j]*c[j];
		}
	}
	N = sqrt(sum);
	for(i=1;i<=4;i++)
		c[i]=c[i]/N;
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

void eigen(double **H, double **S, double d[], double eigv[], int N)
{
int i,j,m,n;
double **L, **Lt, **U;
U = dmatrix(1,N,1,N);
L = dmatrix(1,N,1,N);
Lt = dmatrix(1,N,1,N);
uppertriangular(S,U,N);
double *p;
p = dvector(1,N);
choldc(U,N,p);
lowertriangular(L,U,p,N);
transpose(Lt,L,N);
double **y,**z;
y = dmatrix(1,N,1,N);
z = dmatrix(1,N,1,N);
inverse(y,L,N);   // Inverse of L
inverse(z,Lt,N);  // Inverse of Transpose[L]
double **temp,**C;
C = dmatrix(1,N,1,N);
temp = dmatrix(1,N,1,N); 
matmult(y,H,temp,N);
matmult(temp,z,C,N);
double *e,val;
e = dvector(1,N);
tred2(C,N,d,e);
tqli(d,e,N,C);
eigsrt(d,C,N);
for(i=1;i<=N;i++)
	e[i]=C[i][1];
for(i=1;i<=N;i++){
	val = 0.0;
	for(j=1;j<=N;j++){
		val = val+z[i][j]*e[j];}
	eigv[i]=val;
	}
free_dmatrix(L,1,N,1,N);
free_dmatrix(Lt,1,N,1,N);
free_dmatrix(U,1,N,1,N);
free_dmatrix(y,1,N,1,N);
free_dmatrix(z,1,N,1,N);
free_dmatrix(temp,1,N,1,N);
free_dmatrix(C,1,N,1,N);
free_dvector(p,1,N);
free_dvector(e,1,N);
}

double difference(double a[], double b[])
{
double sum = 0.0;
for(int i=1;i<=4;i++){
	sum = sum + fabs(a[i]-b[i]);
	}
return sum;
}

//.................................................................................
int main()
{
double **S, **h, **F, *alpha, *c, *eigv;
int i,j,k,l;
double dif=1.0, E_0, sum, qpqrs, eps = 1.E-6;
S = dmatrix(1,4,1,4);
h = dmatrix(1,4,1,4);
F = dmatrix(1,4,1,4);
alpha = dvector(1,4);
eigv = dvector(1,4);
c = dvector(1,4);
//alpha[1]= 0.297104; alpha[2]=1.236745; alpha[3]=5.749982; alpha[4]=38.216677;
alpha[1] = 0.298073; alpha[2]= 1.242567; alpha[3]= 5.782948; alpha[4] = 38.474970;
c[1] = 1.0; c[2] = 1.0; c[3] = 1.0; c[4] = 1.0; 
	
for(i=1;i<=4;i++)
{
	for(j=i;j<=4;j++)
	{
		S[i][j] = pow(pi/(alpha[i]+alpha[j]),1.5);
		S[j][i] = S[i][j];
	}
}
for(i=1;i<=4;i++)
{
	for(j=i;j<=4;j++)
	{
		h[i][j] = 3.0*alpha[i]*alpha[j]*pow(pi,1.5)/(pow((alpha[i]+alpha[j]),2.5)) - 4.0*pi/(alpha[i]+alpha[j]); ;
		h[j][i] = h[i][j];
	}
}
normalize(S,c);
while(dif>eps){
for(i=1;i<=4;i++){
	for(j=1;j<=4;j++){
		F[i][j]=h[i][j];}}

for(i=1;i<=4;i++){
	for(j=1;j<=4;j++){
		for(k=1;k<=4;k++){
			for(l=1;l<=4;l++){
				qpqrs = 2*pow(pi,2.5)/((alpha[i] + alpha[k])*(alpha[j] + alpha[l])*sqrt(alpha[i]+alpha[j]+alpha[k]+alpha[l])) ;
				F[i][k] += qpqrs*c[j]*c[l];
			}
		}
	}
}

double *d;
d = dvector(1,4);
eigen(F, S, d, eigv, 4);
E_0 = d[1];
free_dvector(d,1,4);
normalize(S,eigv);
dif = difference(c,eigv);
for(i=1;i<=4;i++)
	c[i]=eigv[i];

}

sum = 0.0; double sum1 = 0.0;
for(i=1;i<=4;i++){
	for(j=1;j<=4;j++){
		for(k=1;k<=4;k++){
			for(l=1;l<=4;l++){
				qpqrs = 2*pow(pi,2.5)/((alpha[i] + alpha[k])*(alpha[j] + alpha[l])*sqrt(alpha[i]+alpha[j]+alpha[k]+alpha[l])) ;
				sum += qpqrs*c[i]*c[j]*c[k]*c[l];
			}
		}
	}
}
for(i=1;i<=4;i++){
	for(j=1;j<=4;j++){
		sum1 += c[i]*c[j]*h[i][j];}} 
printf("Ground state energy = %.8f\n",2*sum1+sum);

free_dmatrix(F,1,4,1,4);
free_dmatrix(S,1,4,1,4);
free_dmatrix(h,1,4,1,4);
free_dvector(alpha,1,4);
free_dvector(eigv,1,4);
free_dvector(c,1,4);
return 0;
}
#undef NRANSI
