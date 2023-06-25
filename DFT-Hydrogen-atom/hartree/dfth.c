/* Driver for routine rtsec */

#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 100
#define NBMAX 20
#define X1 -1.5
#define X2 -0.2
#define NVAR 2
#define NSTEP 1000
#define pi acos(-1)


// Solving the Kohn-Sham differential equation without the potential terms


float E;
float Z = 2.0;
float *vh;

void derivs(float x, float yb[], float dydx[])
{
	int i;
	float dstep = (10 - 0.001)/NSTEP;
	i = (int)(x/dstep) + 1;
	dydx[1] = yb[2];
	dydx[2] = - 2.0*Z*yb[1]/x - 2.0*E*yb[1] + 2.0*vh[i]*yb[1];
}

extern float *xx, **y;

float fx(float temp)
{
	E = temp; FILE *fp;
 	float x2 = 0.001, x1 = 10.0, *vstart;
	float value;
    	vstart=vector(1,NVAR);
    	vstart[1]= 0.0; 
    	vstart[2]= -0.001;
	xx=vector(1,NSTEP+1);
	y=matrix(1,NVAR,1,NSTEP+1);
  	rkdumb(vstart,NVAR,x1,x2,NSTEP,derivs);
    	value = y[1][NSTEP+1];
	free_matrix(y,1,NVAR,1,NSTEP+1);
	free_vector(xx,1,NSTEP+1);
	free_vector(vstart,1,NVAR);
    return value;
}

float u(float ua[])
{
	float x2 = 0.001, x1 = 10.0, *vstart;
	vstart=vector(1,NVAR);
    	vstart[1]= 0.0; 
    	vstart[2]= -0.001;
	xx=vector(1,NSTEP+1);
	y=matrix(1,NVAR,1,NSTEP+1);
  	rkdumb(vstart,NVAR,x1,x2,NSTEP,derivs);
	for(int i=1;i<=NSTEP+1;i++)
		ua[i]=y[1][i];
	float sum = 0.0;
	float dstep = (10.0-0.001)/NSTEP;
	for(int i=1;i<=NSTEP;i++)
		sum += (ua[i]*ua[i] + ua[i+1]*ua[i+1])*0.5*dstep;
	for(int i=1;i<=NSTEP+1;i++)
		ua[i] = ua[i]/sqrt(sum);
	sum = 0.0;
	for(int i=1;i<=NSTEP;i++)
		sum += (ua[i]*ua[i] + ua[i+1]*ua[i+1])*0.5*dstep;
	return sum;
	free_matrix(y,1,NVAR,1,NSTEP+1);
	free_vector(xx,1,NSTEP+1);
	free_vector(vstart,1,NVAR);
}

float *uarray, *vh;
float rmax = 10.0; float rmin = 0.001;

void derivsU(float x, float yb[], float dydx[])
{
	int i;
	float dstep = (10 - 0.001)/NSTEP;
	i = (int)((10.0-x)/dstep) + 1;
	dydx[1] = yb[2];
	dydx[2] = -uarray[i]*uarray[i]/x;
}

int main(void)
{
	int i,nb=NBMAX;
	float xacc,root,*xb1,*xb2; float r,alpha;
	int j;

	xb1=vector(1,NBMAX);
	xb2=vector(1,NBMAX);
	zbrak(fx,X1,X2,N,xb1,xb2,&nb);
	for (i=1;i<=nb;i++) {
		xacc=1.0e-3;
		root=rtsec(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %14.6f\n",i,root,fx(root));
	}
	E = root;
	free_vector(xb2,1,NBMAX);
	free_vector(xb1,1,NBMAX);
	float dstep = (0.001-10.0)/NSTEP;
	uarray = vector(1,NSTEP+1);
	float sum = u(uarray);

	float x1=rmin,x2=rmax,*vstart;
	vstart=vector(1,NVAR);
	// Note: The arrays xx and y must have indices up to NSTEP+1 
	vh=vector(1,NSTEP+1);
	xx=vector(1,NSTEP+1);
	y=matrix(1,NVAR,1,NSTEP+1);
	vstart[1] = 0.0;
	vstart[2] = 1.0;
	
	rkdumb(vstart,NVAR,x1,x2,NSTEP,derivsU);
	alpha = (1.0 - y[1][NSTEP+1])/rmax;
	for (j=1;j<=NSTEP+1;j++)
		vh[j]= (y[1][j]+alpha*xx[j])/xx[j];
	FILE *fp;
	fp = fopen("data.dat","w");
	for(int i=1;i<=NSTEP+1;i++)
		fprintf(fp,"%f\t%f\n", 0.001 - (i-1)*dstep, vh[i]);
	//printf("sum = %f\n",sum);
	fclose(fp);
	free_vector(uarray,1,NSTEP+1);
	free_matrix(y,1,NVAR,1,NSTEP+1);
	free_vector(xx,1,NSTEP+1);
	free_vector(vstart,1,NVAR);
	free_vector(vh,1,NSTEP+1);
	return 0;
}
#undef NRANSI
