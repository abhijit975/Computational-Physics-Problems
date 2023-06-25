#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 2
#define pi acos(-1.0)

// defining parameters............................................................
float g = 0.0, a = 1.0, b = 1.0, F0 = 0.0, w = 0.0;
//................................................................................

float dxsav,*xp,**yp;  
int kmax,kount;
int nrhs; 

void derivs(float x,float y[],float dydx[])
{
	nrhs++;
	dydx[1] = y[2];
	dydx[2] = -g*y[2] + a*y[1] - b*y[1]*y[1]*y[1] + F0*cos(w*x);

}

int main(void)
{
	printf("Enter value of gamma = \n");
	scanf("%f",&g);
	printf("Enter value of F0 = \n");
	scanf("%f",&F0);
	printf("Enter value of omega = \n");
	scanf("%f",&w);
	FILE *fp;
	fp = fopen("duffing.dat","w");
	int i,nbad,nok;
	float eps=1.0e-4,h1,hmin,x1,x2,*ystart;
	ystart=vector(1,N);
	ystart[1] = 0.0;
	ystart[2] = 0.0;	

	for(i=1;i<25000;i++)
	{	
		h1 = 0.1;
		hmin = 0.0;
		x1 = (i-1)*2*pi/w; 
		x2 = i*2*pi/w;      	

		xp=vector(1,200);
		yp=matrix(1,10,1,200);

		nrhs=0;
		kmax=100;
		dxsav=(x2-x1)/20.0;
		odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs);

		fprintf(fp,"%lf\t%lf\n",yp[1][kount],yp[2][kount]);

		ystart[1] = yp[1][kount];
		ystart[2] = yp[2][kount];
		free_matrix(yp,1,10,1,200);
		free_vector(xp,1,200);
		
	}
	printf("pi = %lf\n",pi);
	free_vector(ystart,1,N);
	fclose(fp);
	return 0;
}
#undef NRANSI
