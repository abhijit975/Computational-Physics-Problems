#include <stdio.h>
#include<stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include <math.h>
#include <time.h>
#define pi acos(-1.0)

static int long seed;
double rcut = 2.5;
int np;
double Q;

double sign(double x)
{
double value;
if (x>0)
	value = -1.0;
else if (x<0)
	value = 1.0;
return value;
}

double potential(double x)
{
double value;
value = 48.0/pow(x,14.0) - 24.0/pow(x,8.0);
return value;
}

void forcecal(double **xdist, double **ydist, double **zdist, double **rdist, double **force, double **x, double l)
{
int i,j;
double value;
for(i=1;i<=np;i++){
	for(j=1;j<=np;j++){
		xdist[i][j] = x[i][1] - x[j][1];
		ydist[i][j] = x[i][2] - x[j][2];
		zdist[i][j] = x[i][3] - x[j][3];
		if (fabs(xdist[i][j])>l/2.0)
			xdist[i][j] = sign(xdist[i][j])*(l-fabs(xdist[i][j]));
		if (fabs(ydist[i][j]) >l/2.0)
			ydist[i][j] = sign(ydist[i][j])*(l-fabs(ydist[i][j]));
		if(fabs(zdist[i][j])>l/2.0)
			zdist[i][j] = sign(zdist[i][j])*(l-fabs(zdist[i][j]));
			
		if(i<j){
			rdist[i][j] = sqrt(xdist[i][j]*xdist[i][j] + ydist[i][j]*ydist[i][j] + zdist[i][j]*zdist[i][j]);
			rdist[j][i] = rdist[i][j];
			}
		}
	rdist[i][i]=0.0;
	}
for(i=1;i<=np;i++)
	for(j=1;j<=3;j++)
		force[i][j] = 0.0;
for(i=1;i<=np;i++){
	for(j=1;j<=np;j++){
		if (i!=j){
			value = potential(rdist[i][j]);
			force[i][1] += value*xdist[i][j];
			force[i][2] += value*ydist[i][j];
			force[i][3] += value*zdist[i][j];
			}
		}
	}
}

void pbc(double **x, double l)
{
int i,j;
double temp;
for(i=1;i<=np;i++){
	for(j=1;j<=3;j++){
		if(x[i][j]>l){
			temp = fmod(temp,l);
			x[i][j] = temp;}
		else if (x[i][j]<0){
			temp = fabs(x[i][j]);
			temp = fmod(temp,l);
			x[i][j] = l-temp;}
		}
	}
}

void momentascaleinit(double **p, double T)
{
int i;
double pxsum = 0.0, pysum = 0.0, pzsum = 0.0;
for(i=1;i<=np;i++){
	pxsum += p[i][1];
	pysum += p[i][2];
	pzsum += p[i][3];
	}

for(i=1;i<=np;i++){
	p[i][1] -= pxsum/np; 
	p[i][2] -= pysum/np;
	p[i][3] -= pzsum/np;
	}

double sumvsq = 0.0, scale;
for(i=1;i<=np;i++)
	sumvsq += p[i][1]*p[i][1] + p[i][2]*p[i][2] + p[i][3]*p[i][3];
scale = sqrt(3.0*np*T/sumvsq);

for(i=1;i<=np;i++){
	p[i][1] *= scale; 
	p[i][2] *= scale;
	p[i][3] *= scale;
	}
}

void momentascale(double **p, double T)
{
int i;
double sumvsq = 0.0, scale;
for(i=1;i<=np;i++)
	sumvsq += p[i][1]*p[i][1] + p[i][2]*p[i][2] + p[i][3]*p[i][3];
scale = sqrt(3.0*np*T/sumvsq);

for(i=1;i<=np;i++){
	p[i][1] *= scale; 
	p[i][2] *= scale;
	p[i][3] *= scale;

	}
}

double verletnh (double **xdist, double **ydist, double **zdist, double **rdist, double **force, double **x, double **p, double l, double zeta, double T)
{
int i,j;
double h = 0.004;
double sumv = 0.0,zetatemp;
for(i=1;i<=np;i++){
	for(j=1;j<=3;j++){
		sumv += p[i][j]*p[i][j];}}
zetatemp = zeta + 0.25*h*(sumv - (3.0*np+1.0)*T)/Q;
for(i=1;i<=np;i++){
	for(j=1;j<=3;j++){
		p[i][j] = p[i][j] + 0.5*h*(force[i][j] - zeta*p[i][j]);
		x[i][j] = x[i][j] + h*p[i][j]; 
	}
}
pbc(x,l);
forcecal(xdist,ydist,zdist,rdist,force,x,l);
sumv = 0.0;
for(i=1;i<=np;i++){
	for(j=1;j<=3;j++){
		sumv += p[i][j]*p[i][j];}}
zeta = zetatemp + 0.25*h*(sumv - (3.0*np+1.0)*T)/Q;
//zeta = 0.0;
for(i=1;i<=np;i++)
	for(j=1;j<=3;j++)
		p[i][j] = (p[i][j] + force[i][j]*h*0.5)/(1.0+0.5*h*zeta);
return zeta;
}

double temperature(double **p)
{
int i;
double sumvsq = 0.0,temp;
for(i=1;i<=np;i++)
	sumvsq += p[i][1]*p[i][1] + p[i][2]*p[i][2] + p[i][3]*p[i][3];
temp = sumvsq/(3.0*np);
return temp;
}

int main()
{
seed = time(NULL);
double **p, **x;
FILE *fp;
int M;
int i,j;
double T, l, rho;
printf("Enter the value of M:\n");
scanf("%d",&M);
printf("Enter the desired temperature:\n");
scanf("%lf",&T);
printf("Enter number density = \n");
scanf("%lf",&rho);
printf("Enter the value for Q:\n");
scanf("%lf",&Q);
np = (int)4*pow(M,3);
l = pow(np*1.0/rho, 0.333);
printf("Length of box = %lf\n",l);
//initialising momenta................................................................................
p = dmatrix(1,np,1,3); 
for(i=1;i<=np;i++)
	for(j=1;j<=3;j++)
		p[i][j] = gasdev(&seed);
//....................................................................................................
//scaling momenta to get fixed temperature............................................................
momentascaleinit(p,T);
//....................................................................................................
//initialising position...............................................................................
x = dmatrix(1,np,1,3);
int ix,iy,iz,nc,index,m;
double cell,cell2;
nc = (int)(pow(np/4,1./3.) + 0.1);
cell = 1./(double)nc;
cell2 = 0.5*cell;
x[1][1] = 0.0;
x[2][1] = cell2;
x[3][1] = 0.0;
x[4][1] = cell2;
x[1][2] = 0.0;
x[2][2] = cell2;
x[3][2] = cell2;
x[4][2] = 0.0;
x[1][3] = 0.0;
x[2][3] = 0.0;
x[3][3] = cell2;
x[4][3] = cell2;
m=0;
printf("nc = %d\t cell = %lf\t cell2 = %lf\n",nc, cell, cell2);
for (iz=1; iz<=nc; iz++) {
	for (iy=1; iy<=nc; iy++) {
		for (ix=1; ix<=nc; ix++){
			for (index=1; index<=4;index++){
				x[index+m][1] = x[index][1] + cell * 1.0* (ix-1);	
				x[index+m][2] = x[index][2] + cell * 1.0* (iy-1);
				x[index+m][3] = x[index][3] + cell * 1.0* (iz-1);
			}
		m +=4;
		}
	}
}
for(i=1;i<=np;i++){
	x[i][1] *= l;
	x[i][2] *= l;
	x[i][3] *= l;
}
//....................................................................................................
//calculating force for each particle.................................................................
double **xdist, **ydist, **zdist, **rdist;
xdist = dmatrix(1,np,1,np);
ydist = dmatrix(1,np,1,np);
zdist = dmatrix(1,np,1,np);
rdist = dmatrix(1,np,1,np);
double **force;
double value, h = 0.004;
force = dmatrix(1,np,1,3);
forcecal(xdist,ydist,zdist,rdist,force,x,l);

//equilibration..................................................................................
int iter;
double zeta = 0.0;
for(iter=1; iter<=50000; iter++){
zeta = verletnh(xdist,ydist,zdist,rdist,force,x,p,l,zeta,T);
if(iter%10==0)
	momentascale(p,T);
}
fp = fopen("Temperature2.dat","w");
for(iter=1; iter<=100; iter++){
zeta = verletnh(xdist,ydist,zdist,rdist,force,x,p,l,zeta,T);
fprintf(fp,"%lf\n",temperature(p));
//if(iter%10==0)
//	momentascale(p,T);
}
fclose(fp);
free_dmatrix(x,1,np,1,3);
free_dmatrix(p,1,np,1,3);
free_dmatrix(force,1,np,1,3);
free_dmatrix(xdist,1,np,1,np);
free_dmatrix(ydist,1,np,1,np);
free_dmatrix(zdist,1,np,1,np);
free_dmatrix(rdist,1,np,1,np);
return 0;
}
#undef NRANSI
