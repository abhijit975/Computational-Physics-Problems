#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 2
#define pi acos(-1.0)

//defining global constants..........................................
float E;
float e0 = 5.9, rho = 3.57;
float dxsav,*xp,**yp;  
int kmax,kount;
int nrhs,l; 
float con = sqrt(5.9*6.12/25.0);
float rmax = 35.7, rmin = 1.785;
//...................................................................
// Differential equation.............................................
void derivs(float x,float y[],float dydx[])
{
	nrhs++;
	float c = 6.12/(rho*rho);
	dydx[1] = y[2];
	if(x<=rmax){	
		dydx[2] = (c*(e0*(pow((rho/x),12)-2.0*pow((rho/x),6))-E) + l*(l+1)/(x*x)) * y[1];}      
	else{
		dydx[2] = (c*(-E) + 1.0*l*(l+1)/(x*x)) * y[1];}
}
//...................................................................
int main(void)
{
	int lmax; 
	printf("Enter max value of l = \n"); //asks for user input for maximum value of l.
	scanf("%d",&lmax);
	FILE *fp;
	fp = fopen("test.dat","w");    //to write energy and cross-section in file
	int j;
	float temp, K, sindl, tandl, u1, u2, r1, r2, wave, sum, cross;
	for(j=1; j<101; j++)
		{
		E = 0.3+j*0.08;         // Energy changing inside the loop, starting from 0.3.
		wave = sqrt(6.12*E/(rho*rho));  //wave vector
		sum = 0.0;
		for(l=0; l<lmax+1; l++)
			{
			int i,nbad,nok;
			float eps=1.0e-8,h1,hmin,x1,x2,*ystart;
			h1 = 0.01;
			hmin = 0.0;
			x1 = rmin;		//integration lower limit with LJ potential
			x2 = rmax;		//integration upper limit with LJ potential
			temp = pow(x1,5);
			ystart=vector(1,N);
			xp=vector(1,200);
			yp=matrix(1,10,1,200);
			ystart[1] = exp(-con/temp);  //initial condition u(r)
			ystart[2] = 5*con*ystart[1]/pow(rmin,6);  //initial condition u'(r)
			nrhs=0;
			kmax=1000;
			dxsav=(x2-x1)/100.0;
			odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs); //Solving ODE 
			u1 = yp[1][kount];	// u1 and r1.
			r1 = rmax;
			ystart[1] = u1;   	//updating initial conditions for next integrations
			ystart[2] = yp[2][kount]; 
			free_matrix(yp,1,10,1,200);   // freeing xp and yp and redefining them
			free_vector(xp,1,200);
			xp=vector(1,200);
			yp=matrix(1,10,1,200);	
			nrhs=0;
			kmax=1000;h1 = 0.01;
			hmin = 0.0;
			x1 = r1;     		// integration limits without LJ potential
			x2 = r1 + pi/wave;
			dxsav=(x2-x1)/100.0;
			odeint(ystart,N,x1,x2,eps,h1,hmin,&nok,&nbad,derivs,rkqs); //Solving ODE without LJ potential
			u2 = yp[1][kount];  // storing u2 and r2.
			r2 = xp[kount];
			free_matrix(yp,1,10,1,200);
			free_vector(xp,1,200);
			free_vector(ystart,1,N);
			K = r1*u2/(r2*u1); 	// Finding K
			float sj1, sy1, sj2, sy2, sjp1, syp1, sjp2, syp2;
			sphbes(l,wave*r1,&sj1,&sy1,&sjp1,&syp1); // finding jl(kr1) and nl(kr1)
			sphbes(l,wave*r2,&sj2,&sy2,&sjp2,&syp2); // finding jl(kr2) and nl(kr2)
			tandl = (K*sj1 - sj2)/(K*sy1 - sy2);   // calculating tan(delta_l)
			sindl = (tandl*tandl)/(1.0+tandl*tandl);  //calculating sin(delta_l)
			sum = sum + sindl*(2*l+1); 	// calculating sum over l for one E.
			}
		cross = 4*pi*sum/(wave*wave*rho*rho);  // calculating cross-section for particular E
		fprintf(fp,"%f\t%f\n",E,cross);
		}
	fclose(fp);
	return 0;
}
#undef NRANSI	
