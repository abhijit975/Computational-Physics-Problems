#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

float correlation(int t, int tot, float M[])
{
int iter;
float sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
float corr;
for(iter=1;iter<=tot-t;iter++){
	sum1 += M[iter]*M[iter+t];
	sum2 += M[iter];
	sum3 += M[t+iter];
	}
 
corr = (sum1*(tot-t)*1.0 - sum2*sum3)/(1.0*(tot-t)*(tot-t));
return corr;
}

int main()
{
FILE *fp,*fptr;
fp = fopen("mag2.dat","r");
fptr = fopen("correlation5.dat","w");
float *M,x,corr,corr1,tol;
tol = 1.0e-8;
int i,j,tot;
float sum, corrprev, T = 1.5;
tot = 1000000;
int n = 11000000;
int idum,datapts = 11000000*10;
M = vector(1,tot); 
/*for(idum = 1; idum<=36;idum++){
	sum = 0.0;
	corrprev = 0.0;
for(i=1;i<=datapts;i++){
	fscanf(fp,"%f\n",&x);
	if(i>tot && i<=2*tot)
		M[i-tot-1] = x;
	}
for(i=1;i<=tot;i++){
	corr = correlation(i,tot,M);
	if(corr<tol)
		break;
	if(i==1)
		corr1 = corr;
	sum += corrprev + corr/corr1;
	corrprev = corr/corr1;
	fprintf(fptr,"%f\t%d\t%f\t%f\n",T,i,corr/corr1,log(corr/corr1));
	}
	//fprintf(fptr, "%f\t%f\n", T,(sum-1.0)/2.0);
	T = T + 0.1;
}*/
int count = 0;
for(idum = 1;idum<=datapts;idum++)
{
	fscanf(fp,"%d\t%f\n",&j,&x);
	if(idum%n>11000000-tot){
		count += 1;
		M[count] = x;}
	if(idum%n==0){
		count = 0;
		for(i=1;i<=tot/100;i++){
			corr = correlation(i,tot,M);
			if(i==1)
				corr1 = corr;
			if(corr<tol){
				fprintf(fptr,"%f\t%d\t%f\n",T,i,0.0);
				break;
			}
			fprintf(fptr,"%f\t%d\t%f\n",T,i,corr/corr1);			
		}
		T = T + 0.1;
	}
}	
free_vector(M,1,1000000);
fclose(fp);
fclose(fptr);
return 0;
}
#undef NRANSI
