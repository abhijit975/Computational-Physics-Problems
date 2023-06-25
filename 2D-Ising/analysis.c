#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

int main()
{
float m,div,T=0.5;
div = 10000000.0;
int i,j,iter,total,ct=0;
float avgm = 0.0, avgd = 0.0;
FILE *fp,*fp1,*fp2;
fp = fopen("mag.dat","r");
fp1 = fopen("avgm.dat","w");
fp2 = fopen("avgmsq.dat","w");
float *averagem,*averagemsq;
averagem = vector(1,36);
averagemsq = vector(1,36);
total = 11000000*36;
for(iter = 1; iter<= total; iter++){
	fscanf(fp,"%d\t%f\n",&i,&m);
	if(iter%11000000>1000000){
		avgm += m;
		avgd += m*m;}
	else if(iter%11000000==0){
		avgm += m;
		avgd += m*m;
		ct++;
		averagem[ct] = avgm/div;
		averagemsq[ct] = avgd/div;
		avgm = 0.0;
		avgd = 0.0;}
	}
fseek(fp,0,SEEK_SET);
ct =1;
float sum1 = 0.0, sum2 = 0.0;
for(iter=1;iter<=total;iter++){
	fscanf(fp,"%d\t%f\n",&i,&m);
	if(iter%11000000 > 1000000){
        	sum1 +=(m-averagem[ct])*(m-averagem[ct]);
        	sum2 +=(m*m-averagemsq[ct])*(m*m-averagemsq[ct]);}
	if(iter%11000000 == 0){
    		fprintf(fp1, "%f\t%f\t%f\n",T,averagem[ct],sqrt(sum1/div));
    		fprintf(fp2, "%f\t%f\t%f\n",T,averagemsq[ct],sqrt(sum2/div));
    		T = T + 0.1;
    		ct += 1;
    		sum1 = 0.0;
    		sum2 = 0.0;
    	    }
     }
fclose(fp);
fclose(fp1);
fclose(fp2);
free_vector(averagem,1,36);
free_vector(averagemsq,1,36);	

return 0;
}
#undef NRANSI
