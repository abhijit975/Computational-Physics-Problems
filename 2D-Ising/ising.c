#include <stdio.h>
#include<stdlib.h>
#include<math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"
#include<time.h>

static long int seed;

float initialize(float **x)
{
int i,j;
float m = 0.0;
for(i=1;i<=40;i++){
	for(j=1;j<=40;j++){
		if(ran2(&seed)<0.5){
			x[i][j] = 1.0;
			m += x[i][j];}
		else{
			x[i][j] = -1.0;
			m += x[i][j];}
		}
	}
return m;
}

float spinflip(float **s, float expval[], float m)
{
int i,j;
float x;
int xm,xp,ym,yp;
//float m=0.0;
for(i=1;i<=40;i++){
		xp = i+1;
		xm = i-1;
		if(i==1)
			xm = 40;
		else if(i==40)
			xp = 1;
	for(j=1;j<=40;j++){
		yp = j+1;
		ym = j-1;
		if(j==1)
			ym = 40;
		else if(j==40)
			yp = 1;
		
		x = s[i][j]*(s[i][yp]+s[i][ym]+s[xp][j]+s[xm][j]);
		if(x>=0 && ran2(&seed)<expval[3-(int)(x/2.0)]){
			s[i][j]= -s[i][j];
			m += 2.0*s[i][j];
			}
		else if(x<0){
			s[i][j] = -s[i][j];
			m += 2.0*s[i][j];
			}
		}
	}

return m;
}

int main()
{
seed = time(NULL);
clock_t start, end;
start = clock();
FILE *fp;
float **s;
float T,M;
s = matrix(1,40,1,40);

float *expval;
expval = vector(1,5);
int i,j,itemp,iter;
fp = fopen("mag4.dat","w");
for(itemp=1;itemp<=5;itemp++){
T = 3.5 + (itemp-1)*0.1;
for(i=1;i<=5;i++){
	expval[i] = exp(-(8.0-(i-1)*4.0)/T);}

M = initialize(s);
printf("Magnetization = %f\n",M/1600.0);
for(iter =1;iter<11000000;iter++){
M = spinflip(s,expval,M);
fprintf(fp,"%d\t%f\n",iter,M/1600.0);}
}
fclose(fp);
free_matrix(s,1,40,1,40);
free_vector(expval,1,5);
end = clock();
printf("Time taken = %f\n",((float) (end - start)) / CLOCKS_PER_SEC);
return 0;
}
#undef NRANSI
