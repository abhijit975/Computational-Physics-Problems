#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

int main()
{
time_t begin = time(NULL);
FILE *fp,*fptr;
fp = fopen("duffing.dat","r");
float x,p;
int i=0, j, k, len=0;
while(feof(fp)==0)
	{
	fscanf(fp,"%f\t%f",&x,&p);
	len += 1;
	}
fclose(fp);
fp = fopen("duffing.dat","r");
fptr = fopen("dimension.dat","w");
float *position = (float *)malloc(len*sizeof(float));
float *velocity = (float *)malloc(len*sizeof(float));
while(feof(fp)==0)
	{
	fscanf(fp,"%f\t%f",&position[i],&velocity[i]);
	i += 1;
	}

int l = 2;
for(int iter=1;iter<8;iter++)
{
float x2 = 1.5, x1 = -1.5, y2 = 1.5, y1 = -1.5;
l = l*2;
float dx = (x2-x1)/l;
float dy = (y2-y1)/l;
int **a = (int **)malloc(l*sizeof(int *));
for(i=0;i<l;i++)
	a[i] = (int *)malloc(l* sizeof(int));

for(i=0;i<l;i++)
	for(j=0;j<l;j++)
		a[i][j]=0;

float intx1=x1,intx2,inty1=y1,inty2;

for(i=0;i<l;i++)
	{intx2 = intx1;
	intx1 = x1+(i+1)*dx;
	for(j=0;j<l;j++)
		{inty2 = inty1;
		inty1 = y1+(j+1)*dy;
		for(k=0;k<len;k++)
			{if(position[k]<intx1 && position[k]>=intx2 && velocity[k]<inty1 && velocity[k]>=inty2)
				{a[i][j]=1;
				break;}
			}
		}
	}

int count = 0;
for(i=0;i<l;i++)
	for(j=0;j<l;j++)
		if(a[i][j]!=0)
			count += 1;
free(a);
float logcount = log(count);
float logb = log((x2-x1)/l);
fprintf(fptr,"%f\t%f\n",logb,logcount);
}
fclose(fptr);
fclose(fp);
free(position);
free(velocity);

time_t end = time(NULL);
printf("The time elapsed = %ld seconds\n",end-begin);
}
