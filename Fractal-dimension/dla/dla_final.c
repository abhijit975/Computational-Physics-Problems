#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

int L, ll;
long seed;
int **lattice;			// a pointer to the pointer to lattice

void initialize_lattice();
void addparticle();


int main(){

  int n, npart;
  FILE *fp,*fptr;
  fp = fopen("data.dat","w");
  fptr = fopen("gyration.dat","w");
  float meanx,meany,sqrdx,sqrdy,mean,meansqrd,radgyr;
  printf("The lattice will be L x L. How large is L? (Should be odd integer.): ");
  scanf("%d",&L);
  printf("Input number of particles added to the seed: ");
  scanf("%d",&npart);
  printf("Input random number seed. (Should be large negative odd integer.): ");
  scanf("%ld",&seed);

  ll = (L-1)/2;
  lattice = imatrix(-ll,ll,-ll,ll);	// declares lattice as matrix

  initialize_lattice();		// initializes empty lattice with seed in middle

  for (n=1; n<=npart;++n) {
    addparticle();		// adds new particle to seed
	if(n%200==0)		// Calculates radius of gyration after adding every 200 particles 
		{
		meanx=0.0;meany=0.0;sqrdx=0.0,sqrdy=0.0;
		 for (int y=-ll; y<=ll; y++){
    			for (int x=-ll; x<=ll; x++){
      				if (lattice[y][x] != 0)  {
					meanx += x;
					meany += y;
					sqrdx += x*x;
					sqrdy += y*y;
				}
    			}
  		}
		mean = (meanx/n)*(meanx/n) + (meany/n)*(meany/n);
		meansqrd = (sqrdx + sqrdy)/n;
		radgyr = sqrt(meansqrd - mean);
		fprintf(fptr,"%f\t%f\n",log(radgyr),log(n)); //prints to file log(rg) and log(n)
		}
  }
 for (int y=-ll; y<=ll; y++){		//prints to file tha data points x and y for plotting dla.
    for (int x=-ll; x<=ll; x++){
      if (lattice[y][x] != 0)  fprintf(fp,"%d\t%d\n",x,y);
    }
  }
  free_imatrix(lattice,-ll,ll,-ll,ll);	// frees memory
  fclose(fp);
  fclose(fptr);
  return 0;
}

void addparticle()
{
  int x,y,roll,side,site,done,nsum;

  done = 0;


     roll = 4*(L-1)*ran2(&seed); // random int to determine walker's initial location: (x,y)
     side = roll/(L-1);
     site = roll%(L-1);
     if (side==0) {x=-ll+site; y=ll;}
     else if (side==1) {x=ll; y=ll-site;}
     else if (side==2) {x=ll-site; y=-ll;}
     else if (side==3) {x=-ll; y=-ll+site;}

	while(1>0)
		{
		if(lattice[x][y]==0)
			break;
		else
		{
     		roll = 4*(L-1)*ran2(&seed); // random int to determine walker's initial location: (x,y)
     		side = roll/(L-1);
     		site = roll%(L-1);
     		if (side==0) {x=-ll+site; y=ll;}
     		else if (side==1) {x=ll; y=ll-site;}
     		else if (side==2) {x=ll-site; y=-ll;}
     		else if (side==3) {x=-ll; y=-ll+site;}
		}}

     while (done == 0)  // do unless done or dead   
     {
	if((abs(x)<ll)&&(abs(y)<ll)) 
		{
		roll = 4*ran2(&seed);		// random int to determine walk direction
       		if (roll==0) --x;
       		else if (roll==1) ++x;
       		else if (roll==2) --y;
       		else if (roll==3) ++y;		
		}
	else if(abs(x)==ll || abs(y)==ll)     //what happens at the boundary? either 2 neighbor or 3 neighbor.
		{
		if(abs(x)==abs(y))
			{
			roll = 2*ran2(&seed);
			if(roll==0) 
				{
				if(x>0){x = x-1;}
				else{x=x+1;}
				}
			else if(roll==1) 
				{
				if(y>0){y = y-1;} 
				else{y=y+1;}
				}
			}
		else
		    {
		    roll = 3*ran2(&seed);
		    if (x==-ll) {
			if(roll==0) x = x+1;
			else if (roll==1) y = y+1;
			else if (roll==2) y = y-1;
			}
		   else if (x== ll) {
			if(roll==0) x = x-1;
			else if (roll==1) y = y+1;
			else if (roll==2) y = y-1;
			}
		   else if (y== ll) {
			if(roll==0) y = y-1;
			else if (roll==1) x = x+1;
			else if (roll==2) x = x-1;
			}
		   else if (y==-ll) {
			if(roll==0) y = y+1;
			else if (roll==1) x = x+1;
			else if (roll==2) x = x-1;
			}
		    }
		}
	if(abs(x)<ll && abs(y) < ll)
		nsum = lattice[y-1][x]+lattice[y+1][x]+lattice[y][x-1]+lattice[y][x+1];
 	else if(y==-ll && abs(x)<ll)		
		nsum = lattice[y+1][x] + lattice[y][x+1] + lattice[y][x-1];
	else if (y==ll && abs(x)<ll)
		nsum = lattice[y-1][x] + lattice[y][x+1] + lattice[y][x-1];
	else if(x==-ll && abs(y)<ll) 
		nsum = lattice[y+1][x] + lattice[y-1][x] + lattice[y][x+1];     
	else if (x==ll && abs(y)<ll)
		nsum = lattice[y+1][x] + lattice[y-1][x] + lattice[y][x-1];
	else if(x==-ll && y==-ll)
		nsum = lattice[y+1][x] + lattice[y][x+1];
	else if(x==ll && y==-ll)
		nsum = lattice[y+1][x] + lattice[y][x-1];
	else if(x==-ll && y==ll)
		nsum = lattice[y-1][x] + lattice[y][x+1];
	else if(x==ll && y==ll)
		nsum = lattice[y-1][x] + lattice[y][x-1];
				
       if (nsum > 0) done=1;		// determine if next to seed
     } 

  lattice[y][x] = 1;		// add walker to seed

  return;
}

void initialize_lattice()
{
  int x,y;

  for (y=-ll; y<=ll; ++y){
    for (x=-ll; x<=ll; ++x){
      lattice[y][x] = 0;		// initialize empty lattice
    }
  }
  lattice[0][0] = 1;		// start seed as single particle at center of lattice
  return;
}

#undef NRANSI
