//Calcula el numero de voids a partir de las posiciones
//usando spherical-overdensity>200
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"

const float rho_cr=2.77526627;		//	[10¹¹ h² Ms / Mpc³]
const float pi=4.*atan(1.);		//	numero pi



int main(int argc, char **argv){



  char inFile[400];
  char outFile[400];
  char outFile2[400];
  sprintf(inFile,"%s",halos_file);
  sprintf(outFile,"dnh_%s",halos_file);
  sprintf(outFile2,"nh_%s",halos_file);


  int i, j, k, l;
  int NumHalos=0;
  int id;

  float x,y,z,masa,radio;
  float min,max;

  float nh[bin]={0.};
  float Masa[bin]={0.};


///////////////////////////////////		masa 	maxima y minima

  printf("Leyendo el catalogo para establecer radio max y min\n");fflush(stdout);
  FILE * L;
  FILE * E;
  FILE * F;

  L=fopen(inFile,"r");
  max=0.;
  min=1000000.;
  while(fscanf(L,"%f %f %f %f %f\n",&x, &y, &z, &masa, &radio)!=EOF){
    if(masa<min)
      min=masa;
    if(masa>max)
      max=masa;
    NumHalos++;}
  fclose(L);
  printf("Se encontraron %d halos en el catálogo\n",NumHalos);fflush(stdout);


printf("%e %e\n",max,min);
max=1.89e+04;
min=3.8;

///////////////////////		definicion de los bines radiales
  float Mmax =max*1.01;
  float m    =min*0.99;
  float c    =pow(Mmax/m,1./bin);	//	constante de proporcionalidad log
  float aux1 =1./log(c);
  float aux2;



//////////////////////		lleno los bines, contador y radio
  L=fopen(inFile,"r");
  for(i=0;i<NumHalos;i++){
    fscanf(L,"%f %f %f %f %f\n",&x, &y, &z, &masa, &radio);
    aux2=aux1*log(masa/m);
    id=floor(aux2);		//	es el bin radial de este void
    nh[id]+=1.;			//	contador entero del número de voids
    Masa[id]+=masa;}
  fclose(L);





//////////////////////////		nh
  F=fopen(outFile2,"w+");
  aux1=m*pow(c,1*(bin-1));
  aux2=nh[bin-1];
  fprintf(F,"%f %e\n",aux1,aux2*pow(Lbox,-3.));
  for(i=bin-2;i>=0;i--){
    aux1=m*pow(c,i*1.);
    aux2+=nh[i];
    fprintf(F,"%f %e\n",aux1*100000000000,aux2*pow(Lbox,-3.));
  }
  fclose(F);



//////////////////////////		dnh
  E=fopen(outFile,"w+");
  float factor;
  for(i=0;i<bin;i++){
    Masa[i]/=nh[i];		//	media de radios en el bin i
    aux1=m*pow(c,1.*i+0.5)*100000000000;
    aux2=m*0.5*(pow(c,1.*i)+pow(c,1.*(i+1)))*100000000000;
    factor=Masa[i]/(m*pow(c,i*1.)*(c-1.)*pow(Lbox,3.));	//funcion de masa
//    fprintf(E,"%e %e %e %e %e\n",Masa[i]*100000000000,nh[i]*factor,sqrt(nh[i])*factor,aux1,aux2);
    fprintf(E,"%e %e %e\n",Masa[i]*100000000000,nh[i]*factor,sqrt(nh[i])*factor);
  }
  fclose(E);



return 0;
}

