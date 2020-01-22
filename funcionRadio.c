//Calcula el numero de voids a partir de las posiciones
//usando spherical-overdensity>200
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "parameters.h"
#include "RANDOM.h"



const float rho_cr=2.77526627;          //	[10¹¹ h² Ms / Mpc³]
const float pi=4.*atan(1.);             //	numero pi


const int JN=10;

int main(int argc, char **argv){

  Crandom ran2(102);

  char inFile[560];
  char outFile[560];
  char outFile2[560];

  sprintf(inFile,"%s",voids_file); printf("%s\n",inFile);
  sprintf(outFile,"dnv_%s",voids_file);
  sprintf(outFile2,"nv_%s",voids_file);




  int i, j, k, l;
  int NumVoids=0;
  int id;

  float x,y,z,radio;
  float min,max;

  float nv[JN+1][bin];
  float snv[JN+1][bin];
  float Rad[JN+1][bin];
  float sRad[JN+1][bin];
  float sigma_nv[bin];

  for(i=0;i<bin;i++){
    sigma_nv[i]=0.;
    for(j=0;j<=JN;j++){
      nv[j][i]=0.;
      snv[j][i]=0.;
      Rad[j][i]=0.;
      sRad[j][i]=0.;
    }
  }


///////////////////////////////////		masa 	maxima y minima

  printf("Leyendo el catalogo para establecer radio max y min\n");
  FILE * L;
  FILE * E;
  FILE * F;
  L=fopen(inFile,"r");
  max=0.;
  min=Lbox;
  while(fscanf(L,"%f %f %f %f\n",&x, &y, &z, &radio)!=EOF){
    //  float den1,den2;
    //  while(fscanf(L,"%f %f %f %f %f %f\n",&x, &y, &z, &radio, &den1, &den2)!=EOF){
    if(radio<min)
      min=radio;
    if(radio>max)
      max=radio;
    NumVoids++;}
  printf("Se encontraron %d voids en el catálogo\n",NumVoids);
  printf("max=%le min=%le\n",max,min);



///////////////////////		definicion de los bines radiales
  float Rmax =max*1.01;
  float r    =min*0.99;

//////////////////////  la misma escala para todas las teorias
  r=radio_min_todos;
  Rmax=radio_max_todos;

  printf("Rmax=%le r=%le\n",Rmax,r);

  float c    =pow(Rmax/r,1./bin);	//	constante de proporcionalidad log
  float aux1 =1./log(c);
  float aux2,aux3;

  printf("c=%le aux1=%le\n",c,aux1);

//////////////////////		lleno los bines, contador y radio
  rewind(L);
  for(i=0;i<NumVoids;i++){
    fscanf(L,"%f %f %f %f\n",&x, &y, &z, &radio);
    //    fscanf(L,"%f %f %f %f %f %f\n",&x, &y, &z, &radio, &den1, &den2);
    if((radio>r)&&(radio<Rmax)){
      aux2=aux1*log(radio/r);
      id=floor(aux2);		//	es el bin radial de este void
      l=floor(JN*ran2.r());		//	número entero entre 0 y 9
      if(l==JN) l=0;
      for(j=0;j<JN;j++){		//	submuestras
        if(j!=l){			//	todos las submuestras menos l
          nv[j][id]+=1.;		//	contador entero del número de voids
          Rad[j][id]+=radio;
        }
      }
    }
    nv[JN][id]+=1.;			//	total
    Rad[JN][id]+=radio; 
  }
  fclose(L);

  printf("Se encontraron %d voids en el catálogo\n",NumVoids);


	//	como nv es una densidad, y para j\in[0,JN) solo tube en cuenta una fraccion de
	//	objetos igual a float(JN-1)/JN, entonces debo multiplicarlos por JN/float(JN-1)
  float cor=float(JN)/(JN-1);


  if(acumulada==1){
  	//	errores abundancia acumulada: jacknife
    for(i=0;i<bin;i++){		//	bines
      for(j=0;j<JN;j++){
        sigma_nv[i]+=pow(nv[j][i]*cor-nv[JN][i],2);
      }
    }
    for(i=0;i<bin;i++)
      sigma_nv[i]*=float(JN-1)/JN;

    //	imprine abundancia acumulada
    F=fopen(outFile2,"w+");
    aux1=r*pow(c,1*(bin-1));
    aux2=nv[JN][bin-1];
    aux3=sigma_nv[bin-1];
    fprintf(F,"%f %e %e\n",aux1,aux2*pow(Lbox,-3.),sqrt(aux2)*pow(Lbox,-3.));
    for(i=bin-2;i>=0;i--){
      aux1=r*pow(c,i*1.);
      aux2+=nv[JN][i];
      aux3+=sigma_nv[i];
      fprintf(F,"%f %e %e\n",aux1,aux2*pow(Lbox,-3.),sqrt(aux2)*pow(Lbox,-3.));
    }
    fclose(F);
  }



  else{
    for(i=0;i<bin;i++){
      for(j=0;j<=JN;j++){
        if(nv[j][i]>0)
          Rad[j][i]/=nv[j][i];		//	media de radios en el bin i
        nv[j][i]*=Rad[j][i]/(r*pow(c,i*1.)*(c-1.));	//funcion de radio
      }
    }

    printf("Se encontraron %d voids en el catálogo\n",NumVoids);

  	//	errores abundancia: jacknife
    for(i=0;i<bin;i++){		//	bines
      for(j=0;j<JN;j++){
        sigma_nv[i]+=pow(nv[j][i]*cor-nv[JN][i],2);
      }
    }
    for(i=0;i<bin;i++)
      sigma_nv[i]*=float(JN-1)/JN;

    printf("Se encontraron %d voids en el catálogo\n",NumVoids);

    //	imprine abundancia
    E=fopen(outFile,"w+");  
    for(i=bin-1;i>=0;i--){
      fprintf(E,"%f %e %e\n",
          Rad[JN][i],nv[JN][i]*pow(Lbox,-3.),sqrt(sigma_nv[i])*pow(Lbox,-3.));
    }
/*    for(i=bin-1;i>=0;i--){
      fprintf(E,"%f %e %e ",
          Rad[JN][i],nv[JN][i]*pow(Lbox,-3.),sqrt(sigma_nv[i])*pow(Lbox,-3.));
      for(j=0;j<JN;j++)
        fprintf(E,"%e ",nv[j][i]*pow(Lbox,-3.)*cor);
      fprintf(E,"\n");
    }*/
    fclose(E);
  }



return 0;
}

