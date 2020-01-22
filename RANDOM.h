#include<iostream>
#include<cmath>
using namespace std;

//Constantes del generador aleatorio
#define NTAB 32
const int IM1=2147483563;//numeros magicos
const int  IM2=2147483399;//""
const double AM=(1.0/IM1);
const int IMM1=(IM1-1);
const int IA1=40014;
const int IA2=40692;
const int IQ1=53668;
const int IQ2=52774;
const int IR1=12211;
const int IR2=3791;
const int  NDIV=(1+IMM1/NTAB);
const double EPS=1.2e-7;
const double RNMX=(1.0-EPS);


class Crandom{
  long int idum,idum2,iy,iv[NTAB];
public:
  Crandom(int seed);//arrancar(aca hay una funcion que sellama igual a la clases:constructor)
  ~Crandom(){};//destructor
  float r(void);//generar un número
  double exponencial(float tau);
  double gauss(float mu,float sigma);
};


Crandom::Crandom(int seed){
  int k,j;
  iy=0;
  if (seed<1) idum=1; else idum=seed;//ambos generadores aleatorios estan empezando con la misma semilla
  idum2=idum;
  for (j=NTAB+7;j>=0;j--){
    //idum=(idum*IA1) % IM1
    k=idum/IQ1;idum=IA1*(idum-k*IQ1)-k*IR1;if (idum<0) idum+=IM1;
    if(j<NTAB) iv[j]=idum;//guarda en la columna de datos
  }
  iy=iv[0];
}
float Crandom:: r(void){
  int k,j; float rout;
  //idum=(idum*ia1)%m1
  k=idum/IQ1;idum=IA1*(idum-k*IQ1)-k*IR1; if(idum<0) idum+=IM1;           
  //idum2=(idum2*ia2)%m2
  k=idum2/IQ2;idum2=IA2*(idum2-k*IQ2)-k*IR2;if(idum2<0) idum2+=IM2;
  j=iy/NDIV;iy=iv[j]-idum2;iv[j]=idum;
  if(iy<1) iy+=IMM1;
  rout=AM*iy; if (rout>RNMX) return RNMX; else return rout;
}
double Crandom::exponencial(float tau){
  return -tau*log(r());
}
double Crandom:: gauss(float mu,float sigma){
  return sigma*sqrt(-2*log(r()))*cos(2*M_PI*r())+mu;
}


/*int main(){

  
  return 0;
  }*/
