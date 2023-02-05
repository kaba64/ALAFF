#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Gauss_Seidel.h"
#include "Jacobi.h"

const int n = 200;
const double h = 1.0/(n-1);
const double alpha = 1.0;
const double beta = 1.0;
const double pi = acos(-1.);

void initialize_u(double * f){
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      f[j*n+i]=0.;
    }
  }
}

void initialize_f(double* f){
  for(int j=1;j<n-1;j++){
    for(int i=1;i<n-1;i++){
      f[j*n+i]= (alpha*alpha+beta*beta)*(pi*pi)*sin(alpha*pi*h*i)*cos(beta*pi*h*j);
    }
  }
}

void writing_to_file(double* f){

  FILE *ofile;
  char filename[] = "u.txt";
  ofile = fopen(filename, "w");
  for(int j=0;j<n;j++){
    for(int i=0;i<n;i++){
      fprintf(ofile, "%d\t%d\t%lf\t%lf\t%lf\n",i,j,i*h,j*h,f[j*n+i]);
    }
  }
  fclose(ofile);
}

void writing_residual(char name, double res, int num){

  FILE *ofile;
  char filename[] = "residual.txt";
  char name1[3] = {name ,'_','\0'};
  ofile = fopen(strcat(name1,filename), "a");
  fprintf(ofile, "%d\t%.10f\n",num,res);
  fclose(ofile);
}

void compute_residual(double* u, double* f, double* res){
  double temp;
  temp = fabs(u[n+1]-0.25*(h*h*f[n+1]+u[n+1-n]+u[n+1+n]+u[n+1+1]+u[n+1-1]));
  *res=temp;
  for(int j=1;j<n-1;j++){
    for(int i=1;i<n-1;i++){
      temp=fabs(u[j*n+i]-0.25*(h*h*f[j*n+i]+u[j*n+i-n]+u[j*n+i+n]+u[j*n+i+1]+u[j*n+i-1]));
      if(temp>=*res){
	*res=temp;
      }
    }
  }
}

int main(int argc, char *argv[]){
  double u[n][n], un[n][n],f[n][n];
  double residual;
  char name;
  int count=0;
  initialize_u((double *)u);
  initialize_u((double *)un);
  initialize_u((double *)f);
  initialize_f((double *)f);
  residual = 1.0;
  if(argc!=2){
    printf("Please enter the name of the iterative method.\n");
  }else{
    name = *argv[1];
  }
  if(name=='G'){
    while(residual>=1.e-6){
      Gauss_Seidel((double*)u,(double *)un,(double*)f,h,n);
      count++;
      compute_residual((double *)u,(double *)f,&residual);
      writing_residual('G',residual,count);
    }
  }else if(name=='J'){
    while(residual>=1.e-6){
      Jacobi((double*)u,(double *)un,(double*)f,h,n);
      count++;
      compute_residual((double *)u,(double *)f,&residual);
      writing_residual('J',residual,count);
    }
  }
  writing_to_file((double *)u);
  return 0;
}
