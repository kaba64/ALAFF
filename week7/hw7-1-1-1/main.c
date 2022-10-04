#include <stdio.h>
#include <math.h>

const int n = 80;
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

void compute_residule(double* u, double* f, double* res){
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

int main(void){
  double u[n][n],f[n][n];
  double residule;
  
  initialize_u((double *)u);
  initialize_u((double *)f);
  initialize_f((double *)f);
  residule = 1.0;
  while(residule>=1.e-8){
    for(int j=1;j<n-1;j++){
      for(int i=1;i<n-1;i++){
	u[i][j] = 0.25*(h*h*f[i][j]+u[i][j-1]+u[i][j+1]+u[i+1][j]+u[i-1][j]);
      }
    }
    compute_residule((double *)u,(double *)f,&residule);
  }
  writing_to_file((double *)u);
  return 0;
}
