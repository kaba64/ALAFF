void Jacobi(double* u, double * uold, double* f, double h, int n){

  for(int j=1;j<n-1;j++){
    for(int i=1;i<n-1;i++){
      u[n*j+i] = 0.25*(h*h*f[j*n+i]+uold[j*n+i-n]+uold[j*n+i+n]+uold[j*n+i+1]+uold[j*n+i-1]);
      uold[n*j+i]=u[n*j+i];
    }
  }
}
