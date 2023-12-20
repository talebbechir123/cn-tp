/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void eig_poisson1D(double* eigval, int *la){
  // Compute the eigenvalues of the 1D Poisson matrix
  // eigval: array of size la
  // la: size of the matrix
  // PI: 3.14159265358979323846
  double PI=3.14159265358979323846;
  int i;
  for(i=0;i<*la;i++){
    eigval[i]=2*(1-cos((i+1)*PI/(*la+1)));
  }
}

double eigmax_poisson1D(int *la){
  // Compute the maximum eigenvalue of the 1D Poisson matrix
  // la: size of the matrix
  // PI: 3.14159265358979323846
  //use eig_poisson1D
  double* eigval;
  eigval=(double *) malloc(sizeof(double)*(*la));
  eig_poisson1D(eigval,la);
  //look for the maximum eigenvalue
  double eigmax=eigval[0];
  int i;
  for(i=1;i<*la;i++){
    if(eigval[i]>eigmax){
      eigmax=eigval[i];
    }
  }
  free(eigval);
  return eigmax;
}

double eigmin_poisson1D(int *la){
  // Compute the minimum eigenvalue of the 1D Poisson matrix
  // la: size of the matrix
  // PI: 3.14159265358979323846
  //use eig_poisson1D
  double* eigval;
  eigval=(double *) malloc(sizeof(double)*(*la));
  eig_poisson1D(eigval,la);
  //look for the minimum eigenvalue
  double eigmin=eigval[0];
  int i;
  for(i=1;i<*la;i++){
    if(eigval[i]<eigmin){
      eigmin=eigval[i];
    }
  }
  free(eigval);
  return eigmin;
}

double richardson_alpha_opt(int *la){
  // Optimal value of alpha for Richardson method
  // la: size of the matrix
  // PI: 3.14159265358979323846
  //use eigmax_poisson1D and eigmin_poisson1D
  double eigmax=eigmax_poisson1D(la);
  double eigmin=eigmin_poisson1D(la);
  double alpha_opt=2/(eigmax+eigmin);
  return alpha_opt;
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
}

