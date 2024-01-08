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
  // Richardson method with alpha given
  // AB: matrix of size lab*la
  // RHS: right hand side of size la
  // X: solution of size la
  // alpha_rich: value of alpha
  // lab: number of bands of AB
  // la: size of the matrix
  // ku: number of upper bands of AB
  // kl: number of lower bands of AB
  // tol: tolerance for the stopping criterion
  // maxit: maximum number of iterations
  // resvec: vector of size maxit containing the residual at each iteration
  // nbite: number of iterations
  // PI: 3.14159265358979323846
  //use richardson_alpha_opt
   
  //initialization
  int i;
  for(i=0;i<*la;i++){
    X[i]=0;
  }

  //iteration
  double *R;
  R=(double *) malloc(sizeof(double)*(*la));
  double *AX;
  AX=(double *) malloc(sizeof(double)*(*la));
  double *res;
  res=(double *) malloc(sizeof(double)*(*la));

  double alpha=*alpha_rich;
  double norm_res=1;
  int it=0;
  while(norm_res>*tol && it<*maxit){
    //compute residual
    //dgbmv_("N",la,la,kl,ku,1.0,AB,lab,X,1,0.0,AX,1);
    
    //dgbmv_("N",la,la,kl,ku,1.0,AB,lab,X,1,0.0,AX,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,1.0,AB,*lab,X,1,0.0,AX,1);
   // mat_vec_prod(AB,X,AX,lab,la,ku,kl);
    for(i=0;i<*la;i++){
      R[i]=RHS[i]-AX[i];
    }
    //compute new solution
    for(i=0;i<*la;i++){
      X[i]=X[i]+alpha*R[i];
    }
    //compute new residual
    // use blas function
    //dgbmv_("N",la,la,kl,ku,1.0,AB,lab,X,1,0.0,AX,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,1.0,AB,*lab,X,1,0.0,AX,1);
    //mat_vec_prod(AB,X,AX,lab,la,ku,kl);

    for(i=0;i<*la;i++){
      res[i]=RHS[i]-AX[i];
    }
    //compute norm of the residual without auxiliary function
    norm_res=0;
    for(i=0;i<*la;i++){
      norm_res=norm_res+res[i]*res[i];
    }
    //store the norm of the residual
    resvec[it]=norm_res;
    //update the number of iterations
    it++;
  }
  *nbite=it;
  free(R);
  free(AX);
  free(res);


 

  

}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv) {
  
  int i, j, k,km;
 //int labm=(*ku)+(*kl)+(*kv)+1;
  for (j=0;j<(*la);j++){
    k = j*(*lab);
   // km = j*labm;
    if (*kv>=0){
      for (i=0;i<*kv;i++){
	      MB[k+i]=0.0;
      }
    }
    MB[k+*kv]=AB[k+*kv];
  }
  MB[0]=0.0;
 MB[(*lab)*(*la)-1]=0.0;
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
 
 // X(k+1)=X(k)+MB^-1*(RHS-AB*X(k))
  // MB^-1=diag(MB)^-1

  //initialization
  int i;
  for(i=0;i<*la;i++){
    X[i]=0;
  }

  //iteration
  double *R;
  R=(double *) malloc(sizeof(double)*(*la));
  double *AX;
  AX=(double *) malloc(sizeof(double)*(*la));
  double *res;
  res=(double *) malloc(sizeof(double)*(*la));

  double norm_res=1;
  int it=0;

  //R[0]=RHS[0]-AB[0]*X[0]

  while(norm_res>*tol && it<*maxit){
    //AX=AB*X
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,1.0,AB,*lab,X,1,0.0,AX,1);
    //R=RHS-AX
    for(i=0;i<*la;i++){
      R[i]=RHS[i]-AX[i];
    }
    //X=X+MB^-1*R
    for(i=0;i<*la;i++){
      X[i]=X[i]+R[i]/MB[i*(*lab)+*kl];
    }
    //compute new residual
    // use blas function
    //dgbmv_("N",la,la,kl,ku,1.0,AB,lab,X,1,0.0,AX,1);
    cblas_dgbmv(CblasColMajor,CblasNoTrans,*la,*la,*kl,*ku,1.0,AB,*lab,X,1,0.0,AX,1);

    for(i=0;i<*la;i++){
      res[i]=RHS[i]-AX[i];
    }
    //compute norm of the residual without auxiliary function
    norm_res=0;
    for(i=0;i<*la;i++){
      norm_res=norm_res+res[i]*res[i];
    }
    //store the norm of the residual
    resvec[it]=norm_res;
    //update the number of iterations
    it++;

    
  }

   *nbite=it;
  free(R);
  free(AX);
  free(res);
  



  

}

