/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=-1.0;
    AB[kk+ *kv+1]=2.0;
    AB[kk+ *kv+2]=-1.0;
  }
  AB[0]=0.0;
  if (*kv == 1) {AB[1]=0;}
  
  AB[(*lab)*(*la)-1]=0.0;
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  int ii, jj, kk;
  for (jj=0;jj<(*la);jj++){
    kk = jj*(*lab);
    if (*kv>=0){
      for (ii=0;ii< *kv;ii++){
	AB[kk+ii]=0.0;
      }
    }
    AB[kk+ *kv]=0.0;
    AB[kk+ *kv+1]=1.0;
    AB[kk+ *kv+2]=0.0;
  }
  AB[1]=0.0;
  AB[(*lab)*(*la)-1]=0.0;
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  int jj;
  RHS[0]= *BC0;
  RHS[(*la)-1]= *BC1;
  for (jj=1;jj<(*la)-1;jj++){
    RHS[jj]=0.0;
  }
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  int jj;
  double h, DELTA_T;
  DELTA_T=(*BC1)-(*BC0);
  for (jj=0;jj<(*la);jj++){
    EX_SOL[jj] = (*BC0) + X[jj]*DELTA_T;
  }
}  

void set_grid_points_1D(double* x, int* la){
  int jj;
  double h;
  h=1.0/(1.0*((*la)+1));
  for (jj=0;jj<(*la);jj++){
    x[jj]=(jj+1)*h;
  }
}

// double relative_forward_error(double* x, double* y, int* la){
//   // Compute the relative forward error
//   // x: computed solution
//   // y: exact solution
//   // la: size of the solution
//   int jj;
//   double err, normy;
//   err=0.0;
//   normy=0.0;
//   for (jj=0;jj<(*la);jj++){
//     err=err+pow(x[jj]-y[jj],2);
//     normy=normy+pow(y[jj],2);
//   }
//   err=sqrt(err);
//   normy=sqrt(normy);
//   return err/normy;

// }


double relative_forward_error(double* x, double* y, int* la) {
    // Compute the norm of the difference between x and y
    double* diff = (double*)malloc((*la) * sizeof(double));
    cblas_dcopy(*la, y, 1, diff, 1);
    cblas_daxpy(*la, -1.0, x, 1, diff, 1);

    double norm_diff = cblas_dnrm2(*la, diff, 1);

    // Compute the norm of y
    double norm_y = cblas_dnrm2(*la, y, 1);

    // Compute and return the relative forward error
    double rel_forward_error = norm_diff / norm_y;

    // Free allocated memory
    free(diff);

    return rel_forward_error;
}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

  
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info) {
    // Check for valid input
    if (*la <= 0 || *n <= 0 || *kl < 0 || *ku < 0 || *lab < 2 * *kl + *ku + 1 || !AB || !ipiv || !info) {
        *info = -1;  // Invalid input
        return *info;
    }

    // LAPACK function call for LU factorization of a tridiagonal band matrix
    dgbtrf_(n, n, kl, ku, AB, lab, ipiv,info);

    return *info;
}

