#include "sor.h"
#include <math.h>
#include "parallel.h"
#include "boundary_val.h"

void sor(
        double omg,
        double dx,
        double dy,
        int    imax,
        int    jmax,
        double **P,
        double **RS,
        double *res,
        int il, 
        int ir, 
        int jt,
        int jb, 
        int rank_l, 
        int rank_r, 
        int rank_t, 
        int rank_b,
        double* bufSend,
        double* bufRecv,
        MPI_Status status,
        int myrank) 
{

  int i,j;
  double rloc;
  double coeff = omg/(2.0*(1.0/(dx*dx)+1.0/(dy*dy)));

  /* SOR iteration */
  for(i = 1; i <= imax; i++) {
    for(j = 1; j<=jmax; j++) {
      P[i][j] = (1.0-omg)*P[i][j]
              + coeff*(( P[i+1][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }

  pressure_comm(P, il, ir, jt,jb, rank_l, rank_r, rank_t, rank_b, bufSend, bufRecv, &status, myrank);

  /* compute the residual */
  rloc = 0;
  for(i = 1; i <= imax; i++) {
    for(j = 1; j <= jmax; j++) {
      rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
              ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + ( P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
    }
  }
  //rloc = rloc/(imax*jmax);
  //rloc = sqrt(rloc);

  /* set residual */
  *res = rloc;

  // Setting the boundary values for pressure based on its neighbours
  boundaryvalues_p(rank_l, rank_r, rank_t, rank_b , il , ir , jt , jb , P);
}

