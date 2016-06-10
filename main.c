#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include <mpi.h>

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed using the matrix() command
 * - create the initial setup init_uvp(), init_flag(), output_uvp()
 * - perform the main loop
 * - trailer: destroy memory allocated and do some statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two dimensional degree of freedom layout, so they are not stored in
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop the following big steps are done (for some of the 
 * operations a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */

int main(int argn, char** args)
{
   double Re;               
   double UI;                
   double VI;               
   double PI;           
   double GX;                
   double GY;                
   double t_end;            
   double xlength;          
   double ylength;           
   double dt;                
   double dx;             
   double dy;               
   int  imax;                
   int  jmax;               
   double alpha;            
   double omg;               
   double tau;              
   int  itermax;             
   double eps;              
   double dt_value;
   int iproc;
   int jproc;
   double** U;
   double** V;
   double** P;          
   double** RS;
   double** F;
   double** G;
   int myrank;
   int num_proc;
   int *omg_i;
   int *omg_j;
   int *il;
   int *ir;
   int *jb;
   int *jt;
   int *rank_l;
   int *rank_r;
   int *rank_b;
   int *rank_t;
   int basic_info[6];   // Is of size 8 because 1st 2 will save the position of the block wrt to the global system, next 2 the size of that block and last 4 will save the left,right,bottom and top neighbour 


//setting the parameters
read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,&itermax, &eps, &dt_value,&iproc,
&jproc);


// Initializing MPI
	MPI_Init(&argn, &args);
	MPI_Comm_size( MPI_COMM_WORLD, &num_proc ); 
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
  {
    init_parallel (iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t, 
                   &omg_i , &omg_j, num_proc );
  }
	



	}


































// Creating the arrays U,V and P
	U = matrix ( 0 , imax+1 , 0 , jmax+1 );
        V = matrix ( 0 , imax+1 , 0 , jmax+1 );
        P = matrix ( 0 , imax+1 , 0 , jmax+1 );
        

// Creating arrays for right side of pressure poissons equation (RS) and F and G
	RS = matrix ( 0,imax+1,0,jmax+1);
        F = matrix (0,imax+1,0,jmax+1);
        G = matrix (0,imax+1,0,jmax+1);

// Initializing the arrays U,V,P,RS,F and G
        init_uvp( UI, VI,PI,imax, jmax,U,V,P);
        init_matrix(RS,0,imax+1,0,jmax+1,0);
        init_matrix(F,0,imax+1,0,jmax+1,0);
        init_matrix(G,0,imax+1,0,jmax+1,0);


        double t=0;   // initialize the time
        int n = 0;    // number of time steps

while (t<t_end)
  {
      calculate_dt(Re,tau,&dt,dx,dy,imax,jmax,U,V);
      boundaryvalues(imax, jmax, U, V,P, G, F);
      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, imax, jmax, U, V, F, G);
      calculate_rs(dt,dx,dy, imax,jmax, F, G, RS);
      int it = 0;
      double res = 1000;

      while(it<itermax && res > eps) 
          {
            sor(omg, dx,dy,imax,jmax, P, RS, &res);
            it++; 
          }
	if (it>=itermax-1 && res > eps){
	printf("Not converged in %d iterations  Residual = %f \n",it, res);
	}
	else
	printf("Converged in %d iterations  Residual = %f \n",it, res);

      calculate_uv(dt,dx, dy,imax,jmax,U,V,F,G,P);
      t = t+dt;
      n = n+1;
  }


write_vtkFile("szProblem.vtk", n, xlength, ylength, imax, jmax,dx, dy, U, V, P);
free_matrix(U,0,imax+1,0,jmax+1);
free_matrix(V,0,imax+1,0,jmax+1);
free_matrix(P,0,imax+1,0,jmax+1);
free_matrix(RS,0,imax+1,0,jmax+1);
free_matrix(F,0,imax+1,0,jmax+1);
free_matrix(G,0,imax+1,0,jmax+1);


  return 0;
}
