#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"

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
   int *omg_i = NULL;
   int *omg_j = NULL;
   int *il = NULL;
   int *ir = NULL;
   int *jb = NULL;
   int *jt = NULL;
   int *rank_l = NULL;
   int *rank_r = NULL;
   int *rank_b = NULL;
   int *rank_t = NULL;
   MPI_Status status;
  
//setting the parameters
read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,&itermax, &eps, &dt_value,&iproc,
&jproc);


// Initializing MPI
      	MPI_Init(&argn, &args);
      	MPI_Comm_size( MPI_COMM_WORLD, &num_proc ); 
      	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


if (myrank == 0){

        omg_i = malloc(iproc*sizeof(int));
        omg_j = malloc(jproc*sizeof(int));
        il = malloc(iproc*jproc*sizeof(int));
        ir = malloc(iproc*jproc*sizeof(int));
        jt = malloc(iproc*jproc*sizeof(int));
        jb = malloc(iproc*jproc*sizeof(int));
        rank_l = malloc(iproc*jproc*sizeof(int));
        rank_r = malloc(iproc*jproc*sizeof(int));
        rank_t = malloc(iproc*jproc*sizeof(int));
        rank_b = malloc(iproc*jproc*sizeof(int));
        
        init_parallel (iproc, jproc, imax, jmax, &myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i , omg_j, num_proc );
        printf("main after init p\n");
        printf("Rank in main %d %d %d %d \n", *rank_l, *rank_r, *rank_t, *rank_b);  //1st element
        printf("omg_i in main %d %d %d %d \n", omg_i[0], omg_i[1], omg_j[0], omg_j[1]);

        int* array_pos = malloc(2*sizeof(int));
        int* array_size = malloc(4*sizeof(int));
        int* array_neighbours = malloc(4*sizeof(int));

        
        array_pos[0] = 0;
        array_pos[1] = 0;
 
        array_size[0] = il[0];
        array_size[1] = ir[0];
        array_size[2] = jt[0];
        array_size[3] = jb[0];
                

        array_neighbours[0] = rank_l[0];
        array_neighbours[1] = rank_r[0];
        array_neighbours[2] = rank_t[0];
        array_neighbours[3] = rank_b[0];  


/*// Creating the arrays U,V and P
        U = matrix ( il[0] - 2 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 );
        V = matrix ( il[0] - 1 , ir[0] + 1 , jb[0] - 2 , jt[0]+1 );
        P = matrix ( il[0] - 1 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 );
          

// Creating arrays for right side of pressure poissons equation (RS) and F and G
        RS = matrix ( il[0],ir[0],jb[0],jt[0]);
        F = matrix ( il[0] - 2 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 );
        G = matrix ( il[0] - 1 , ir[0] + 1 , jb[0] - 2 , jt[0]+1 );

// Initializing the arrays U,V,P,RS,F and G
        init_matrix(U , il[0] - 2 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 , UI);
        init_matrix(V , il[0] - 1 , ir[0] + 1 , jb[0] - 2 , jt[0]+1 , VI);
        init_matrix(P , il[0] - 1 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 , PI);
        init_matrix(RS , il[0] , ir[0] , jb[0] , jt[0] , 0);
        init_matrix(F , il[0] - 2 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 , 0);
        init_matrix(G , il[0] - 1 , ir[0] + 1 , jb[0] - 2 , jt[0]+1 , 0);*/

        U = matrix ( array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        V = matrix ( array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 );
        P = matrix ( array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        

// Creating arrays for right side of pressure poissons equation (RS) and F and G
        RS = matrix ( array_size[0],array_size[1],array_size[3],array_size[2]);
        F = matrix ( array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        G = matrix ( array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 );

// Initializing the arrays U,V,P,RS,F and G
        init_matrix(U , array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 , UI);
        init_matrix(V , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 , VI);
        init_matrix(P , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 , PI);
        init_matrix(RS , array_size[0] , array_size[1] , array_size[3] , array_size[2] , 0);
        init_matrix(F , array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 , 0);
        init_matrix(G , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 , 0);
        Program_Message("we are here \n \n \n \n");
// initialize the time
        double t=0;
// number of time steps
        int n = 0;   


while (t<t_end)
  {
      calculate_dt(Re , tau , &dt , dx , dy , ir[0] - il[0] + 1, jt[0] - jb[0] +1 , U , V );
      boundary_values(rank_l[0], rank_r[0], rank_t[0], rank_b[0] ,il[0], ir[0], jt[0],jb[0], U , V , P);
      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, ir[0] - il[0] + 1, jt[0] - jb[0] +1 , U, V, F, G);
      calculate_rs(dt,dx,dy, ir[0] - il[0] + 1 , jt[0] - jb[0] +1 , F , G , RS);
      int it = 0;
      double res = 1000;

      while(it<itermax && res > eps) 
          {
            sor(omg, dx,dy,ir[0] - il[0] + 1,jt[0] - jb[0] +1 , P, RS, &res);
            int compRes;
            for (int i = 1; i < num_proc; ++i)
            {
              MPI_Recv(&compRes, 1, MPI_DOUBLE, i , 20 , MPI_COMM_WORLD, &status);
                res += compRes;             
            }
            //this command sends the maximum res to all processes. To receive value same command has to be given in other processes
              MPI_Bcast(&res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
              it++;   //we do not need to check iteration number all the processors will have same iteration number
          }
      if (it>=itermax-1 && res > eps)
      printf("Not converged in %d iterations  Residual = %f \n",it, res);
      else
      printf("Converged in %d iterations  Residual = %f \n",it, res);

      calculate_uv( dt , dx , dy , ir[0] - il[0] + 1,jt[0] - jb[0] +1 , U , V , F , G , P);
      t = t+dt;
      n = n+1;
  }
        write_vtkFile("szProblem.vtk", n, xlength, ylength, iproc, jproc,dx, dy, U, V, P,myrank);
        free(omg_j);
        free(omg_i);
        free(il);
        free(ir);
        free(jt);
        free(jb);
        free(rank_l);
        free(rank_r);
        free(rank_t);
        free(rank_b);
        free_matrix(U , il[0] - 2 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 );
        free_matrix(V , il[0] - 1 , ir[0] + 1 , jb[0] - 2 , jt[0]+1 );
        free_matrix(P , il[0] - 1 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 );
        free_matrix(RS , il[0] , ir[0] , jb[0] , jt[0] );
        free_matrix(F , il[0] - 2 , ir[0] + 1 , jb[0] - 1 , jt[0]+1 );
        free_matrix(G , il[0] - 1 , ir[0] + 1 , jb[0] - 2 , jt[0]+1 );



}
        




	
else{

// Creating the arrays U,V and P
        //int* array_pos = malloc(2*sizeof(int));
        int* array_size = malloc(4*sizeof(int));
        int* array_neighbours = malloc(4*sizeof(int));

        
        MPI_Recv(array_neighbours,4,MPI_INT,0,2,MPI_COMM_WORLD,&status);
        printf("processes %d \n", myrank);
        printf("Neighbours %d %d %d %d\n", *(array_neighbours + 0),*(array_neighbours + 1),*(array_neighbours + 2),*(array_neighbours + 3));
        MPI_Recv(array_size,4,MPI_INT,0,3,MPI_COMM_WORLD,&status);
        printf("Size %d %d %d %d\n", *(array_size + 0),*(array_size + 1),*(array_size + 2),*(array_size + 3));
        //MPI_Recv(array_pos,2,MPI_INT,0,1,MPI_COMM_WORLD,&status);
        //printf("Position %d %d \n", *(array_pos + 0),*(array_pos + 1));

        


        U = matrix ( array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        V = matrix ( array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 );
        P = matrix ( array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        

// Creating arrays for right side of pressure poissons equation (RS) and F and G
        RS = matrix ( array_size[0],array_size[1],array_size[3],array_size[2]);
        F = matrix ( array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        G = matrix ( array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 );

// Initializing the arrays U,V,P,RS,F and G
        init_matrix(U , array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 , UI);
        init_matrix(V , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 , VI);
        init_matrix(P , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 , PI);
        init_matrix(RS , array_size[0] , array_size[1] , array_size[3] , array_size[2] , 0);
        init_matrix(F , array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 , 0);
        init_matrix(G , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 , 0);
        

// initialize the time
        double t=0;
// number of time steps
        int n = 0;

while (t<t_end)
  {
      calculate_dt(Re , tau , &dt , dx , dy , array_size[1] - array_size[0] + 1, array_size[2] - array_size[3] +1 , U , V );
      boundary_values(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3], U , V , P);
      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, array_size[1] - array_size[0] + 1, array_size[2] - array_size[3] +1 , U, V, F, G);
      calculate_rs(dt,dx,dy, array_size[1] - array_size[0] + 1 , array_size[2] - array_size[3] +1 , F , G , RS);
      int it = 0;
      double res = 1000;


      while(it<itermax && res > eps) 
          {
            sor(omg, dx,dy,array_size[1] - array_size[0] + 1,array_size[2] - array_size[3] +1 , P, RS, &res);
            MPI_Send(&res, 1 , MPI_DOUBLE, 0, 20 , MPI_COMM_WORLD);
            MPI_Bcast(&res, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            it++; 
          }
	/*if (it>=itermax-1 && res > eps){
	printf("Not converged in %d iterations  Residual = %f \n",it, res);
	}
	else
	printf("Converged in %d iterations  Residual = %f \n",it, res);*/

      calculate_uv( dt , dx , dy , array_size[1] - array_size[0] + 1,array_size[2] - array_size[3] +1 , U , V , F , G , P);
      t = t+dt;
      n = n+1;
  }


        write_vtkFile("szProblem.vtk", n, xlength, ylength, iproc, jproc,dx, dy, U, V, P,myrank);
        free_matrix(U , array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        free_matrix(V , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 );
        free_matrix(P , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        free_matrix(RS , array_size[0] , array_size[1] , array_size[3] , array_size[2] );
        free_matrix(F , array_size[0] - 2 , array_size[1] + 1 , array_size[3] - 1 , array_size[2]+1 );
        free_matrix(G , array_size[0] - 1 , array_size[1] + 1 , array_size[3] - 2 , array_size[2]+1 );

Programm_Sync("Synchronizing all the processors");
Programm_Stop("Stoping Parallel Run");
return 0;
}
}
