#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"

int main(int argn, char** args){
  // Initializing MPI
  MPI_Init(&argn, &args);

  // Defining sufficient variables.
    double Re;               
    double UI;                
    double VI;               
    double PI;           
    double GX;                
    double GY;                
    double t_end;            
    double xlength;          
    double ylength;                        
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
    double dt; 
    double** U;
    double** V;
    double** P;          
    double** RS;
    double** F;
    double** G;
    int myrank;
    int num_proc;
    int omg_i ;
    int omg_j ;
    int il ;
    int ir ;
    int jb;
    int jt ;
    int rank_l ;
    int rank_r ;
    int rank_b ;
    int rank_t ;
    MPI_Status status;

    MPI_Comm_size( MPI_COMM_WORLD, &num_proc ); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);          

  // Reading necessary parameters from the file problem.dat
    read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,&itermax, 
    &eps, &dt_value,&iproc,&jproc);

  // Finding out the left,right,top,bottom neighbours and the size by providing the rank, iproc, jproc, and num_proc to the function.
    init_parallel (iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t, &omg_i , &omg_j, num_proc );

  //printf("\n %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n" , myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i , omg_j);

  // Creating and then initializing all the necessary mmatrices
    U = matrix ( 0 , ir - il  + 3 , 0 , jt-jb+2 );
    V = matrix ( 0 , ir - il  + 2 , 0 , jt-jb+3 );
    P = matrix ( 0 , ir - il  + 2 , 0 , jt-jb+2 );
    RS = matrix ( 0 , ir - il  , 0 , jt-jb);
    F = matrix (0 , ir - il  + 3 , 0 , jt-jb+2);
    G = matrix (0 , ir - il  + 2 , 0 , jt-jb+3);

    init_matrix(U , 0 , ir - il  + 3 , 0 , jt-jb+2 , UI);
    init_matrix(V , 0 , ir - il  + 2 , 0 , jt-jb+3 , VI);
    init_matrix(P , 0 , ir - il  + 2 , 0 , jt-jb+2 , PI);
    init_matrix(RS, 0 , ir - il  , 0 , jt-jb  , 0 );
    init_matrix(F , 0 , ir - il  + 3 , 0 , jt-jb+2 , 0 );
    init_matrix(G , 0 , ir - il  + 2 , 0 , jt-jb+3 , 0 );

// Allocating sufficient memory to the send and receive buffer.
    int biggest;
    if (ir-il>jt-jb)
    {
    biggest = ir-il+4;
    }   
    else
    {
    biggest = jt-jb+4;
    }

    double* bufSend = malloc(biggest*sizeof(double));
    double* bufRecv =  malloc(biggest*sizeof(double));

    // initialize the time
    double t=0;

    // number of time steps
    int n = 0; 
  
    int vtk_file_num = 1;

    while (t<=t_end)
    {
    // Setting boundary values for U,V, F and G  
      boundary_values(rank_l, rank_r , rank_t , rank_b , il , ir , jt , jb , U , V ,F, G);
    
    // Calculating the F and G matrices  
      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, ir - il + 2 , jt - jb + 2  , U, V, F, G, rank_l, rank_r, rank_t, rank_b );

    // Calculating the right hand side of the pressure poissions equation
      calculate_rs(dt,dx,dy, ir - il  , jt - jb , F , G , RS);

    // Iteration number
      int it = 0;

    // Residual
      double res = 1000;

            
      while(it<itermax && res > eps) 
      {
      // Setting boundary values for pressure
        boundaryvalues_p(rank_l, rank_r, rank_t, rank_b , il, ir, jt,jb , P);

      // Doing successive over relaxation
        sor(omg, dx , dy , ir - il , jt - jb   , P, RS, &res, il, ir, jt,jb, rank_l, rank_r, rank_t, rank_b, bufSend, bufRecv, status, myrank);

      //Taking the residual from each processor and then summing it up to find the residual of the whole problem and then broadcasting it.
        double res_new ;
        MPI_Allreduce(&res , &res_new , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ); 
	      res_new = res_new/(imax*jmax);
	      res = sqrt(res_new);

        it++;
      }

      // Calculaing the matrices U and V
        calculate_uv( dt , dx , dy , ir - il , jt - jb   , U , V , F , G , P);

      // Communicating the values of U and V  
        uv_comm(U, V, il, ir, jt,jb, rank_l, rank_r, rank_t, rank_b, bufSend, bufRecv, &status, myrank); 

      // Calculating the time step for each processor
        calculate_dt(Re , tau , &dt , dx , dy , ir - il +1 , jt - jb + 1  , U , V );

      // collecting dt from each processor and then finding out the minimum of it and finally broadcasting it.
        double dt_new ;

        MPI_Allreduce(&dt , &dt_new , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

        dt = dt_new;

        t = t+dt;

        n = n+1;

      // Printing the vtk files after the required interval.
        if ((int)t/(int)dt_value == vtk_file_num )
        {
          output_uvp(U, V, P, il, ir, jb, jt, omg_i, omg_j,"szProblem1.vtk",n,myrank);
          vtk_file_num++;
        }

        printf("Processor Num = %d Converged in iteration = %d dt value = %f residual = %f \n ", myrank , it , dt , res);

   }
    
    //Freeing all the momory allocated on the processor.
    free_matrix ( U, 0 , ir - il  + 3 , 0 , jt-jb+2 );
    free_matrix ( V, 0 , ir - il  + 2 , 0 , jt-jb+3);
    free_matrix ( P, 0 , ir - il  + 2 , 0 , jt-jb+2);
    free_matrix ( RS, 0 , ir - il  , 0 , jt-jb);
    free_matrix ( F, 0 , ir - il  + 3 , 0 , jt-jb+2);
    free_matrix ( G, 0 , ir - il  + 2 , 0 , jt-jb+3);
    free(bufSend);
    free(bufRecv);

    Programm_Sync("Synchronizing all the processors");
    Programm_Stop("Stoping Parallel Run");
    return 0;
}               
