#include "helper.h"
#include "visual.h"
#include "init.h"
#include "boundary_val.h"
#include "uvp.h"
#include "sor.h"
#include <stdio.h>
#include <mpi.h>
#include "parallel.h"


int main(int argn, char** args)
{
  // Initializing MPI
  MPI_Init(&argn, &args);

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

  read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,&itermax, 
  &eps, &dt_value,&iproc,&jproc);

// Running the process for master processor

    init_parallel (iproc, jproc, imax, jmax, &myrank, &il, &ir, &jb, &jt, &rank_l, &rank_r, &rank_b, &rank_t, &omg_i , &omg_j, num_proc );

    printf("\n %d  %d  %d  %d  %d  %d  %d  %d  %d  %d  %d\n" , myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i , omg_j);

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
  
    while (t<t_end)
    {
      boundary_values(rank_l, rank_r , rank_t , rank_b , il , ir , jt , jb , U , V ,F, G);
      
      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, ir - il +2 , jt - jb + 2  , U, V, F, G, rank_l, rank_r, rank_t, rank_b );

      calculate_rs(dt,dx,dy, ir - il  , jt - jb , F , G , RS);

      int it = 0;

      double res = 1000;

            
      while(it<itermax && res > eps) 
      {
        boundaryvalues_p(rank_l, rank_r, rank_t, rank_b , il, ir, jt,jb , P);

        sor(omg, dx , dy , ir - il , jt - jb   , P, RS, &res, il, ir, jt,jb, rank_l, rank_r, rank_t, rank_b, bufSend, bufRecv, status, myrank);
        
        printf("old residual = %f\n", res );

        double res_new ;

        MPI_Allreduce(&res , &res_new , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ); 

        res = res_new;

        printf("new residual = %f\n", res );

        it++;
      }

        calculate_uv( dt , dx , dy , ir - il , jt - jb   , U , V , F , G , P);

        uv_comm(U, V, il, ir, jt,jb, rank_l, rank_r, rank_t, rank_b, bufSend, bufRecv, &status, myrank); 

        calculate_dt(Re , tau , &dt , dx , dy , ir - il +1 , jt - jb + 1  , U , V );

        double dt_new ;

        MPI_Allreduce(&dt , &dt_new , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

        dt = dt_new;

        t = t+dt;

        n = n+1;
   }



    write_vtkFile("szProblem.vtk", n , xlength/2.0 , ylength/2.0 , ir - il , jt - jb,dx, dy, U, V, P,il,ir,jt,jb,myrank);
    
    //Freeing all the momory allocated on the processor.

    free_matrix ( U, 0 , ir - il  + 3 , 0 , jt-jb+2 );
    free_matrix ( V, 0 , ir - il  + 2 , 0 , jt-jb+3);
    free_matrix ( P, 0 , ir - il  + 2 , 0 , jt-jb+2);
    free_matrix ( RS, 0 , ir - il  , 0 , jt-jb);
    free_matrix ( F, 0 , ir - il  + 3 , 0 , jt-jb+2);
    free_matrix ( G, 0 , ir - il  + 2 , 0 , jt-jb+3);


    Programm_Sync("Synchronizing all the processors");
    Programm_Stop("Stoping Parallel Run");
    return 0;
}               