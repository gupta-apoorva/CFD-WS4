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
  double* bufSend = NULL;
  double* bufRecv = NULL;
  MPI_Status status;

  MPI_Comm_size( MPI_COMM_WORLD, &num_proc ); 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);          

  read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,&itermax, 
  &eps, &dt_value,&iproc,&jproc);

// Running the process for master processor
  if (myrank == 0)
  {
    omg_i = malloc(iproc*jproc*sizeof(int));
    omg_j = malloc(iproc*jproc*sizeof(int));
    il = malloc(iproc*jproc*sizeof(int));
    ir = malloc(iproc*jproc*sizeof(int));
    jt = malloc(iproc*jproc*sizeof(int));
    jb = malloc(iproc*jproc*sizeof(int));
    rank_l = malloc(iproc*jproc*sizeof(int));
    rank_r = malloc(iproc*jproc*sizeof(int));
    rank_t = malloc(iproc*jproc*sizeof(int));
    rank_b = malloc(iproc*jproc*sizeof(int));

    init_parallel (iproc, jproc, imax, jmax, &myrank, il, ir, jb, jt, rank_l, rank_r, rank_b, rank_t, omg_i , omg_j, num_proc );
    /*printf("main after init p\n");
    printf("Rank in main %d %d %d %d \n", *rank_l, *rank_r, *rank_t, *rank_b);  //1st element
    printf("omg_i in main %d %d %d %d \n", omg_i[0], omg_i[1], omg_j[0], omg_j[1]);*/

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

    // Creating arrays for U,V,P RS, F and G   
    U = matrix (  0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 );
    V = matrix (  0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );
    P = matrix (  0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 );
    RS = matrix ( 0 , array_size[1] - array_size[0] , 0 , array_size[2] - array_size[3]);
    F = matrix (  0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] +1 );
    G = matrix (  0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );

    // Initializing the arrays U,V,P,RS,F and G
    init_matrix(U , 0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 , UI);
    init_matrix(V , 0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 , VI);
    init_matrix(P , 0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 , PI);
    init_matrix(RS , 0 , array_size[1] - array_size[0] , 0 , array_size[2] - array_size[3] , 0);
    init_matrix(F ,  0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] + 1 , 0);
    init_matrix(G , 0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] + 2 , 0);

    // initialize the time
    double t=0;

    // number of time steps
    int n = 0;   

    while (t<t_end)
    {
      boundary_values(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3], U , V ,F, G);

      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, array_size[1] - array_size[0] +2 , array_size[2] - array_size[3] + 2  , U, V, F, G, array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] );

      calculate_rs(dt,dx,dy, array_size[1] - array_size[0]  , array_size[2] - array_size[3] , F , G , RS);

      int it = 0;

      double res = 1000;

      while(it<itermax && res > eps) 
      {
        boundaryvalues_p(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3] , P);

        sor(omg, dx,dy,array_size[1] - array_size[0] ,array_size[2] - array_size[3] , P, RS, &res);

        pressure_comm(P, array_size[0], array_size[1], array_size[2],array_size[3], array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3], bufSend, bufRecv, &status, myrank);
        
        Program_Message("In the second while loop"); printf("%f ", res);
        
        double res_new ;

        MPI_Allreduce(&res , &res_new , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD ); 
        printf("%f ", res_new);
        res = res_new;

        it++;
      }

        calculate_uv( dt , dx , dy , array_size[1] - array_size[0] ,array_size[2] - array_size[3]  , U , V , F , G , P);

        uv_comm(U, V, array_size[0], array_size[1], array_size[2],array_size[3], array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3], bufSend, bufRecv, &status, myrank); 

        calculate_dt(Re , tau , &dt , dx , dy , array_size[1] - array_size[0] + 1, array_size[2] - array_size[3] +1 , U , V );

        double dt_new ;

        MPI_Allreduce(&dt , &dt_new , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

        dt = dt_new;

        t = t+dt;

        n = n+1;
    }

    //write_vtkFile("szProblem.vtk", n , xlength/2.0 , ylength/2.0 , array_size[1] - array_size[0], array_size[2] - array_size[3],dx, dy, U, V, P,myrank,array_pos[0],array_pos[1]);
    
    // Freeing all the momory allocated on the processor master processor.
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
    free_matrix(U , 0, array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 );
    free_matrix(V , 0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );
    free_matrix(P , 0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 );
    free_matrix(RS ,0, array_size[1] - array_size[0] , 0 ,array_size[2] - array_size[3] );
    free_matrix(F , 0, array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] +1);
    free_matrix(G ,0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2  );
    free(array_pos);
    free(array_size);
    free(array_neighbours);
  } 

// Running the process for all other threads except the master processer.
  else
  {
    //int* array_pos = malloc(2*sizeof(int));
    int* array_size = malloc(4*sizeof(int));
    int* array_neighbours = malloc(4*sizeof(int));

    // Reciving the size of the array along with its neighbours.
    MPI_Recv(array_neighbours,4,MPI_INT,0,2,MPI_COMM_WORLD,&status);    
    MPI_Recv(array_size,4,MPI_INT,0,3,MPI_COMM_WORLD,&status);

    // Creating arrays for U,V,P RS, F and G   
    U = matrix (  0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 );
    V = matrix (  0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );
    P = matrix (  0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 );
    RS = matrix ( 0 , array_size[1] - array_size[0] , 0 ,array_size[2] - array_size[3]);
    F = matrix (  0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] +1 );
    G = matrix (  0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );

    // Initializing the arrays U,V,P,RS,F and G
    init_matrix(U , 0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 , UI);
    init_matrix(V , 0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 , VI);
    init_matrix(P , 0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 , PI);
    init_matrix(RS , 0 , array_size[1] - array_size[0] , 0 ,array_size[2] - array_size[3] , 0);
    init_matrix(F ,  0 , array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] +1 , 0);
    init_matrix(G , 0 , array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 , 0);

    // initialize the time
    double t=0;

    // number of time steps
    int n = 0;

    while (t<t_end)
    {
    
      boundary_values(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3], U , V ,F, G);

      calculate_fg(Re,GX, GY, alpha, dt, dx, dy, array_size[1] - array_size[0] +2 , array_size[2] - array_size[3] + 2  , U, V, F, G, array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] );

      calculate_rs(dt,dx,dy, array_size[1] - array_size[0]  , array_size[2] - array_size[3] , F , G , RS);

      int it = 0;

      double res = 1000;

      while(it<itermax && res > eps) 
      {
        boundaryvalues_p(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3] , P);

        sor(omg, dx,dy,array_size[1] - array_size[0] ,array_size[2] - array_size[3] , P, RS, &res);

        pressure_comm(P, array_size[0], array_size[1], array_size[2],array_size[3], array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3], bufSend, bufRecv, &status, myrank);

        double res_new ;

        MPI_Allreduce(&res , &res_new , 1 , MPI_DOUBLE , MPI_SUM , MPI_COMM_WORLD );

        res =  res_new;

        it++;
      }

        calculate_uv( dt , dx , dy , array_size[1] - array_size[0] ,array_size[2] - array_size[3]  , U , V , F , G , P);

        uv_comm(U, V, array_size[0], array_size[1], array_size[2],array_size[3], array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3], bufSend, bufRecv, &status, myrank); 

        calculate_dt(Re , tau , &dt , dx , dy , array_size[1] - array_size[0] + 1, array_size[2] - array_size[3] +1 , U , V );

        double dt_new ;

        MPI_Allreduce(&dt , &dt_new , 1 , MPI_DOUBLE , MPI_MIN , MPI_COMM_WORLD );

        dt = dt_new;

        t = t+dt;

        n = n+1;
    }

    // write_vtkFile("szProblem.vtk", n , xlength/2.0 , ylength/2.0, array_size[1] - array_size[0], array_size[2] - array_size[3],dx, dy, U, V, P,myrank);
    
    // Freeing all the momory allocated on the processor.
    free(array_size);
    free(array_neighbours);
    free_matrix(U , 0, array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 );
    free_matrix(V , 0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );
    free_matrix(P , 0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 );
    free_matrix(RS ,0, array_size[1] - array_size[0] , 0 ,array_size[2] - array_size[3] );
    free_matrix(F , 0, array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] +1);
    free_matrix(G ,0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2  );
  }

  Programm_Sync("Synchronizing all the processors");
  Programm_Stop("Stoping Parallel Run");
  return 0;
}               