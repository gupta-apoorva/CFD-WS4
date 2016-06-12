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
            
// Initializing MPI
          	MPI_Init(&argn, &args);
          	MPI_Comm_size( MPI_COMM_WORLD, &num_proc ); 
          	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);


if (myrank == 0){

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
        
        //setting the parameters
         read_parameters( "problem.dat", &Re , &UI , &VI, &PI, &GX, &GY, &t_end, &xlength, &ylength, &dt, &dx, &dy, &imax, &jmax, &alpha, &omg, &tau,&itermax, 
         &eps, &dt_value,&iproc,&jproc);
        
        MPI_Bcast(&Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&UI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&VI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&GX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&GY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&imax, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&jmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);

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
        //Program_Message("we are here \n \n \n \n");
// initialize the time
        double t=0;
// number of time steps
        int n = 0;   


  while (t<t_end)
   {

        calculate_dt(Re , tau , &dt , dx , dy , array_size[1] - array_size[0] + 1, array_size[2] - array_size[3] +1 , U , V );
        boundary_values(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3], U , V , P);
        calculate_fg(Re,GX, GY, alpha, dt, dx, dy, array_size[1] - array_size[0] , array_size[2] - array_size[3] , U, V, F, G);
                Program_Message("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
        calculate_rs(dt,dx,dy, array_size[1] - array_size[0]  , array_size[2] - array_size[3] , F , G , RS);

        int it = 0;
        double res = 1000;

        while(it<itermax && res > eps) 
        {
            sor(omg, dx,dy,array_size[1] - array_size[0] + 1,array_size[2] - array_size[3] +1 , P, RS, &res);
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

        calculate_uv( dt , dx , dy , array_size[1] - array_size[0] + 1,array_size[2] - array_size[3] +1 , U , V , F , G , P);
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
      	
else{
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

        MPI_Bcast(&Re, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&UI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&VI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&PI, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&GX, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&GY, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&t_end, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xlength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ylength, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dx, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&imax, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&jmax, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&alpha, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&omg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&tau, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&itermax, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&eps, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&dt_value, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&iproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&jproc, 1, MPI_INT, 0, MPI_COMM_WORLD);
        //printf("printing value of dt............. %f\n", dt );

// Creating the arrays U,V and P
        //int* array_pos = malloc(2*sizeof(int));
        int* array_size = malloc(4*sizeof(int));
        int* array_neighbours = malloc(4*sizeof(int));

        
        MPI_Recv(array_neighbours,4,MPI_INT,0,2,MPI_COMM_WORLD,&status);
       // printf("processes %d \n", myrank);
        //printf("Neighbours %d %d %d %d\n", *(array_neighbours + 0),*(array_neighbours + 1),*(array_neighbours + 2),*(array_neighbours + 3));
        MPI_Recv(array_size,4,MPI_INT,0,3,MPI_COMM_WORLD,&status);
        //printf("Size %d %d %d %d\n", *(array_size + 0),*(array_size + 1),*(array_size + 2),*(array_size + 3));
        //MPI_Recv(array_pos,2,MPI_INT,0,1,MPI_COMM_WORLD,&status);
        //printf("Position %d %d \n", *(array_pos + 0),*(array_pos + 1));

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
        //calculate_dt(Re , tau , &dt , dx , dy , array_size[1] - array_size[0] + 1, array_size[2] - array_size[3] +1 , U , V );
        boundary_values(array_neighbours[0], array_neighbours[1], array_neighbours[2], array_neighbours[3] ,array_size[0], array_size[1], array_size[2],array_size[3], U , V , P);
        
        calculate_fg(Re,GX, GY, alpha, dt, dx, dy, array_size[1] - array_size[0] , array_size[2] - array_size[3]  , U, V, F, G);
        calculate_rs(dt,dx,dy, array_size[1] - array_size[0]  , array_size[2] - array_size[3]  , F , G , RS);
                Program_Message("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
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
        free_matrix(U , 0, array_size[1] - array_size[0] + 2 ,  0 , array_size[2] - array_size[3] +1 );
        free_matrix(V , 0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2 );
        free_matrix(P , 0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +1 );
        free_matrix(RS ,0, array_size[1] - array_size[0] , 0 ,array_size[2] - array_size[3] );
        free_matrix(F , 0, array_size[1] - array_size[0] + 2 ,  0 , array_size[2] + array_size[3] +1);
        free_matrix(G ,0, array_size[1] - array_size[0] + 1 ,  0 , array_size[2] - array_size[3] +2  );

        Programm_Sync("Synchronizing all the processors");
        Programm_Stop("Stoping Parallel Run");
        return 0;
}
}
