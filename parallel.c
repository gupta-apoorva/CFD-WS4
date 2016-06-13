#include "parallel.h"


void Program_Message(char *txt)
/* produces a stderr text output  */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
}


void Programm_Sync(char *txt)
/* produces a stderr textoutput and synchronize all processes */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                             /* synchronize output */  
   fprintf(stderr,"-MESSAGE- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
}


void Programm_Stop(char *txt)
/* all processes will produce a text output, be synchronized and finished */

{
   int myrank;

   MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
   MPI_Barrier(MPI_COMM_WORLD);                           /* synchronize output */
   fprintf(stderr,"-STOP- P:%2d : %s\n",myrank,txt);
   fflush(stdout);
   fflush(stderr);
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   exit(1);
}

void init_parallel (int iproc, int jproc, int imax, int jmax, int *myrank, int *il, int *ir, int *jb, int *jt, int *rank_l, int *rank_r, int *rank_b, int *rank_t, 
                   int *omg_i , int *omg_j, int num_proc )
{
   int d1 = imax/iproc;
   int d2 = imax%iproc;
   int d3 = jmax/jproc;
   int d4 = jmax%jproc;


    int x = *myrank/jproc;
    int y = *myrank - x*jproc;

    *omg_i = x;
    *omg_j = y;

    if(x == 0){
    *rank_l = MPI_PROC_NULL;
    }
    else {
    *rank_l = jproc*(x-1)+y;
    }

    if (x == iproc-1){
    *rank_r = MPI_PROC_NULL;
    }
    else{
    *rank_r = jproc*(x+1)+y;
    }

    if( y == 0){
    *rank_b = MPI_PROC_NULL;
    }
    else {
    *rank_b = jproc*(x)+y-1;
    }

    if (y == jproc-1){
    *rank_t = MPI_PROC_NULL;
    }
    else{
    *rank_t = jproc*(x)+y+1;
    }

   


      // finding the size of each block and saving them at the correct location.
      if (x==0){
         *il = 0;
      }
      else{ 
        *il = d1*x;
      }
      
      *ir = d1*(x+1);

      if(d2>x)
      { if (x == 0)
        *ir = *ir + x +1;
         else
        {
          *il = *il + x;
          *ir = *ir + x +1;
        }
       }
       else
       {
        *il = *il + d2;
        *ir = *ir + d2; 
       }

      
      if (y==0){
         *jb = 0;      
      }
      else{
         *jb = y*d3;
      }
      
      *jt = d3*(y+1);

      if(d4>y)
      { if (y ==0)
         *jt = *jt + y+1;
         else
        {
          *jb = *jb + y;
          *jt = *jt + y+1;
        }
      }
      else
      {
          *jb = *jb + d4;
          *jt = *jt + d4;
      }

  }

  

void pressure_comm(double **P, int il, int ir , int jt, int jb, int rank_l, int rank_r, int rank_t, int rank_b, double *bufSend, double* bufRecv, MPI_Status *status, int chunk)
{
  int x = ir - il +2;
  int y = jt - jb +2;

  if (rank_l != MPI_PROC_NULL)
   {
      int i =1;
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = P[i][j];
      }

      MPI_Sendrecv(bufSend, y, MPI_DOUBLE, rank_l, 21 , bufRecv, y ,MPI_DOUBLE, rank_r , 22 ,MPI_COMM_WORLD, status);


     i =0;
      for (int j = 0;j<=jt-jb+1;j++)
      {
        //printf("bufRecv[j]  %f \n" , bufRecv[j]);
          P[i][j] = bufRecv[j];
      }
   }

   if (rank_r != MPI_PROC_NULL)
   {
      int i =ir-il;
      for (int j = 0;j<=jt-jb+1;j++)
      {
          bufSend[j] = P[i][j];
      }
      
      MPI_Sendrecv(bufSend, y, MPI_DOUBLE, rank_r, 22 , bufRecv, y ,MPI_DOUBLE, rank_l  , 21 ,MPI_COMM_WORLD, status);;
      i = ir-il+1;
      for (int j = 0;j<=jt-jb+1;j++)
      {
          P[i][j] = bufRecv[j];
      }
   }
  

   if (rank_t != MPI_PROC_NULL)
   {
      int j = jt-jb;
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = P[i][j];
      }
      MPI_Sendrecv(bufSend, x, MPI_DOUBLE, rank_t, 0, bufRecv, x ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, status);
      
      j = jt-jb+1;
      for (int i = 0;i<=ir-il+1;i++)
      {
          P[i][j] = bufRecv[i];
      }
   }
   

   if (rank_b != MPI_PROC_NULL)
   {
      int j =1;
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = P[i][j];
      }
       MPI_Sendrecv(bufSend, x, MPI_DOUBLE, rank_b, 0, bufRecv, x ,MPI_DOUBLE, rank_t , 0 ,MPI_COMM_WORLD, status);
      j = 0;
      for (int i = 0;i<=ir-il+1;i++)
      {
          P[i][j] = bufRecv[i];
      }
   }
    
}



void uv_comm(double **U,double **V, int il, int ir , int jt, int jb, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double* bufRecv, MPI_Status *status, int chunk)
{
   int biggest;
   if (ir-il>jt-jb)
   {
      biggest = ir-il;
   }   
   else
   {
      biggest = jt-jb;
   }

   bufSend = (double*) malloc((biggest+3)*sizeof(double));
   bufRecv = (double*) malloc((biggest+3)*sizeof(double));


// Communicaing the values of U

   if (rank_l != MPI_PROC_NULL)
   {
      int i =2;
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = U[i][j];
      }
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_r , 0 ,MPI_COMM_WORLD, status);

      i =0;
      for (int j = 0;j<=jt-jb+1;j++)
      {
          U[i][j] = bufRecv[j];
      }
   }
  

   if (rank_r != MPI_PROC_NULL)
   {
      int i =ir-il;
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = U[i][j];
      }
    
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, status);
      
      i = ir-il+2;
      for (int j = 0;j<=jt-jb+1;j++)
      {
          U[i][j] = bufRecv[j];
      }
   }


   if (rank_t != MPI_PROC_NULL)
   {
      int j = jt-jb;
      for (int i = 0;i<=ir-il+2;i++)
      {
         bufSend[i] = U[i][j];
      }
  
      MPI_Sendrecv(bufSend, ir-il+3, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+3 ,MPI_DOUBLE, rank_b , 0 ,MPI_COMM_WORLD, status);
      j = jt-jb+1;
      for (int i = 0;i<=ir-il+2;i++)
      {
          U[i][j] = bufRecv[i];
      }
      
   }
   

   if (rank_b != MPI_PROC_NULL)
   {
      int j =1;
      for (int i = 0;i<=ir-il+2;i++)
      {
         bufSend[i] = U[i][j];
      }
      
      MPI_Sendrecv(bufSend, ir-il+3, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+3 ,MPI_DOUBLE, rank_t , 0 ,MPI_COMM_WORLD, status);
      j = 0;
      for (int i = 0;i<=ir-il+2;i++)
      {
          U[i][j] = bufRecv[i];
      }
   }


// communicating the values of V

   if (rank_l != MPI_PROC_NULL)
   {
      int i =1;
      for (int j = 0;j<=jt-jb+2;j++)
      {
         bufSend[j] = V[i][j];
      }
      MPI_Sendrecv(bufSend, jt-jb+3, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+3 ,MPI_DOUBLE, rank_r, 0 ,MPI_COMM_WORLD, status);

      i =0;
      for (int j = 0;j<=jt-jb+2;j++)
      {
          V[i][j] = bufRecv[j];
      }
   }

   if (rank_r != MPI_PROC_NULL)
   {
      int i =ir-il;
      for (int j = 0;j<=jt-jb+2;j++)
      {
         bufSend[j] = V[i][j];
      }
      
      MPI_Sendrecv(bufSend, jt-jb+3, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+3 ,MPI_DOUBLE, rank_l  , 0 ,MPI_COMM_WORLD, status);
      i = ir-il+1;
      for (int j = 0;j<=jt-jb+2;j++)
      {
          V[i][j] = bufRecv[j];
      }
   }

   if (rank_t != MPI_PROC_NULL)
   {
      int j = jt-jb;
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = V[i][j];
      }
      MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, status);
      j = jt-jb+2;
      for (int i = 0;i<=ir-il+1;i++)
      {
          V[i][j] = bufRecv[i];
      }
   }

   if (rank_b != MPI_PROC_NULL)
   {
      int j =2;
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = V[i][j];
      }
       MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, status);
      j = 0;
      for (int i = 0;i<=ir-il+1;i++)
      {
          V[i][j] = bufRecv[i];
      }
   }
   free(bufSend);
   free(bufRecv);
}

