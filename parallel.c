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

   for (int i = 0;i < iproc ; i++){
         omg_i[i] = i*d1+d1;
         printf("%d\n",omg_i[i] );
         if (d2-- > 0)
            omg_i[i]++;
       }
   for (int j = 0;j < jproc ; j++){
         omg_j[j] = j*d3+d3;
         printf("%d\n",omg_j[j] );
         if (d4-- > 0)
            omg_j[j]++;
      }


  int k =0;
  for (int i = 0 ; i < iproc ; i++)
  {
    for (int j = 0; j < jproc; j++)
    {

      // finding the size of each block and saving them at the correct location.
      if (i==0){
         il[k] = 0;
      }
      else{ 
        il[k] = omg_i[i-1] ;
      }
      
      ir[k] = omg_i[i];
      
      if (j==0){
         jb[k] = 0;      
      }
      else{
         jb[k] = omg_j[j-1];
      }
      
      jt[k] = omg_j[j];


      printf("1st %d %d %d %d\n",il[k], ir[k] ,jt[k], jb[k] );
      int s[4] = {il[k], ir[k] ,jt[k], jb[k]};
      
      // finding the neighbours of each block and saving them at relevent positions.  
      if(i == 0){
        rank_l[k] = MPI_PROC_NULL;
      }
      else {
        rank_l[k] = jproc*(i-1)+j;
      }

      if (i == iproc-1){
        rank_r[k] = MPI_PROC_NULL;
      }
      else{
        rank_r[k] = jproc*(i+1)+j;
      }

      if(j == 0){
        rank_b[k] = MPI_PROC_NULL;
      }
      else {
        rank_b[k] = jproc*(i)+j-1;
      }

      if (j == jproc-1){
        rank_t[k] = MPI_PROC_NULL;
      }
      else{
        rank_t[k] = jproc*(i)+j+1;
      }
      int n[4] = {rank_t[k], rank_l[k], rank_r[k], rank_b[k]};
      int p[2] = {omg_i[i], omg_j[j]};

      printf("%d %d %d %d\n",rank_l[k], rank_r[k] ,rank_t[k], rank_b[k] );
      k = k+1;
      if (k<num_proc)
      {
         MPI_Send(n, 4, MPI_INT, k, 2, MPI_COMM_WORLD);
         MPI_Send(s, 4, MPI_INT, k, 3, MPI_COMM_WORLD);
         MPI_Send(p, 2, MPI_INT, k, 1, MPI_COMM_WORLD);
      }
      
    }
              printf("hurrY.............................................................\n"); 
  }

  }

void pressure_comm(double **P, int il, int ir , int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double* bufRecv, MPI_Status *status, int chunk)
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

   bufSend = (double*) malloc((biggest+2)*sizeof(double));
   bufRecv = (double*) malloc((biggest+2)*sizeof(double));

   if (rank_l != MPI_PROC_NULL)
   {
      int i =1;
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = P[i][j];
      }
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, status);

      i =0;
      for (int j = 0;j<=jt-jb+1;j++)
      {
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
      
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_r  , 0 ,MPI_COMM_WORLD, status);
      
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
      MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, status);
      
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
       MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, status);
      j = 1;
      for (int i = 0;i<=ir-il+1;i++)
      {
          P[i][j] = bufRecv[i];
      }
   }
   free(bufSend);
   free(bufRecv);
  
}



void uv_comm(double **U,double **V, int il, int ir , int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double* bufRecv, MPI_Status *status, int chunk)
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
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, status);

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
    
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_r  , 0 ,MPI_COMM_WORLD, status);
      
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
  
      MPI_Sendrecv(bufSend, ir-il+3, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+3 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, status);
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
      
      MPI_Sendrecv(bufSend, ir-il+3, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+3 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, status);
      j = 1;
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
      MPI_Sendrecv(bufSend, jt-jb+3, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+3 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, status);

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
      
      MPI_Sendrecv(bufSend, jt-jb+3, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+3 ,MPI_DOUBLE, rank_r  , 0 ,MPI_COMM_WORLD, status);
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
      MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, status);
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
       MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, status);
      j = 2;
      for (int i = 0;i<=ir-il+1;i++)
      {
          V[i][j] = bufRecv[i];
      }
   }
   free(bufSend);
   free(bufRecv);
}

