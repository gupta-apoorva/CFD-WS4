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
   
   int* array_size = malloc(2*sizeof(int));
   int* array_pos = malloc(4*sizeof(int));
   int* array_neighbours = malloc(4*sizeof(int));

   omg_i = malloc(iproc*sizeof(int));
   omg_j = malloc(jproc*sizeof(int));

   int d1 = imax/iproc;
   int d2 = imax%iproc;
   int d3 = jmax/jproc;
   int d4 = jmax%jproc;

   
Program_Message("statting" );
   for (int i = 0;i < iproc ; i++)
      {
         omg_i[i] = d1;
         printf("%d\n",omg_i[i] );
         if (d2-- > 0)
            omg_i[i]++;

      }
   for (int j = 0;j < jproc ; j++)
      {
         omg_j[j] = d3;
         printf("%d\n",omg_j[j] );
         if (d4-- > 0)
            omg_j[j]++;
      }

  for (int i = 0 ; i<iproc ; i++)
  {
    for (int j = 0; j < jproc; ++j)
    {

      // finding the size of each block and saving them at the correct location.
      if (i==0)
      {
         int dummy = 0;
         il = &dummy;
         

      }

      else
      {
         *il = omg_i[i-1] ;
      }
      
      ir = (omg_i+i);
      
      if (j==0)
      {
         int dummy = 0;
         jb = &dummy;
          
      }
      else
      {
         jb = (omg_j+j-1);
      }
      
      jt = (omg_j+j);

      printf("%d %d %d %d\n",*il, *ir ,*jt, *jb );
      
      // finding the neighbours of each block and saving them at relevent positions.  
      if(i == 0){
         
         int dummy = MPI_PROC_NULL;
        rank_l = &dummy;
        printf("rank_l = %d\n",*rank_l);
     }

      else {
         int dummy = i-1;
        rank_l = &dummy;}

      if (i == iproc-1){
         int dummy = MPI_PROC_NULL;
        rank_r = &dummy;}
      else{
         int dummy = i+1;
        rank_r = &dummy;}

      if(j == 0){
         int dummy = MPI_PROC_NULL;
        rank_b = &dummy;}
      else {
         int dummy = j-1;
        rank_b = &dummy;}


      if (j == jproc-1){
         int dummy = MPI_PROC_NULL;
        rank_t = &dummy;}
      else{
         int dummy = j+1;
        rank_t = &dummy;}

      printf("%d %d %d %d\n",*rank_l, *rank_r ,*rank_t, *rank_b );

      *(array_pos+0) = i;
      *(array_pos+1) = j;

      *(array_size+0) = *il;
      *(array_size+1) = *ir;
      *(array_size+2) = *jb;
      *(array_size+3) = *jt;

      *(array_neighbours+0) = *rank_l;
      *(array_neighbours+1) = *rank_r;
      *(array_neighbours+2) = *rank_b;
      *(array_neighbours+3) = *rank_t;

      MPI_Send(&array_pos, 2, MPI_INT, (i+j)%num_proc , 1, MPI_COMM_WORLD);
      MPI_Send(&array_size, 4, MPI_INT, (i+j)%num_proc , 2, MPI_COMM_WORLD);
      MPI_Send(&array_neighbours, 4, MPI_INT, (i+j)%num_proc , 3, MPI_COMM_WORLD);

    }
  }
    free(omg_i);
    free(omg_j);
    free(array_pos);
    free(array_size);
    free(array_neighbours);
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

