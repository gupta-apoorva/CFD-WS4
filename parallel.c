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
   int array_size[4];
   int array_pos[2];
   int array_neighbours[4];
   omg_i = malloc(iproc*sizeof(int));
   omg_j = malloc(jproc*sizeof(int));

   int d1 = imax/iproc;
   int d2 = imax%iproc;
   int d3 = jmax/jproc;
   int d4 = jmax%jproc;
   for (int i = 0;i < iproc ; i++)
      {
         omg_i[i] = d1;
         if (d2-- > 0)
            omg_i[i]++;
      }
   for (int j = 0;j < jproc ; j++)
      {
         omg_j[j] = d3;
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
         il = 0;
      }
      else
      {
         il = omg_i[i-1] ;
      }
      ir = omg_i[i];
      if (j==0)
      {
         ib = 0;d   
      }
      else
      {
         ib = omg_j[j-1] +1;
      }
      it = omg_j[j];
      
      // finding the neighbours of each block and saving them at relevent positions.  
      if(i == 0)
        rank_l = MPI_PROC_NULL;
      else 
        rank_l = i-1;

      if (i == iproc-1)
        rank_r = MPI_PROC_NULL;
      else
        rank_r = i+1;

      if(j == 0)
        rank_b = MPI_PROC_NULL;
      else 
        rank_b = j-1;

      if (j == jproc-1)
        rank_t = MPI_PROC_NULL;
      else
        rank_t = j+1; 

      array_pos[0] = i;
      array_pos[1] = j;

      array_size[0] = il;
      array_size[1] = ir;
      array_size[2] = jb;
      array_size[3] = jt;

      array_neighbours[0] = rank_l;
      array_neighbours[1] = rank_r;
      array_neighbours[2] = rank_b;
      array_neighbours[3] = rank_t;



      MPI_Send(array_pos, 2, MPI_INT, (i+j)%num_proc , 0, MPI_COMM_WORLD);
      MPI_Send(array_size, 4, MPI_INT, (i+j)%num_proc , 0, MPI_COMM_WORLD);
      MPI_Send(array_neighbours, 4, MPI_INT, (i+j)%num_proc , 0, MPI_COMM_WORLD);

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
      for (int i =1)   
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = P[i][j];
      }
      }
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, &status);

      for (int i =0)
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
          P[i][j] = bufRecv[j];
      }
      }
   }
  else
   {
    void boundary_p_left(int il,int ir,int it,int ib,double** P);
   }

   if (rank_r != MPI_PROC_NULL)
   {
      for (int i =ir-il)
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = P[i][j];
      }
      }
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_r  , 0 ,MPI_COMM_WORLD, &status);
      for (int i = ir-il+1)
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
          P[i][j] = bufRecv[j];
      }
      }
   }
   else
   {
    void boundary_p_right(int il,int ir,int it,int ib,double** P);
   }

   if (rank_t != MPI_PROC_NULL)
   {
      for (int j = it-ib)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = P[i][j];
      }
      }
      MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, &status);
      for (int j = it-ib+1)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
          P[i][j] = bufRecv[i];
      }
      }
   }
   else
   {
    void boundary_p_top(int il,int ir,int it,int ib,double** P);
   }

   if (rank_b != MPI_PROC_NULL)
   {
      for (int j =1)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = P[i][j];
      }
      }
       MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, &status);
       for (int j = 1)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
          P[i][j] = bufRecv[i];
      }
      }
   }
   else
   {
    void boundary_p_bottom(int il,int ir,int it,int ib,double** P);
   }
}



void uv_comm(double **U,double **V, int il, int ir , int jb, int jt, int rank_l, int rank_r, int rank_b, int rank_t, double *bufSend, double* bufRecv, MPI_Status *status, int chunk);
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
      for (int i =2)   
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = U[i][j];
      }
      }
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, &status);

      for (int i =0)
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
          U[i][j] = bufRecv[j];
      }
      }
   }
  else
   {
    void boundary_u_left(int il,int ir,int it,int ib,double** U);
   }

   if (rank_r != MPI_PROC_NULL)
   {
      for (int i =ir-il)
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
         bufSend[j] = U[i][j];
      }
      }
      MPI_Sendrecv(bufSend, jt-jb+2, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+2 ,MPI_DOUBLE, rank_r  , 0 ,MPI_COMM_WORLD, &status);
      for (int i = ir-il+2)
      {
      for (int j = 0;j<=jt-jb+1;j++)
      {
          U[i][j] = bufRecv[j];
      }
      }
   }
   else
   {
    void boundary_u_right(int il,int ir,int it,int ib,double** U);
   }

   if (rank_t != MPI_PROC_NULL)
   {
      for (int j = it-ib)
      {
      for (int i = 0;i<=ir-il+2;i++)
      {
         bufSend[i] = U[i][j];
      }
      }
      MPI_Sendrecv(bufSend, ir-il+3, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+3 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, &status);
      for (int j = it-ib+1)
      {
      for (int i = 0;i<=ir-il+2;i++)
      {
          U[i][j] = bufRecv[i];
      }
      }
   }
   else
   {
    void boundary_u_top(int il,int ir,int it,int ib,double** U);
   }

   if (rank_b != MPI_PROC_NULL)
   {
      for (int j =1)
      {
      for (int i = 0;i<=ir-il+2;i++)
      {
         bufSend[i] = U[i][j];
      }
      }
       MPI_Sendrecv(bufSend, ir-il+3, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+3 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, &status);
       for (int j = 1)
      {
      for (int i = 0;i<=ir-il+2;i++)
      {
          U[i][j] = bufRecv[i];
      }
      }
   }
   else
   {
    void boundary_u_bottom(int il,int ir,int it,int ib,double** U);
   }

// communicating the values of V

   if (rank_l != MPI_PROC_NULL)
   {
      for (int i =1)   
      {
      for (int j = 0;j<=jt-jb+2;j++)
      {
         bufSend[j] = V[i][j];
      }
      }
      MPI_Sendrecv(bufSend, jt-jb+3, MPI_DOUBLE, rank_l, 0, bufRecv, jt-jb+3 ,MPI_DOUBLE, rank_l , 0 ,MPI_COMM_WORLD, &status);

      for (int i =0)
      {
      for (int j = 0;j<=jt-jb+2;j++)
      {
          V[i][j] = bufRecv[j];
      }
      }
   }
  else
   {
    void boundary_v_left(int il,int ir,int it,int ib,double** V);
   }

   if (rank_r != MPI_PROC_NULL)
   {
      for (int i =ir-il)
      {
      for (int j = 0;j<=jt-jb+2;j++)
      {
         bufSend[j] = V[i][j];
      }
      }
      MPI_Sendrecv(bufSend, jt-jb+3, MPI_DOUBLE, rank_r, 0, bufRecv, jt-jb+3 ,MPI_DOUBLE, rank_r  , 0 ,MPI_COMM_WORLD, &status);
      for (int i = ir-il+1)
      {
      for (int j = 0;j<=jt-jb+2;j++)
      {
          V[i][j] = bufRecv[j];
      }
      }
   }
   else
   {
    void boundary_v_right(int il,int ir,int it,int ib,double** V);
   }

   if (rank_t != MPI_PROC_NULL)
   {
      for (int j = it-ib)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = V[i][j];
      }
      }
      MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_t, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_t  , 0 ,MPI_COMM_WORLD, &status);
      for (int j = it-ib+2)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
          V[i][j] = bufRecv[i];
      }
      }
   }
   else
   {
    void boundary_v_top(int il,int ir,int it,int ib,double** V);
   }

   if (rank_b != MPI_PROC_NULL)
   {
      for (int j =2)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
         bufSend[i] = V[i][j];
      }
      }
       MPI_Sendrecv(bufSend, ir-il+2, MPI_DOUBLE, rank_b, 0, bufRecv, ir-il+2 ,MPI_DOUBLE, rank_b  , 0 ,MPI_COMM_WORLD, &status);
       for (int j = 2)
      {
      for (int i = 0;i<=ir-il+1;i++)
      {
          V[i][j] = bufRecv[i];
      }
      }
   }
   else
   {
    void boundary_v_bottom(int il,int ir,int it,int ib,double** V);
   }

}

