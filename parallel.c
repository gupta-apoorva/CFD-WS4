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
   
   //int* array_size = malloc(4*sizeof(int));
   int* array_pos = malloc(2*sizeof(int));
   int* array_neighbours = malloc(4*sizeof(int));
   int *array_size = NULL;
   //int *array_pos = NULL;
   //int *array_neighbours = NULL;
   int array_size1[4] ;
   //int array_pos1[2] = NULL;
   //int *array_neighbours = NULL;

   int dummy5;
   int dummy6;
   omg_i = malloc(iproc*sizeof(int));
   omg_j = malloc(jproc*sizeof(int));

   int d1 = imax/iproc;
   int d2 = imax%iproc;
   int d3 = jmax/jproc;
   int d4 = jmax%jproc;

   
Program_Message("statting" );
   for (int i = 0;i < iproc ; i++)
      {
         omg_i[i] = i*d1+d1;
         printf("%d\n",omg_i[i] );
         if (d2-- > 0)
            omg_i[i]++;

      }
   for (int j = 0;j < jproc ; j++)
      {
         omg_j[j] = j*d3+d3;
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
         dummy5 = 0;
         il = &dummy5;
         

      }

      else
      {
         
         il = (omg_i+i-1) ;
      }
      
      ir = (omg_i+i);
      
      if (j==0)
      {
         dummy6 = 0;
         jb = &dummy6;
          
      }
      else
      {
         jb = (omg_j+j-1);
      }
      
      jt = (omg_j+j);
      int dummy1;
      int dummy2;
      int dummy3;
      int dummy4;

      printf("1st %d %d %d %d\n",*il, *ir ,*jt, *jb );
      
      // finding the neighbours of each block and saving them at relevent positions.  
      if(i == 0){
         
        dummy1 = MPI_PROC_NULL;
        rank_l = &dummy1;
        //printf("%d\n",*rank_l );
     }

      else {
         dummy1 = jproc*(i-1)+j;
        rank_l = &dummy1;
        //printf("%d\n",*rank_l );
     }


      if (i == iproc-1){
         dummy2 = MPI_PROC_NULL;
        rank_r = &dummy2;
      //printf("%d\n",*rank_r );
     }
      else{
         dummy2 = jproc*(i+1)+j;
        rank_r = &dummy2;
      //printf("%d\n",*rank_r );
     }
      //printf("rank_l = %d\n",*rank_l );
      if(j == 0){
        dummy3 = MPI_PROC_NULL;
        rank_b = &dummy3;
      //printf("%d\n",*rank_b );
     }
      else {
         dummy3 = jproc*(i)+j-1;
        rank_b = &dummy3;
      //printf("%d\n",*rank_b );
     }


      if (j == jproc-1){
         dummy4 = MPI_PROC_NULL;
        rank_t = &dummy4;
      //printf("%d\n",*rank_t );
     }
      else{
        dummy4 = jproc*(i)+j+1;
        rank_t = &dummy4;
        //printf("%d\n",*rank_t );
     }

      printf("%d %d %d %d\n",*rank_l, *rank_r ,*rank_t, *rank_b );

      *(array_pos+0) = i;
      *(array_pos+1) = j;
      printf("array_pos = %d\n", *(array_pos+0));

      /**(array_size+0) = *il;
      *(array_size+1) = *ir;
      *(array_size+3) = *jb;
      *(array_size+2) = *jt;*/
      printf("2nd %d %d %d %d\n", *il, *ir, *jt, *jb);
      array_size1[0] = *il;
      array_size1[1] = *ir;
      array_size1[3] = *jb;
      array_size1[2] = *jt;
      array_size = array_size1;

      printf("array_size il = %d\n", array_size[0]);
      printf("array_size ir = %d\n", array_size[1]);
      printf("array_size jt = %d\n", array_size[2]);
      printf("array_size jb = %d\n", array_size[3]);

      *(array_neighbours+0) = *rank_l;
      *(array_neighbours+1) = *rank_r;
      *(array_neighbours+2) = *rank_b;
      *(array_neighbours+3) = *rank_t;

      MPI_Send(&array_size, 4, MPI_INT, (i+j)%num_proc , 2, MPI_COMM_WORLD);
      printf("MPI_Send\n");
      MPI_Send(&array_pos, 2, MPI_INT, (i+j)%num_proc , 1, MPI_COMM_WORLD);
      printf("MPI_2ndsned\n");
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

