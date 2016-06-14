#include "boundary_val.h"
#include <mpi.h>
#include "parallel.h"

// Boundary values for pressure..

void boundary_values(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int jt,int jb,double** U,double** V , double** F, double** G)
{
	//boundaryvalues_p(rank_l , rank_r , rank_t , rank_b , il , ir, jt , jb , P);
	boundaryvalues_u(rank_l , rank_r , rank_t , rank_b , il , ir, jt , jb , U);
	boundaryvalues_v( rank_l , rank_r , rank_t , rank_b , il , ir, jt , jb , V);
	boundaryvalues_F(rank_l , rank_r , rank_t , rank_b , il , ir, jt , jb , F , U);
	boundaryvalues_G(rank_l , rank_r , rank_t , rank_b , il , ir, jt , jb , G , V);
}

void boundaryvalues_p(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int jt,int jb,double** P)
{

	if (rank_l == MPI_PROC_NULL )
	{
		int i = 0; //left
		for (int j = 0 ; j<=jt-jb+1 ; j++)
				P[i][j] = P[i+1][j];
		
	}

	if (rank_r == MPI_PROC_NULL )
	{
		int i = ir-il+1; // right
		for (int j = 0 ; j<=jt-jb+1 ; j++)
				P[i][j] = P[i-1][j];
		
	}

	
	if (rank_t == MPI_PROC_NULL )
	{
		int j = jt-jb+1; //top
		for (int i = 0 ; i<=ir-il+1 ; i++)
				P[i][j] = P[i][j-1];

	}


	if (rank_b == MPI_PROC_NULL )
	{
		int j = 0;  // bottom
		for (int i = 0 ; i<=ir-il+1 ; i++)
				P[i][j] = P[i][j+1];

	}

}


//Boundary values for U velocjty

void boundaryvalues_u(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int jt,int jb,double** U)
{
	if (rank_l == MPI_PROC_NULL )
	{
		int i = 1;
		for (int j = 0 ; j<jt-jb +2 ; j++){
				U[i][j] = 0;
			    U[i-1][j] =0;}
		}
	


	if (rank_r == MPI_PROC_NULL )
	{
		int i = ir-il+1;
		
			for (int j = 0 ; j<jt-jb+2 ; j++){
				U[i][j] = 0;
			    U[i+1][j] =0;
			}
		
	}

	if (rank_t == MPI_PROC_NULL )
	{
		int j = jt-jb+1;
		
			for (int i = 0 ; i<ir-il+3 ; i++)
			    U[i][j] = 2.0 - U[i][j-1];
		
	}


	if (rank_b == MPI_PROC_NULL )
	{
		int j = 0;
		
			for (int i = 0 ; i< ir-il+3 ; i++)
				U[i][j] = -U[i][j+1];
		
	}

}


//Boundary values for V velocjties

void boundaryvalues_v(int rank_l, int rank_r, int rank_t, int rank_b , int il, int ir, int jt,int jb,double** V)
{
	if (rank_l == MPI_PROC_NULL )
	{
		int i = 0;
		
			for (int j = 0 ; j<=jt-jb+2 ; j++)
			V[i][j] = -V[i+1][j];
		
	}
	
	if (rank_r == MPI_PROC_NULL )
	{
		int i = ir-il+1;
		
			for (int j = 0 ; j<=jt-jb+2 ; j++)
			V[i][j] = -V[i-1][j];
		
	}

	if (rank_t == MPI_PROC_NULL )
	{
		int j = jt-jb+1;
		
			for (int i = 0 ; i<=ir-il+1 ; i++){
			V[i][j] = 0;
			V[i][j+1] = 0;}
		
	}

	if (rank_b == MPI_PROC_NULL )
	{
		int j = 0;
		
			for (int i = 0 ; i<=ir-il+1 ; i++){
			V[i][j] = 0;
			V[i][j+1] = 0;}
		
	}
}


void boundaryvalues_F(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int jt,int jb,double** F, double **U)
{
	if (rank_l == MPI_PROC_NULL )
	{
		int i = 1;
		for (int j = 0 ; j<=jt-jb+1 ; j++){
				F[i][j] = U[i][j];
			    F[i-1][j] = U[i][j];}
		}
	


	if (rank_r == MPI_PROC_NULL )
	{
		int i = ir-il+1;
		
			for (int j = 0 ; j<=jt-jb+1 ; j++){
				F[i][j] = U[i][j];
			    F[i+1][j] = U[i+1][j];}
		
	}
}

void boundaryvalues_G(int rank_l, int rank_r, int rank_t, int rank_b , int il, int ir, int jt,int jb,double** G, double** V)
{
	if (rank_t == MPI_PROC_NULL )
	{
		int j = jt-jb+1;
		
			for (int i = 0 ; i<=ir-il+1 ; i++){
			G[i][j] = V[i][j];
			G[i][j+1] = V[i][j+1];}
		
	}

	if (rank_b == MPI_PROC_NULL )
	{
		int j = 0;
		
			for (int i = 0 ; i<=ir-il+1 ; i++){
			G[i][j] = V[i][j];
			G[i][j+1] = V[i][j+1];}
		
	}
}
