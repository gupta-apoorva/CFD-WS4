#include "boundary_val.h"

// Boundary values for pressure..

void boundary_p_left(int il, int ir, int it,int ib,double** P)
{
	for (int i = 0)
	{
		for (int j = 0 ; j<=it-ib+1 ; j++)
			P[i][j] = P[i+1][j];
	}
}

void boundary_p_right(int il, int ir, int it,int ib,double** P)
{
	for (int i = ir-il+1)
	{
		for (int j = 0 ; j<=it-ib+1 ; j++)
			P[i][j] = P[i-1][j];
	}
}

void boundary_p_top(int il, int ir, int it,int ib,double** P)
{
	for (int j = it-ib+1)
	{
		for (int i = 0 ; i<=ir-il+1 ; i++)
			P[i][j] = P[i][j-1];
	}
}

void boundary_p_bottom(int il, int ir, int it,int ib,double** P)
{
	for (int j = 0)
	{
		for (int i = 0 ; i<=ir-il+1 ; i++)
			P[i][j] = P[i][j+1];
	}
}

//Boundary values for U velocity

void boundary_u_left(int il, int ir, int it,int ib,double** V)
{
	for (int i = 1)
	{
		for (int j = 0 ; j<=it-ib+1 ; j++)
			U[i][j] = 0;
		    U[i-1][j] =0;
	}
}

void boundary_u_right(int il, int ir, int it,int ib,double** V)
{
	for (int i = ir-il+1)
	{
		for (int j = 0 ; j<=it-ib+1 ; j++)
			U[i][j] = 0;
		    U[i+1][j] =0;
	}
}

void boundary_u_top(int il, int ir, int it,int ib,double** U)
{
	for (int j = it-ib+1)
	{
		for (int i = 0 ; i<=ir-il+2 ; i++)
		    U[i][j] = 2.0 - U[i][j-1];
	}
}

void boundary_u_bottom(int il, int ir, int it,int ib,double** U)
{
	for (int j = 0)
	{
		for (int i = 0 ; i<=ir-il+2 ; i++)
			U[i][j] = -U[i][j+1];
	}
}

//Boundary values for V velocities

void boundary_v_left(int il, int ir, int it,int ib,double** V)
{
	for (int i = 0)
	{
		for (int j = 0 ; j<=it-ib+2 ; j++)
			V[i][j] = -V[i+1][j];
	}
}

void boundary_v_right(int il, int ir, int it,int ib,double** V)
{
	for (int i = ir-il+1)
	{
		for (int j = 0 ; j<=it-ib+2 ; j++)
			V[i][j] = -V[i-1][j];
	}
}

void boundary_v_top(int il, int ir, int it,int ib,double** V)
{
	for (int j = it-ib+1)
	{
		for (int i = 0 ; i<=ir-il+1 ; i++)
			V[i][j] = 0;
			V[i+1][j] = 0;
	}
}

void boundary_v_bottom(int il, int ir, int it,int ib,double** V)
{
	for (int j = 0)
	{
		for (int i = 0 ; i<=ir-il+1 ; i++)
			V[i][j] = 0;
			V[i+1][j] = 0;
	}
}


void boundaryvalues(int imax,int jmax,double **U,double **V, double** P, double** G, double** F)
{

    
	for(int j=1; j<=jmax; j++)
	{
		U[0][j] = 0;
		U[imax][j] = 0;
        V[0][j] = -V[1][j];
		V[imax+1][j] = -V[imax][j];
        P[0][j] = P[1][j];
		P[imax+1][j] = P[imax][j];
        F[0][j] = U[0][j];
        F[imax][j] = U[imax][j];
	}

	for(int i=1; i<=imax; i++)
	{
				
		U[i][0] = -U[i][1];
		U[i][jmax+1] = 2.0 - U[i][jmax];
        V[i][0] = 0;
		V[i][jmax] = 0;
        P[i][0] = P[i][1];
		P[i][jmax+1] = P[i][jmax];
        G[i][0] = V[i][0];
        G[i][jmax] = V[i][jmax]; 
    }
}
