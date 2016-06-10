#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundary_p_left(int il,int ir,int it,int ib,double** P);

void boundary_p_right(int il,int ir,int it,int ib,double** P);

void boundary_p_top(int il,int ir,int it,int ib,double** P);

void boundary_p_bottom(int il,int ir,int it,int ib,double** P);

void boundary_u_left(int il,int ir,int it,int ib,double** U);

void boundary_u_right(int il,int ir,int it,int ib,double** U);

void boundary_u_top(int il,int ir,int it,int ib,double** U);

void boundary_u_bottom(int il,int ir,int it,int ib,double** U);

void boundary_v_left(int il,int ir,int it,int ib,double** V);

void boundary_v_right(int il,int ir,int it,int ib,double** V);

void boundary_v_top(int il,int ir,int it,int ib,double** V);

void boundary_v_bottom(int il,int ir,int it,int ib,double** V);

#endif
