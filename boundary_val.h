#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */

void boundary_values(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int jt,int jb,double** U,double** V , double** F, double** G);

void boundaryvalues_p(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int it,int ib,double** P);

void boundaryvalues_u(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int it,int ib,double** V);

void boundaryvalues_v(int rank_l, int rank_r, int rank_t, int rank_b , int il, int ir, int it,int ib,double** V);

void boundaryvalues_G(int rank_l, int rank_r, int rank_t, int rank_b , int il, int ir, int jt,int jb,double** G, double** V);

void boundaryvalues_F(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int jt,int jb,double** F, double **U);

#endif
