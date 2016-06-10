#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues_p(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int it,int ib,double** P);

void boundaryvalues_u(int rank_l, int rank_r, int rank_t, int rank_b ,int il, int ir, int it,int ib,double** V);

void boundaryvalues_v(int rank_l, int rank_r, int rank_t, int rank_b , int il, int ir, int it,int ib,double** V);

#endif
