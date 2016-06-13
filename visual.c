#include "helper.h"
#include "visual.h"
#include <stdio.h>


void write_vtkFile(const char *szProblem,
		 int    timeStepNumber,
		 double xlength,
                 double ylength,
                 int    imax,
                 int    jmax,
		 double dx,
		 double dy,
                 double **U,
                 double **V,
                 double **P,
		int il,
int ir,
int jt,
int jb,
int myrank) {
  
  int i,j;
  char szFileName[80];
  FILE *fp=NULL;
  sprintf( szFileName, "%s.%i.%d.vtk", szProblem, timeStepNumber,myrank );
  fp = fopen( szFileName, "w");
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to open %s", szFileName );
    ERROR( szBuff );
    return;
  }

  write_vtkHeader( fp, imax, jmax, dx, dy);
  write_vtkPointCoordinates(fp, imax, jmax, dx, dy, il ,  ir,  jt, jb);

  fprintf(fp,"POINT_DATA %i \n", (imax+1)*(jmax+1) );
	
  fprintf(fp,"\n");
  fprintf(fp, "VECTORS velocity float\n");
  for(j = 0; j < jmax+1; j++) {
    for(i = 0; i < imax+1; i++) {
      fprintf(fp, "%f %f 0\n", (U[i][j] + U[i][j+1]) * 0.5, (V[i][j] + V[i+1][j]) * 0.5 );
    }
  }

  fprintf(fp,"\n");
  fprintf(fp,"CELL_DATA %i \n", ((imax)*(jmax)) );
  fprintf(fp, "SCALARS pressure float 1 \n"); 
  fprintf(fp, "LOOKUP_TABLE default \n");
  for(j = 1; j < jmax+1; j++) {
    for(i = 1; i < imax+1; i++) {
      fprintf(fp, "%f\n", P[i][j] );
    }
  }

  if( fclose(fp) )
  {
    char szBuff[80];
    sprintf( szBuff, "Failed to close %s", szFileName );
    ERROR( szBuff );
  }
}


void write_vtkHeader( FILE *fp, int imax, int jmax, 
                      double dx, double dy) {
  if( fp == NULL )		       
  {
    char szBuff[80];
    sprintf( szBuff, "Null pointer in write_vtkHeader" );
    ERROR( szBuff );
    return;
  }

  fprintf(fp,"# vtk DataFile Version 2.0\n");
  fprintf(fp,"generated by CFD-lab course output (written by Tobias Neckel) \n");
  fprintf(fp,"ASCII\n");
  fprintf(fp,"\n");	
  fprintf(fp,"DATASET STRUCTURED_GRID\n");
  fprintf(fp,"DIMENSIONS  %i %i 1 \n", imax+1, jmax+1);
  fprintf(fp,"POINTS %i float\n", (imax+1)*(jmax+1) );
  fprintf(fp,"\n");
}


void write_vtkPointCoordinates( FILE *fp, int imax, int jmax, 
                      double dx, double dy,int il , int ir, int jt, int jb) {
  int originX = il;  
  int originY = jb;
  int i ;
  int j ;

  for(j = 0; j <= jt-jb; j++) {
    for(i = 0; i <= ir-il; i++) {
      fprintf(fp, "%d %d 0\n", originX + i, originY + j);
    }
  }

  /*int i = 0;
  int j = 0;

  for(j = 0; j < jmax+1; j++) {
    for(i = 0; i < imax+1; i++) {
      fprintf(fp, "%f %f 0\n", originX+(i*dx), originY+(j*dy) );
    }
  }*/
}

void output_uvp(double **U,double **V,double **P,int il,int ir,
int jb,int jt,int omg_i,int omg_j,char *output_file){

}



void write_vtkPointCoordinates2( FILE *fp, int il,int ir,
int jb,int jt,int omg_i,int omg_j ) {
  //double originX = omg_i*(ir-il);  
  //double originY = omg_j*(jt-jb);


  int i ;
  int j ;

  for(j = il; j <= ir; j++) {
    for(i = jb; i <= jt; i++) {
      fprintf(fp, "%d %d 0\n", i, j);
    }
  }}



