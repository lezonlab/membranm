/*
  Licensed under the MIT License (MIT)

  Copyright (c) 2008-2015 Timothy Lezon

  This file is part of the imANM software package

  Permission is hereby granted, free of charge, to any person obtaining a 
  copy of this software and associated documentation files (the "Software"), 
  to deal in the Software without restriction, including without limitation 
  the rights to use, copy, modify, merge, publish, distribute, sublicense, 
  and/or sell copies of the Software, and to permit persons to whom the 
  Software is furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
  DEALINGS IN THE SOFTWARE.


  ------------------------------------------------------------------------
  Author: Tim Lezon

  This code multiplies the modes of a block Hessian by a projection matrix,
  giving the modes of the full system as approximated by the blocks.
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<ctype.h>
#include "membranmutil.h"

void read_proj_matrix(char *,dSparse_Matrix *,int);
void dsort_PP(dSparse_Matrix *,int);
int filecol(char *file);
void read_dvecs1(char *file,double **VEC,int nrow,int nmod);

int main(int argc,char *argv[])
{
  dSparse_Matrix PP;
  double **VEC,*MAG,dd;
  int nres,nmod,row,elm,kold,i,j,k,p;

  if(argc!=3){
    fprintf(stderr,"\n%s: %s%s\n\n",argv[0],"specify a .prj file ",
	    "and a .vec file for rigid blocks");
    exit(1);}
  nmod=filecol(argv[2]);
  row=filerow(argv[2]);
  elm=filerow(argv[1]);

  fprintf(stderr,"\nConverting %s using %s...\n",argv[2],argv[1]);
  fprintf(stderr,"%d modes\n~%d blocks\n%d non-zero elements\n",nmod,row/6,elm);

  PP.IDX=imatrix(1,elm,1,2);
  PP.X=dvector(1,elm);
  VEC=dmatrix(1,row,1,nmod);

  read_proj_matrix(argv[1],&PP,elm);
  dsort_PP(&PP,elm);
  read_dvecs1(argv[2],VEC,row,nmod); 


  nres=0;
  for(i=1;i<=elm;i++)
    if(PP.IDX[i][1]>nres)
      nres=PP.IDX[i][1];

  fprintf(stderr,"%d residues\n",nres/3);


  /* CALCULATE MAGNITUDE OF EACH VECTOR */
  MAG=dvector(1,nmod);
  for(i=1;i<=nmod;i++) MAG[i]=0.0;
  k=kold=1;
  for(i=1;i<=nres;i++){
    for(j=1;j<=nmod;j++){
      dd=0.0;
      k=kold;
      while(k<=elm && PP.IDX[k][1]==i){
	p=PP.IDX[k][2];
	dd+=PP.X[k]*VEC[p][j];
	k++;
      }
      MAG[j]+=(dd*dd);
    }
    kold=k;
  }
  for(i=1;i<=nmod;i++) MAG[i]=sqrt(MAG[i]);

  /* PRINT NORMALIZED VECTORS */
  k=kold=1;
  for(i=1;i<=nres;i++){
    for(j=1;j<=nmod;j++){
      dd=0.0;
      k=kold;
      while(k<=elm && PP.IDX[k][1]==i){
	p=PP.IDX[k][2];
	dd+=PP.X[k]*VEC[p][j];
	k++;
      }
      printf("%f\t",dd/MAG[j]);
    }
    kold=k;
    printf("\n");
  }
  return 0;
}
  

/* "read_proj_matrix" READS THE PROJECTION MATRIX FROM FILE */
void read_proj_matrix(char *file,dSparse_Matrix *PP,int elm)
{
  FILE *data;
  int i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_proj_matrix: unable to open %s\n\n",file);
    exit(1);}
  for(i=1;i<=elm;i++)
    fscanf(data,"%d%d%lf",&PP->IDX[i][1],&PP->IDX[i][2],&PP->X[i]);
  fclose(data);
}


/* "dsort_PP" SORTS THE PROJECTION MATRIX IN ASCENDING ORDER OF THE INDICES.  
   ADAPTED FROM THE NUMERICAL RECIPES 'HEAPSORT' ROUTINE. */
void dsort_PP(dSparse_Matrix *MM,int n)
{
  double x;
  int i,ir,j,l,hi,i1,i2;
  unsigned long rra,*ra;

  if(n<2) return;

  /* CREATE A VECTOR TO INDEX THE ELEMENTS OF MM */
  hi=0;
  for(i=1;i<=n;i++)
    if(MM->IDX[i][2]>hi)
      hi=MM->IDX[i][2];
  ra=lvector(1,n);
  for(i=1;i<=n;i++)
    ra[i]=(long)hi*(MM->IDX[i][1]-1)+MM->IDX[i][2];


  /* SORT */
  l=(n >> 1)+1;
  ir=n;
  for(;;){
    if(l > 1){
      rra=ra[--l];
      i1=MM->IDX[l][1];
      i2=MM->IDX[l][2];
      x=MM->X[l];
    } 
    else {
      rra=ra[ir];
      i1=MM->IDX[ir][1];
      i2=MM->IDX[ir][2];
      x=MM->X[ir];
      ra[ir]=ra[1];
      MM->IDX[ir][1]=MM->IDX[1][1];
      MM->IDX[ir][2]=MM->IDX[1][2];
      MM->X[ir]=MM->X[1];
      if (--ir == 1) {
	ra[1]=rra;
	MM->IDX[1][1]=i1;
	MM->IDX[1][2]=i2;
	MM->X[1]=x;
	break;
      }
    }
    i=l;
    j=l+l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
	ra[i]=ra[j];
	MM->IDX[i][1]=MM->IDX[j][1];
	MM->IDX[i][2]=MM->IDX[j][2];
	MM->X[i]=MM->X[j];
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
    MM->IDX[i][1]=i1;
    MM->IDX[i][2]=i2;
    MM->X[i]=x;
  }
  free_lvector(ra,1,n);
}


/* "filecol" RETURNS THE NUMBER OF COLUMNS IN 
   THE FIRST LINE OF THE SPECIFIED FILE */
int filecol(char *file)
{
  FILE *data;
  char c;
  int lc,col;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nfilecol: unable to open %s\n\n",file);
    exit(1);}
  lc=1;
  col=0;
  while((c=getc(data))!=EOF){
    if(!isspace(c)){
      if(lc==1)
	col++;
      lc=0;
    }
    else{
      lc=1;
      if(c=='\n'){
	fclose(data);
	return col;
      }
    }
  }
  return col;
}


/* "read_dvecs1" READS A LIST OF EIGENVECTORS FROM A FILE */
void read_dvecs1(char *file,double **VEC,int nrow,int nmod)
{
  FILE *data;
  int i,j;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_dvecs1: unable to open %s\n\n",file);
    exit(1);}
  for(i=1;i<=nrow;i++)
    for(j=1;j<=nmod;j++)
      fscanf(data,"%lf",&VEC[i][j]);
  fclose(data);
}
