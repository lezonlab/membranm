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

  Driver for lapsolve
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "membranmutil.h"
#include "lapsolve.h"

char *extnsn(char *name,char *ext);
void sparsedim(char *file,int *row,int *col);

int main(int argc,char *argv[])
{
  FILE *data;
  double **HH,*VAL,**VEC;
  char *garp,*valfile,*vecfile;
  int neig;
  int row,col,i,j;


  /* Check for proper command-line format */
  if(argc!=2 && argc!=3){
    fprintf(stderr,"\nUsage:\n\n");
    fprintf(stderr,"\t%s mtxfile.sparsemtx [neigs]\n\n",argv[0]);
    fprintf(stderr,
	    "\tmtxfile.sparsemtx -- Sparse-format matrix to be decomposed\n");
    fprintf(stderr,"\tneigs -- Optional number of eigenvalues/vectors\n\n");
    fprintf(stderr,"\n%s: specify hessian matrix\n\n",argv[0]);
    fprintf(stderr,"Output: Prints eigensystem to mtxfile.val/mtxfile.vec\n\n");
    exit(1);}


  /* Read the matrix */
  sparsedim(argv[1],&row,&col);
  if(row!=col){
    fprintf(stderr,"\n%s: matrix '%s' is not symmetric\n\n",argv[0],argv[1]);
    exit(1);}
  HH=dmatrix(1,row,1,col);
  fprintf(stderr,"\nMatrix is %d by %d\n\n",row,col);
  read_sparsemtx(argv[1],HH,row,col);


  /* Allocate memory for eigensystem */
  if(argc==3) sscanf(argv[2],"%d",&neig);
  else neig=row;
  VAL=dvector(1,neig);
  VEC=dmatrix(1,row,1,neig);


  /* -------------------- Solve -------------------- */
  lapsolve(HH,VAL,VEC,row,neig);


  /* Print the output */
  garp=(char *)calloc(99,sizeof(char));
  strncpy(garp,argv[1],strlen(argv[1])-strlen(strrchr(argv[1],'.'))+4);
  valfile=extnsn(garp,"val");
  vecfile=extnsn(garp,"vec");
  fprintf(stderr,"valfile is %s\n",valfile);
  fprintf(stderr,"vecfile is %s\n",vecfile);
  data=fopen(valfile,"w");
  for(i=1;i<=neig;i++)
    fprintf(data,"% 16.7e\n",VAL[i]);
  fclose(data);
  data=fopen(vecfile,"w");
  for(i=1;i<=row;i++){
    for(j=1;j<=neig;j++)
      fprintf(data,"% 16.7e",VEC[i][j]);
    fprintf(data,"\n");
  }
  fclose(data);
  return 0;
}


/* "extnsn" REMOVES THE LAST THREE CHARACTERS FROM THE 
   FILENAME AND APPENDS IT WITH THE EXTENSION PROVIDED */
char *extnsn(char *name,char *ext)
{
  static char *garp;

  garp=(char *)calloc(99,sizeof(char));
  strncpy(garp,name,strlen(name)-3);
  strcat(garp,ext);
  return garp;
}


/* "sparsedim" RETURNS THE DIMENSIONS OF A SPARSE MATRIX */
void sparsedim(char *file,int *row,int *col)
{
  FILE *data;
  int i,j;
  double x;

  (*row)=(*col)=0;
  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nsparsedim: unable to open %s\n\n",file);
    exit(1);}
  while(!feof(data)){
    fscanf(data,"%d%d%lf",&i,&j,&x);
    if(i>(*row)) *row=i;
    if(j>(*col)) *col=j;
  }
  fclose(data);
}
