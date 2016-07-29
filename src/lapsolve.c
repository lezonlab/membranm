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

  Finds eigenvalues and eigenvectors of a matrix using LAPACK:  Takes an 
  nn-by-nn matrix, 'MM', a 'neig'-component vector 'VAL' and a nn-by-neig 
  array 'VEC'.  Places the 'neig' smallest eigenvalues in 'VAL' and their 
  corresponding eigenvectors in 'VEC'.
*/
#include "lapsolve.h"

void lapsolve(double **MM,double *VAL,double **VEC,int nn,int neig)
{
  int maxa,i,j,k;

  /* LAPACK variables */
  char *jobz="Vectors",*uplo="Upper";
  int len_jobz=1,len_uplo=1;
  double *A;
  double *W,*WORK;
  double wkopt;
  int lwork,info;


  /* Allocate memory in FFF (FORTRAN-friendly-format): A(I,J)=A[n*(J-1)+I-1]*/
  maxa=nn*nn;
  A=(double*)malloc(maxa*sizeof(double));
  if(!A){
    fprintf(stderr,"\nfailure allocating A in lapsolve\n\n");
    exit(1);}
  for(i=0;i<maxa;i++) A[i]=0.0;
  for(i=1;i<=nn;i++)
    for(j=1;j<=nn;j++)
      A[nn*(j-1)+i-1]=MM[i][j];
  W=(double*)malloc(nn*sizeof(double));
  if(!W){
    fprintf(stderr,"\nfailure allocating W in lapsolve\n\n");
    exit(1);}


  /* Call 'dsyev' once to find the optimal block size,
   and allocate memory for the WORK array */
  lwork = -1;
  dsyev_(jobz,uplo,&nn,A,&nn,W,&wkopt,&lwork,&info);
  lwork = (int)wkopt;
  WORK = (double*)malloc( lwork*sizeof(double) );
  if(!WORK){
    fprintf(stderr,"\nfailure allocating WORK in lapsolve\n\n");
    exit(1);}


  /* Solve eigenproblem */
  dsyev_(jobz,uplo,&nn,A,&nn,W,WORK,&lwork,&info);

  /* Put the solution into the output matrices */
  for(j=0;j<neig;j++){
    VAL[j+1]=W[j];
    for(i=0;i<nn;i++)
      VEC[i+1][j+1]=A[nn*j+i];
  }

  free(A);
  free(W);
  free(WORK);
}
