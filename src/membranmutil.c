/*
  Licensed under the MIT License (MIT)

  Copyright (c) 2008-2015 Timothy Lezon

  This file is part of the imANM and exANM software packages

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

  Utility functions for membrane ANM calculations.
*/

#include "membranmutil.h"
#include<stdio.h>
#include <stddef.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<ctype.h>
#define NR_END 1

/* "assign_rigid_blocks" ASSIGNS EACH RESIDUE IN THE PDB STRUCTURE 
   TO A RIGID BLOCK AND RETURNS THE TOTAL NUMBER OF RIGID BLOCKS */
int assign_rigid_blocks(PDB_File *PDB,Rigid_Block *BLX,
			int nres,int nrg,int *bmx)
{
  char ch;
  int *UB,**BMP,nb,nbt,ok,kwn,mold,rr,i,j,k;

  /* DETERIMINE THE NUMBER OF UNIQUE BLOCKS PER CHAIN */
  UB=ivector(1,nrg);
  for(i=1;i<=nrg;i++) UB[i]=0;
  nb=0;
  for(i=1;i<=nrg;i++){
    ok=1;
    for(j=1;j<=nb;j++)
      if(BLX[i].blknum==UB[j]){
	ok=0;
	break;
      }
    if(ok==1)
      UB[++nb]=BLX[i].blknum;
  }


  /* ASSIGN EACH RESIDUE TO A BLOCK: BMP[i][1] CONTAINS THE NUMBER OF 
     BLOCK i AS PROVIDED BY THE .blk FILE.  BMP[i][2] CONTAINS THE 
     NUMBER OF BLOCK i THAT WILL BE USED IN THE RTB CALCULATION. */
  BMP=imatrix(1,nb,1,2);
  nbt=kwn=0;
  for(i=1;i<=nb;i++)
    BMP[i][1]=BMP[i][2]=0;
  mold=PDB->atom[1].model;

  for(i=1;i<=nres;i++){
    rr=PDB->atom[i].resnum;
    ch=PDB->atom[i].chain;
    if(PDB->atom[i].model!=mold){
      for(j=1;j<=nb;j++)
	BMP[j][1]=BMP[j][2]=0;
      mold=PDB->atom[i].model;
      kwn=0;
    }
    ok=0;
    for(j=1;j<=nrg;j++){
      if((BLX[j].lochain<ch && ch<BLX[j].hichain) ||
	 (BLX[j].lochain==ch && ch<BLX[j].hichain && rr>=BLX[j].lonum) ||
	 (BLX[j].lochain<ch && ch==BLX[j].hichain && rr<=BLX[j].hinum) ||
	 (BLX[j].lochain==ch && BLX[j].hichain==ch && 
	  rr>=BLX[j].lonum && rr<=BLX[j].hinum)){
	for(k=1;k<=kwn;k++)
	  if(BLX[j].blknum==BMP[k][1]){
	    ok=1;
	    PDB->atom[i].model=BMP[k][2];
	    break;
	  }
	if(ok==0){
	  BMP[++kwn][1]=BLX[j].blknum;
	  BMP[kwn][2]=PDB->atom[i].model=++nbt;
	  ok=1;
	}
	break;
      }
    }
    if(ok==0)
      PDB->atom[i].model=0;
  }
  free_ivector(UB,1,nrg);
  free_imatrix(BMP,1,nb,1,nb);


  /* FIND THE SIZE OF THE LARGEST BLOCK */
  UB=ivector(1,nbt);
  for(i=1;i<=nbt;i++) UB[i]=0;
  for(i=1;i<=nres;i++) 
    if(PDB->atom[i].model!=0)
      UB[PDB->atom[i].model]++;
  (*bmx)=0;
  for(i=1;i<=nbt;i++){
    if(UB[i]>(*bmx))
      (*bmx)=UB[i];
  }
  free_ivector(UB,1,nbt);

  return nbt;
}



/* "bless_from_tensor" transfers the block Hessian 
   from  the tensor HT into the array HB */
int bless_from_tensor(double **HB,double ***HT,int **CT,int nblx)
{
  int *I1,*I2,i,j,p,sb,ii,jj,max,a,b,imx;

  
  max=6*nblx;
  I1=ivector(1,max);
  I2=ivector(1,max);

  /* I1[i]==i iff there is a non-zero element in column i 
     (removes zeroes that are caused by single-node blocks) */
  for(i=1;i<=max;i++){ 
    I1[i]=0;
    for(j=i;j<=max;j++)
      HB[i][j]=HB[j][i]=0.0;
  }
  for(ii=1;ii<=nblx;ii++){
    for(i=1;i<=6;i++){
      for(jj=ii;jj<=nblx;jj++){
	sb=CT[ii][jj];
	if(sb!=0){
	  p = jj==ii ? i : 1;
	  for(j=p;j<=6;j++)
	    if(HT[sb][i][j]!=0)
	      I1[6*(jj-1)+j]=6*(jj-1)+j;
	}
      }
    }
  }

  /* If I1[i]!=0, then I2[i] is a sequential index */
  imx=0;
  for(i=1;i<=max;i++){
    if(I1[i]!=0) imx++;
    I2[i]=imx;
  }

  for(ii=1;ii<=nblx;ii++){
    for(i=1;i<=6;i++){
      for(jj=ii;jj<=nblx;jj++){
	sb=CT[ii][jj];
	if(sb!=0){
	  p = jj==ii ? i : 1;
	  for(j=p;j<=6;j++)
	    if(HT[sb][i][j]!=0){
	      a=I2[6*(ii-1)+i];
	      b=I2[6*(jj-1)+j];
	      HB[a][b]=HB[b][a]=HT[sb][i][j];
	    }
	}
      }
    }
  }
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
  return imx;
}


/* "centroid_vector" ALLOCATES MEMORY FOR A 1D ARRAY OF MASS CENTROIDS */
Centroid *centroid_vector(int lo,int hi)
{
  static Centroid *v;
  int i;

  /* ALLOCATE ARRAY OF Centroids */
  v=(Centroid *)malloc((size_t) ((hi-lo+2)*sizeof(Centroid)));
  if (!v) nrerror("allocation failure in centroid_vector");
  v+=1;
  v-=lo;

  /* ALLOCATE A 3-VECTOR FOR COORDINATES OF EACH CENTROID */
  v[lo].X=(double *) malloc((size_t)(((hi-lo+1)*3+1)*sizeof(double)));
  if (!v[lo].X) nrerror("allocation failure 2 in centroid_vector");
  for(i=lo+1;i<=hi;i++)
    v[i].X=v[i-1].X+3;

  return v;
}


/* allocate a char matrix with subscript range m[nrl..nrh][ncl..nch] */
char **cmatrix(long nrl,long nrh,long ncl,long nch)
{
  long i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  static char **m;

  /* allocate pointers to rows */
  m=(char **) malloc((size_t)((nrow+1)*sizeof(char*)));
  if(!m){
    fprintf(stderr,"\nallocation failure 1 in cmatrix\n\n");
    exit(1);}
  m++;
  m-=nrl;

  /* allocate rows and set pointers to them */
  m[nrl]=(char *) malloc((size_t)((nrow*ncol+1)*sizeof(char)));
  if(!m[nrl]){
    fprintf(stderr,"\nallocation failure 2 in cmatrix\n\n");
    exit(1);}
  m[nrl]++;
  m[nrl]-=ncl;

  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

  /* return pointer to array of pointers to rows */
  return m;
}



/* "copy_dsparse" COPIES ELEMENTS lo THROUGH hi 
   OF SPARSE MATRIX 'A' TO SPARSE MATRIX 'B' */
void copy_dsparse(dSparse_Matrix *A,dSparse_Matrix *B,int lo,int hi)
{
  int i;

  for(i=lo;i<=hi;i++){
    B->IDX[i][1]=A->IDX[i][1];
    B->IDX[i][2]=A->IDX[i][2];
    B->X[i]=A->X[i];
  }
}


/* "cross" TAKES THE 3D CROSS PRODUCT OF ITS ARGUMENTS. */
void cross(double x[], double y[], double z[])
{
  z[0] = x[1]*y[2] - x[2]*y[1];
  z[1] = x[2]*y[0] - x[0]*y[2];
  z[2] = x[0]*y[1] - x[1]*y[0];
}


/* "dsort_PP2" SORTS THE PROJECTION MATRIX IN ASCENDING ORDER OF THE  
   INDEX 'idx'.  ADAPTED FROM THE NUMERICAL RECIPES 'HEAPSORT' ROUTINE. */
void dsort_PP2(dSparse_Matrix *MM,int n,int idx)
{
  double x;
  int i,ir,j,l,hi,i1,i2,ndx;
  unsigned long rra,*ra;

  if(n<2) return;
  if(idx<1 || idx>2){
    fprintf(stderr,"dsort_PP2: bad index value %d\n\n",idx);
    exit(1);}
  ndx = idx==1 ? 2 : 1;

  /* CREATE A VECTOR TO INDEX THE ELEMENTS OF MM */
  hi=0;
  for(i=1;i<=n;i++)
    if(MM->IDX[i][ndx]>hi)
      hi=MM->IDX[i][ndx];
  ra=lvector(1,n);
  for(i=1;i<=n;i++)
    ra[i]=(long)hi*(MM->IDX[i][idx]-1)+MM->IDX[i][ndx];


  /* SORT */
  l=(n >> 1)+1;
  ir=n;
  for(;;){
    if(l > 1){
      rra=ra[--l];
      i1=MM->IDX[l][idx];
      i2=MM->IDX[l][ndx];
      x=MM->X[l];
    } 
    else {
      rra=ra[ir];
      i1=MM->IDX[ir][idx];
      i2=MM->IDX[ir][ndx];
      x=MM->X[ir];
      ra[ir]=ra[1];
      MM->IDX[ir][idx]=MM->IDX[1][idx];
      MM->IDX[ir][ndx]=MM->IDX[1][ndx];
      MM->X[ir]=MM->X[1];
      if (--ir == 1) {
	ra[1]=rra;
	MM->IDX[1][idx]=i1;
	MM->IDX[1][ndx]=i2;
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
	MM->IDX[i][idx]=MM->IDX[j][idx];
	MM->IDX[i][ndx]=MM->IDX[j][ndx];
	MM->X[i]=MM->X[j];
	i=j;
	j <<= 1;
      } else j=ir+1;
    }
    ra[i]=rra;
    MM->IDX[i][idx]=i1;
    MM->IDX[i][ndx]=i2;
    MM->X[i]=x;
  }
  free_lvector(ra,1,n);
}


/* "filerow" RETURNS THE NUMBER OF NEWLINES IN THE SPECIFIED FILE */
int filerow(char *file)
{
  FILE *data;
  char c;
  int row=0;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nfilecol: can't open %s\n\n",file);
    exit(1);}
  while((c=getc(data))!=EOF)
    if(c=='\n')
      row++;
  fclose(data);
  return row;
}



/* "find_contacts" FINDS WHICH BLOCKS ARE IN CONTACT, AND ASSIGNS EACH 
   PAIR OF CONTACTING BLOCKS A UNIQUE INDEX.  IT RETURNS THE TOTAL NUMBER 
   OF CONTACTS BETWEEN BLOCKS. */
int find_contacts1(int **CT,PDB_File *PDB,int nres,int nblx,double cut)
{
  int nc,i,j,k,ii,jj;
  double csq=cut*cut,df,dd;

  for(i=1;i<=nres;i++){
    ii=PDB->atom[i].model;
    for(j=i+1;j<=nres;j++){
      jj=PDB->atom[j].model;
      
      if(ii!=jj && ii!=0 && jj!=0 && CT[ii][jj]==0){
	dd=0.0;
	for(k=0;k<3;k++){
	  df=(double)PDB->atom[i].X[k]-PDB->atom[j].X[k];
	  dd+=df*df;
	}
	if(dd<csq)
	  CT[ii][jj]=CT[jj][ii]=1;
      }
      
    }
  }
  nc=0;
  for(i=1;i<=nblx;i++)
    for(j=i;j<=nblx;j++)
      if(CT[i][j]!=0){
	nc++;
	CT[i][j]=CT[j][i]=nc;
      }
  return nc;
}


/* free a char matrix allocated by cmatrix() */
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
{
  free((FREE_ARG) (m[nrl]+ncl-1));
  free((FREE_ARG) (m+nrl-1));
}


/* "get_param" RETURNS THE PARAMETER OF THE 
   SPECIFIED NAME IN THE SPECIFIED FILE */
char *get_param(char *file,char *param)
{
  FILE *data;
  char LINE[200],*s;
  static char *garp;

  garp=(char *)calloc(99,sizeof(char));
  strcpy(garp,param);
  strcat(garp,"=");
  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nget_param: unable to open %s\n\n",file);
    exit(1);}
  while(!feof(data)){
    fgets(LINE,200,data);
    if(strstr(LINE,garp)!=NULL && LINE[0]!='#'){
      s=strpbrk(LINE,"=")+1;
      strcpy(garp,s);
      garp[strlen(garp)-1]='\0';
      return garp;
    }
  }
  return NULL;
}



/* "init_bst" INITIALIZES THE n-COMPONENT VECTOR 'BST': GIVEN THE 'elm' 
   ELEMENT SPARSE MATRIX 'PP', SORTED BY INDEX 'idx', INITIALIZES 'BST' 
   SUCH THAT FOR ALL j: BST[i]<=j<BST[i+1], PP->IDX[j][idx]=i */
void init_bst(int *BST,dSparse_Matrix *PP,int elm,int n,int idx)
{
  int i;

  for(i=1;i<n;i++) BST[i]=0;
  for(i=elm;i>0;i--) BST[PP->IDX[i][idx]]=i;
  BST[n]=elm+1;
  for(i=n-1;i>0;i--)
    if(BST[i]==0)
      BST[i]=BST[i+1];
}


/* "pdb_field_count" READS THE NUMBER OF TIMES THAT A 
   LINE BEGINNING WITH 'str' OCCURS IN THE PDB FILE */
int pdb_field_count(char *file,char *str)
{
  FILE *data;
  char c,LINE[PDB_MAX_LINE],HED[7];
  int i,n=0;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\npdb_field_count: unable to open %s\n\n",file);
    exit(1);} 

  for(;;){
    i=0;
    do{                                  /* Skip remaining columns */
      c=getc(data);
      LINE[i++]=c;
    }while(c!='\n' && !feof(data));
    LINE[i]='\0';
    for(i=0;i<=5;i++)                       /* Cols 1-6 are field name */
      HED[i]=LINE[i];
    HED[6]='\0';
    if(!strncmp(HED,str,6))
      n++;
    else if((!strncmp(HED,"END",3) && strncmp(HED,"ENDMDL",6)) || feof(data)){
      fclose(data);
      return n;
    }
  }
  fclose(data);
  return n;
}


/* "pdb_hmr" ESTIMATES THE NUMBER OF HEADER 
   LINES, MODELS AND RESIDUES IN A PDB FILE */
void pdb_hmr(char *file,int *nhed,int *nmod,int *nca)
{
  FILE *data;
  char HED[7],ATM[5],c,calt;
  int i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\npdb_hmr: unable to open %s\n\n",file);
    exit(1);}
  (*nca)=(*nhed)=(*nmod)=0;
  for(;;){

    /* Cols 1-6 are field name */
    for(i=0;i<6;i++)
      HED[i]=getc(data);

    /* Check for ATOM lines */
    if(!strncmp(HED,"ATOM  ",6)){  
      for(i=7;i<=12;i++) c = getc(data);
      for(i=13;i<=16;i++) ATM[i-13] = getc(data);
      calt = getc(data);
      if(strstr(ATM,"CA")!=NULL && (calt=='A' || calt==' '))
	(*nca)++;
    }

    /* Check for lines to be included in the header */
    else if(!strncmp(HED,"HEADER",6) || !strncmp(HED,"TITLE ",6) || 
	    !strncmp(HED,"COMPND",6) || !strncmp(HED,"SEQRES",6) || 
	    !strncmp(HED,"HELIX ",6) || !strncmp(HED,"SHEET ",6) || 
	    !strncmp(HED,"TURN  ",6))
      (*nhed)++;

    /* Check for rotation matrices */
    else if(!strncmp(HED,"MODEL ",6))
      (*nmod)++;

    else if(!strncmp(HED,"END   ",6) || feof(data))  /* Identify END of file */
      break;
    do{
      c=getc(data);
    }while(c!='\n');
  }
  fclose(data);
  if(*nmod==0)
    *nmod=1;
}



/* "pdb_init" INITIALIZES A PDB_File STRUCTURE: 'HEADER' IS A CHARACTER
   ARRAY WITH ROWS (1,nhed) AND COLUMNS (0,PDB_MAX_LINE-1); Calpha IS
   AN ARRAY OF Atom_Line STRUCTURES WITH ELEMENTS (1,nres).  */
void pdb_init(PDB_File *PDB,int nhed,int nres)
{

  /* Allocate memory for PDB.HEADER */
  PDB->HEADER=cmatrix(1,nhed,0,PDB_MAX_LINE-1);

  /* Allocate memory for PDB.atom */
  PDB->atom=malloc((size_t)((nres+2)*sizeof(Atom_Line)));
  if(!PDB->atom){
    fprintf(stderr,"\npdb_init: fail to allocate atom\n\n");
    exit(1);}
}



/* "print_prj_ofst" PRINTS THE PROJECTION MATRIX, OFFSETTING THE COLUMN 
   ELEMENTS IF THERE ARE BLOCKS WITH FEWER THAN SIX DEGREES OF FREEDOM */
void print_prj_ofst(dSparse_Matrix *PP,int elm)
{
  int *I1,*I2,max=0,i,j=0;

  for(i=1;i<=elm;i++)
    if(PP->IDX[i][2]>max)
      max=PP->IDX[i][2];
  I1=ivector(1,max);
  I2=ivector(1,max);
  for(i=1;i<=max;i++) I1[i]=0;
  for(i=1;i<=elm;i++)
    I1[PP->IDX[i][2]]=PP->IDX[i][2];
  for(i=1;i<=max;i++){
    if(I1[i]!=0) j++;
    I2[i]=j;
  }
  fprintf(stderr,"Printing projection matrix...\n");
  for(i=1;i<=elm;i++)
    if(PP->X[i]!=0.0)
      printf("%8d%8d% 25.15e\n",PP->IDX[i][1],I2[PP->IDX[i][2]],PP->X[i]);
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}


/* "read_blockfile" READS THE INPUT blk FILE FOR 
   INFORMATION ON THE SOURCE PDB AND THE BLOCK STRUCTURE */
void read_blockfile(char *file,char **pdbfile,Rigid_Block **BLK,int *nrg)
{
  FILE *data;
  Rigid_Block *Btmp;
  char *ptmp,HED[7],BNK[8],c;
  int blk=0,rng,nn,i;
  long int blkst;


  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_blockfile: unable to open %s\n\n",file);
    exit(1);}
  for(;;){

    /* GET THE FIRST CHARACTER OF EACH LINE */
    c=getc(data);


    /* COMMENT OR WHITESPACE: SKIP OVER REST OF LINE */
    if(c=='#' || isspace(c))
      while(c!='\n' && c!=EOF) c = getc(data);


    /* OTHERWISE CHECK FOR READ ROW HEADING */
    else{
      i=0;
      while(!isspace(c) && i<6){
	HED[i++]=c;
	c=getc(data);
      }
      HED[i] = '\0';


      /* PDB FILENAME */
      if(!strncmp(HED,"PDB",3)){
	i=0;
	do{
	  c=getc(data);
	  i++;
	} while(c!='\n' && c!=EOF);

	/* Allocate space for filename, then go back and read it */
	ptmp=(char *)malloc((size_t)(i+1)*sizeof(char));
	fseek(data,-i,SEEK_CUR);
	i=0;
	do{
	  c = getc(data);
	  if(!isspace(c)) ptmp[i++] = c;
	}while(c!='\n' && c!=EOF);
	ptmp[i] = '\0';
	*pdbfile = ptmp;
      }


      /* BLOCK DATA */
      else if(!strncmp(HED,"BLOCK",5)  && blk==0){
	blk=1;
	nn=0;
	blkst=ftell(data)-strlen(HED)-1;

	/* Count the number of ranges */
	while(!strncmp(HED,"BLOCK",5)){
	  while(c!='\n' && c!=EOF)
	    c=getc(data);
	  c=getc(data);
	  if(c!='#' && !isspace(c)){
	    i=0;
	    while(c!='\n' && c!=EOF && i<6){
	      HED[i++]=c;
	      c=getc(data);
	    }
	    HED[i]='\0';
	    nn++;
	  }
	}


	/* Allocate space for ranges, then go back and read the data */
	Btmp=(Rigid_Block *)malloc((size_t)(nn+1)*sizeof(Rigid_Block));
	fseek(data,blkst,SEEK_SET);


	/* --------------- READ ALL THE BLOCKS --------------- */
	rng=1;
	do{
	  c=getc(data);
	  if(c!='#' && !isspace(c) && c!=EOF){
	    i=0;
	    while(!isspace(c) && i<6){
	      HED[i++]=c;
	      c=getc(data);
	    }
	    HED[i]='\0';
	    if(!strncmp(HED,"BLOCK",5)){
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      sscanf(BNK,"%d",&Btmp[rng].blknum);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      strncpy(Btmp[rng].LORES,BNK,3);
	      while(isspace(c)) c=getc(data);
	      Btmp[rng].lochain=c;
	      c=getc(data);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      sscanf(BNK,"%d",&Btmp[rng].lonum);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      strncpy(Btmp[rng].HIRES,BNK,3);
	      while(isspace(c)) c=getc(data);
	      Btmp[rng].hichain=c;
	      c=getc(data);
	      while(isspace(c)) c=getc(data);
	      i=0;
	      while(!isspace(c)){
		BNK[i++]=c;
		c=getc(data);
	      }
	      BNK[i]='\0';
	      sscanf(BNK,"%d",&Btmp[rng].hinum);
	      while(c!='\n' && c!=EOF) c=getc(data);
	      rng++;
	    }
	  }
	  else{
	    while(c!='\n' && c!=EOF) c=getc(data);
	    if(c==EOF){
	      fclose(data);
	      *BLK=Btmp;
	      *nrg=nn;
	      return;
	    }
	  }
	}while(!strncmp(HED,"BLOCK",5));
      }

      /* END OF FILE */
      if(!strncmp(HED,"END",3)){
	fclose(data);
	*BLK=Btmp;
	*nrg=nn;
	return;
      }

      /* EAT UP REST OF LINE */
      while(c!='\n' && c!=EOF) c=getc(data);
      if(c==EOF){
	fclose(data);
	*BLK=Btmp;
	*nrg=nn;
	return;
      }
    } /* <------- END of else{...*/
  } /* <--------- END of for(;;){...*/
  fclose(data);
  free(ptmp);
}


/* "read_centfile" READS THE CENTROIDS FROM 'file', GETS 
   THEIR NAMES AND MASSES FROM 'map', AND STORES THEM IN 'PTS' */
int read_centfile(char *file,Centroid *PTS,Map *map,int np,int nmp)
{
  FILE *data;
  char LINE[PDB_MAX_LINE],HED[7],HOLD[9],s[9],c;
  int ca,ok,i,j;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_centfile: unable to open %s\n\n",file);
    exit(1);}

  ca=1;
  for(;;){
    
    i=0;
    do{                                  /* Skip remaining columns */
      c=getc(data);
      LINE[i++]=c;
    }while(c!='\n' && !feof(data));
    LINE[i]='\0';
    if(i<5){
      fclose(data);
      return 0;
    }
    for(i=0;i<=5;i++)                       /* Cols 1-6 are field name */
      HED[i]=LINE[i];
    HED[6]='\0';


    /* ---------- Check for ATOM or HETATM lines ---------- */
    if(!strncmp(HED,"ATOM  ",6) || !strncmp(HED,"HETATM",6)){

      for(i=6;i<=10;i++)                    /* Cols 7-11 are atom # */
	HOLD[i-6]=LINE[i];
      HOLD[5]='\0';
      sscanf(HOLD,"%d",&PTS[ca].num);
      for(i=30;i<=37;i++)                   /* Cols 31-38 are x-coordinate */
	HOLD[i-30]=LINE[i];
      HOLD[8]='\0';
      sscanf(HOLD,"%lf",&PTS[ca].X[1]);
      for(i=38;i<=45;i++)                   /* Cols 39-46 are y-coordinate */
	HOLD[i-38]=LINE[i];
      HOLD[8]='\0';
      sscanf(HOLD,"%lf",&PTS[ca].X[2]);
      for(i=46;i<=53;i++)                   /* Cols 47-54 are z-coordinate */
	HOLD[i-46]=LINE[i];
      HOLD[8]='\0';
      sscanf(HOLD,"%lf",&PTS[ca].X[3]);
      for(i=12;i<=15;i++)                    /* Cols 13-16 are atom name */
	HOLD[i-12]=LINE[i];
      HOLD[i-12]='\0';
      
      /* FIND PSEUDONYM AND ASSIGN MASS */
      ok=0;
      sscanf(HOLD,"%s",s);
      for(j=1;j<=nmp;j++)
	if(!strcmp(s,map[j].name1)){
	  sprintf(PTS[ca].name,"%s",map[j].name2);
	  PTS[ca].mass=map[j].mass;
	  ok=1;
	  break;
	}
      if(ok==0){
	fprintf(stderr,"\nread_centfile: no mass for %s in\n%s\n\n",HOLD,LINE);
	exit(1);}
      ca++;
    }

    /* ---------- Check for end of file ---------- */
    else if((!strncmp(HED,"END",3) && strncmp(HED,"ENDMDL",6)) || feof(data)){
      fclose(data);
      return 0;
    }
  }
}


/* "read_centroids1" RETURNS A VECTOR TO 'np' CENTROIDS 
   WITH THE PROPERTIES IN THE ARGUMENT FILES. */
Centroid *read_centroids1(char *pdbfile,char *mapfile,int *np)
{
  Map *map;
  static Centroid *PTS;
  int nmp;

  /* COUNT THE NUMBER OF PROTEINS TO BE NAMED */
  nmp=filerow(mapfile)-1;


  /* CREATE A STRUCTURE TO HOLD PROTEIN NAMES AND MASSES */
  map=(Map *)malloc((size_t) ((nmp+1)*sizeof(Map)));


  /* READ THE MAP */
  read_map(mapfile,map,nmp);


  /* COUNT THE NUMBER OF ATOM AND HETATM FIELDS */
  (*np)=pdb_field_count(pdbfile,"ATOM  ");
  (*np)+=pdb_field_count(pdbfile,"HETATM");


  /* INITIALIZE A Centroid ARRAY OF THE CORRECT SIZE */
  PTS=centroid_vector(1,*np);


  /* READ HETATM DATA FROM PDB FILE, SUBSTITUTING PROTEIN 
     NAMES FOR STAND-INS AND INCLUDING THE CENTROID MASSES */
  read_centfile(pdbfile,PTS,map,*np,nmp);


  free(map);
  return PTS;
}


/* "read_map" READS THE MAP FROM THE FILE INTO THE Map STRUCTURE OF SIZE np */
int read_map(char *file,Map *map,int np)
{
  FILE *data;
  char c;
  int i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_map: unable to open %s\n\n",file);
    exit(1);}
  do c=getc(data); 
  while(c!='\n');
  for(i=1;i<=np;i++)
    fscanf(data,"%s%s%lf",map[i].name1,map[i].name2,&map[i].mass);
  fclose(data);
  return 0;
}


/* "read_pdb3" READS HEADER AND Calpha INFORMATION FROM A PDB FILE */
int read_pdb3(char *file,PDB_File *PDB,int nhed,int nres)
{
  FILE *data;
  char LINE[PDB_MAX_LINE],HED[7],ATM[5],HOLD[9],c,calt;
  int ca,hd,mdl,i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_pdb3: can't open %s\n\n",file);
    exit(1);}

  ca=hd=mdl=1;
  for(;;){
    
    i=0;
    do{                                  /* Skip remaining columns */
      c=getc(data);
      LINE[i++]=c;
    }while(c!='\n' && !feof(data));
    LINE[i]='\0';
    for(i=0;i<=5;i++)                       /* Cols 1-6 are field name */
      HED[i]=LINE[i];
    HED[6]='\0';



    /* ---------- Check for ATOM lines ---------- */
    if(!strncmp(HED,"ATOM  ",6)){/* || !strncmp(HED,"HETATM",6)){ */

      for(i=0;i<=5;i++) PDB->atom[ca].HEAD[i]=LINE[i];
      PDB->atom[ca].HEAD[6]='\0';

      for(i=12;i<=15;i++)                     /* Cols 13-16 are atom name*/
	ATM[i-12]=LINE[i];
      ATM[4]='\0';
      calt=LINE[16];                          /* Col 17 is alt. location */

      /* Keep atom if it is alpha carbon */
      if(strstr(ATM,"CA")!=NULL && (calt=='A' || calt==' ')){

	for(i=6;i<=10;i++)                    /* Cols 7-11 are atom # */
	  HOLD[i-6]=LINE[i];
	HOLD[5]='\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].atmnum);

	for(i=12;i<=15;i++)                     /* Cols 13-16 are atom name*/
	  PDB->atom[ca].ATOM[i-12]=LINE[i];
	PDB->atom[ca].ATOM[4]='\0';

	for(i=17;i<=19;i++)                   /* Cols 18-20 are res. name */
	  PDB->atom[ca].RES[i-17]=LINE[i];
	PDB->atom[ca].RES[3]='\0';
	PDB->atom[ca].chain=LINE[21];       /* Col 22 is chain name */

	/* Assign a model number */
	PDB->atom[ca].model=mdl;

	for(i=22;i<=25;i++)                   /* Cols 23-26 are residue # */
	  HOLD[i-22]=LINE[i];
	HOLD[4]='\0';
	sscanf(HOLD,"%d",&PDB->atom[ca].resnum);
	for(i=30;i<=37;i++)                   /* Cols 31-38 are x-coordinate */
	  HOLD[i-30]=LINE[i];
	HOLD[8]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[0]);
	for(i=38;i<=45;i++)                   /* Cols 39-46 are y-coordinate */
	  HOLD[i-38]=LINE[i];
	HOLD[8]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[1]);
	for(i=46;i<=53;i++)                   /* Cols 47-54 are z-coordinate */
	  HOLD[i-46]=LINE[i];
	HOLD[8]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].X[2]);
	for(i=60;i<=65;i++)              /* Columns 61-66 are beta */
	  HOLD[i-60]=LINE[i];
	HOLD[6]='\0';
	sscanf(HOLD,"%f",&PDB->atom[ca].beta);
	for(i=76;i<=77;i++)              /* Columns 77-78 are element */
	  PDB->atom[ca].ELEMENT[i-76] = LINE[i]=='\n' ? '\0' : LINE[i];
	PDB->atom[ca].ELEMENT[2]='\0';
	if(++ca>nres){
	  fclose(data);
	  return 0;
	}
      }
    }/* <---- End of 'if(!strncmp(HED,"ATOM  ",6)){... */


    /* ---------- Check for lines to be included in the header ---------- */
    else if(!strncmp(HED,"HEADER",6) || !strncmp(HED,"TITLE ",6) || 
	    !strncmp(HED,"COMPND",6) || !strncmp(HED,"SEQRES",6) ||  
	    !strncmp(HED,"HELIX ",6) || !strncmp(HED,"SHEET ",6) || 
	    !strncmp(HED,"TURN  ",6)){
      sprintf(PDB->HEADER[hd],"%s",LINE);
      hd++;
    }

    /* ---------- Check for MODEL lines ----------- */
    else if(!strncmp(HED,"ENDMDL",6)) 
      mdl++;

    /* ---------- Check for end of file ---------- */
    else if(!strncmp(HED,"END",3) || feof(data)){
      fclose(data);
      return 0;
    }
  }
}


/* "read_sparsemtx" READS A SPARSE MATRIX FROM FILE */
void read_sparsemtx(char *file,double **MM,int row,int col)
{
  FILE *data;
  double x;
  int i,j;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_sparsemtx: unable to open %s\n\n",file);
    exit(1);}
  for(i=1;i<=row;i++)
    for(j=1;j<=col;j++)
      MM[i][j]=0.0;
  while(!feof(data)){
    fscanf(data,"%d%d%lf",&i,&j,&x);
    if(row==col)
      MM[i][j]=MM[j][i]=x;
    else
      MM[i][j]=x;
  }
  fclose(data);
}



/* "righthand2" MAKES SURE THAT THE EIGENVECTORS 
   FORM A RIGHT-HANDED COORDINATE SYSTEM */
void righthand2(double *VAL,double **VEC,int n)
{
  double A[3],B[3],C[3],CP[3],dot=0.0;
  int i;

  /* FIND THE CROSS PRODUCT OF THE FIRST TWO EIGENVECTORS */
  for(i=0;i<3;i++){
    A[i]=VEC[i+1][1];
    B[i]=VEC[i+1][2];
    C[i]=VEC[i+1][3];}
  cross(A,B,CP);

  /* PROJECT IT ON THE THIRD EIGENVECTOR */
  for(i=0; i<3; i++)
    dot+=C[i]*CP[i];
  if(dot<0.0)
    for(i=1;i<=3;i++)
      VEC[i][3]=-VEC[i][3];
}



/* "unit_imatrix" ALLOCATES MEMORY FOR A UNIT MATRIX */
int **unit_imatrix(long lo,long hi)
{
  static int **M;
  int i,j;

  M=imatrix(lo,hi,lo,hi);
  for(i=lo;i<=hi;i++){
    M[i][i]=1;
    for(j=i+1;j<=hi;j++)
      M[i][j]=M[j][i]=0;
  }
  return M;
}



/* "zero_d3tensor" ALLOCATES MEMORY FOR A DOUBLE 
   3-TENSOR AND INITIALIZES IT TO ZERO */
double ***zero_d3tensor(long nrl,long nrh,long ncl,long nch,long ndl,long ndh)
{
  static double ***T;
  int i,j,k;

  T=d3tensor(nrl,nrh,ncl,nch,ndl,ndh);
  for(i=nrl;i<=nrh;i++)
    for(j=ncl;j<=nch;j++)
      for(k=ndl;k<=ndh;k++)
	T[i][j][k]=0.0;
  return T;
}



/* "zero_dmatrix" ALLOCATES MEMORY FOR A 
   DOUBLE MATRIX AND INITIALIZES IT TO ZERO */
double **zero_dmatrix(long nrl,long nrh,long ncl,long nch)
{
  static double **M;
  int i,j;

  M=dmatrix(nrl,nrh,ncl,nch);
  for(i=nrl;i<=nrh;i++)
    for(j=ncl;j<=nch;j++)
      M[i][j]=0.0;
  return M;
}


/* ------------------ Numerical Recipes Routines ------------------ */
double dpythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void dsvdcmp(double **a, int m, int n, double w[], double **v)
{
	double dpythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;
	static int maxits=100;

	rv1=dvector(1,n);
	g=scale=anorm=0.0;
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) {
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) {
		for (its=1;its<=maxits;its++) {
			flag=1;
			for (l=k;l>=1;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == maxits) nrerror("no convergence in many dsvdcmp iterations");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);
}

void deigsrt(double d[], double **v, int n)
{
	int k,j,i;
	double p;

	for (i=1;i<n;i++) {
		p=d[k=i];
		for (j=i+1;j<=n;j++)
			if (d[j] >= p) p=d[k=j];
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			for (j=1;j<=n;j++) {
				p=v[j][i];
				v[j][i]=v[j][k];
				v[j][k]=p;
			}
		}
	}
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}



double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t=(double ***) malloc((size_t)((nrow+1)*sizeof(double**)));
  if (!t) nrerror("allocation failure 1 in d3tensor()");
  t += 1;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(double **) malloc((size_t)((nrow*ncol+1)*sizeof(double*)));
  if (!t[nrl]) nrerror("allocation failure 2 in d3tensor()");
  t[nrl] += 1;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+1)*sizeof(double)));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in d3tensor()");
  t[nrl][ncl] += 1;
  t[nrl][ncl] -= ndl;

  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}



void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh)
/* free a double d3tensor allocated by d3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-1));
	free((FREE_ARG) (t[nrl]+ncl-1));
	free((FREE_ARG) (t+nrl-1));
}


