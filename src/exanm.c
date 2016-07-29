/*
  Licensed under the MIT License (MIT)

  Copyright (c) 2008-2015 Timothy Lezon

  This file is part of the exANM software package

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

  This code calculates the effective hessian matrix of a protein that is 
  influenced by an explicit environment.                
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include "membranmutil.h"
#define DEFCUT 7.5 /* default cutoff */
#define DEFGAM 1.0 /* default gamma */
#define DELTA 1.5 /* Distance, in \AA, of zero-crossing from membrane surface */
#define EPS 0.2   /* Scale factor in center of membrane */


typedef struct {char **name;double **M;} Sprngmtx;

void read_command_line(int,char *[],int*,int*,int*,double*);
int num_clusters(int **,int);
void radius_contact_sysenv(int **,Centroid *,Centroid *,int,int,double *);
double *read_cutfile2(char *,Centroid *,Centroid *,int,int);
Sprngmtx *read_springfile_sysenv(char *,Centroid *,Centroid *,int,int,int *);
void spring_constants_sysenv(Centroid *,Centroid *,double **,int **,Sprngmtx *,int,int,int);
void mwhess_sysenv(double **,double **,double **,Centroid *,Centroid *,double **,int,int);
void read_contacts(char *,int **,int);
double scalefunc1(double,double);
void membounds(Centroid*,int,double*,double*);

/* Variables for lateral pressure profile */
double MHI,MLO,MSCL;

int main(int argc,char *argv[])
{
  Centroid *SYS,*ENV;
  Sprngmtx *Gamma;
  char *sysfile,*envfile,*massfile,*ctfile,*spgfile;
  double **GG;
  double *CUT,**HS,**HE,**HX,**PH,*P;
  double dd,x,y;
  int **CT,**PIX,nss,nen,nn,ntp,nse,r1,c1,r2,c2,i,j,k;
  int ptmd,rprm,pprm;


  /* Formalities */
  read_command_line(argc,argv,&pprm,&rprm,&ptmd,&MSCL);


  /* Read coordinate and mass data */
  sysfile=get_param(argv[pprm],"syscoords");
  envfile=get_param(argv[pprm],"envcoords");
  massfile=get_param(argv[pprm],"massfile");
  SYS=read_centroids1(sysfile,massfile,&nss);
  ENV=read_centroids1(envfile,massfile,&nen);
  nn=nss+nen;
  fprintf(stderr,"System read from %s: %d centroids\n",sysfile,nss);
  fprintf(stderr,"Environment read from %s: %d centroids\n",envfile,nen);
  fprintf(stderr,"Masses read from %s\n",massfile);
  CT=imatrix(1,nn,1,nn);


  /* Find extent of membrane */
  membounds(ENV,nen,&MHI,&MLO);


  /* Print masses, if called for */
  if(ptmd==2){
    for(i=1;i<=nss;i++)
      for(j=-2;j<=0;j++){
	k=3*i+j;
	printf("%8d%8d% 25.15e\n",k,k,SYS[i].mass);
      }
    for(i=1;i<=nen;i++)
      for(j=-2;j<=0;j++){
	k=3*nss+3*i+j;
	printf("%8d%8d% 25.15e\n",k,k,ENV[i].mass);
      }
    return 0;
  }


  /* ---------------- Assign contacts... ----------------- */
  /* From a contact file... */
  if((ctfile=get_param(argv[pprm],"contactfile"))!=NULL){
    read_contacts(ctfile,CT,nn);
    fprintf(stderr,"Contacts read from %s\n",ctfile);
  }
  /* ...or from a cutoff file... */
  else if((ctfile=get_param(argv[pprm],"cutfile"))!=NULL){
    fprintf(stderr,"Cutoff values read from %s\n",ctfile);
    CUT=read_cutfile2(ctfile,SYS,ENV,nss,nen);
    radius_contact_sysenv(CT,SYS,ENV,nss,nen,CUT);
  }
  /* ...or from default values */
  else{
    CUT=dvector(1,nn);
    for(i=1;i<=nn;i++) CUT[i]=DEFCUT;
    fprintf(stderr,"All cutoff values set to %.3f\n",DEFCUT);
    radius_contact_sysenv(CT,SYS,ENV,nss,nen,CUT);
  }
  fprintf(stderr,"%d clusters\n",num_clusters(CT,nn));

  if(ptmd==1){
    fprintf(stderr,"Printing contacts\n");
    for(i=1;i<=nn;i++)
      for(j=i+1;j<=nn;j++)
	if(CT[i][j]!=0)
	  printf("%d\t%d\n",i,j);
    return 0;
  }

    



  /* ------- Construct the matrix of force constants -------*/
  GG=dmatrix(1,nn,1,nn);
  
  /* Read force constants from file... */
  if((spgfile=get_param(argv[pprm],"springfile"))!=NULL){
    fprintf(stderr,"Reading spring constants from %s\n",spgfile);
    Gamma=read_springfile_sysenv(spgfile,SYS,ENV,nss,nen,&ntp);
    spring_constants_sysenv(SYS,ENV,GG,CT,Gamma,nss,nen,ntp);
  }
  /* ...or else assign the default value to all springs */
  else
    for(i=1;i<=nn;i++)
      for(j=i;j<=nn;j++)
	if(CT[i][j]!=0)
	  GG[i][j]=GG[j][i]=DEFGAM;



  /* Construct the mass-weighted Hessian from 
     coordinates, masses, and potential matrix */
  fprintf(stderr,"Calculating Hessian...\n");
  HS=dmatrix(1,3*nss,1,3*nss);
  HE=dmatrix(1,3*nen,1,3*nen);
  HX=dmatrix(1,3*nss,1,3*nen);
  mwhess_sysenv(HS,HE,HX,SYS,ENV,GG,nss,nen);
      


  /* PRINT THE ENVIRONMENT-ENVIRONMENT SUB-HESSIAN IN SPARSE FORMAT */
  if(ptmd==0 && rprm==-1){
    fprintf(stderr,"\nPrinting env-env sub-hessian...\n\n");
    for(i=1;i<=3*nen;i++)
      for(j=i;j<=3*nen;j++)
	if(fabs(HE[i][j])>1.0e-10)
	  printf("%8d%8d% 25.15e\n",i,j,HE[i][j]);
    return 0;
  }


  /* PRINT THE FULL HESSIAN IN SPARSE FORMAT */
  if(ptmd==3){
    for(i=1;i<=3*nss;i++){
      for(j=i;j<=3*nss;j++)
	if(fabs(HS[i][j])>1.0e-10)
	  printf("%8d%8d% 20.10e\n",i,j,HS[i][j]);
      for(j=1;j<=3*nen;j++)
	if(fabs(HX[i][j])>1.0e-10)
	  printf("%8d%8d% 20.10e\n",i,j+3*nss,HX[i][j]);
    }
    for(i=1;i<=3*nen;i++)
      for(j=i;j<=3*nen;j++)
	if(fabs(HE[i][j])>1.0e-10)
	  printf("%8d%8d% 20.10e\n",i+3*nss,j+3*nss,HE[i][j]);
    return 0;
  }


  /* READ INVERSE OF ENVIRONMENTAL HESSIAN, OR INVERT HE */
  free_imatrix(CT,1,nn,1,nn);
  free_dmatrix(GG,1,nn,1,nn);
  if(rprm!=-1){
    fprintf(stderr,"Reading matrix from %s...\n",argv[rprm]);
    read_sparsemtx(argv[rprm],HE,3*nen,3*nen);
  }
  else{
    fprintf(stderr,"\nWell...How did I get here?\n\n");
    exit(1);}


  /* ---------------- CALCULATE AND PRINT THE PSEUDOHESSIAN ---------------- */

  /* COUNT THE NUMBER OF NON-ZERO TERMS IN THE PROJECTION MATRIX */
  nse=0;
  for(i=1;i<=3*nss;i++)
    for(j=1;j<=3*nen;j++)
      if(fabs(HX[i][j])>1.0e-9) nse++;
  fprintf(stderr,"%d non-zero projection elements\n",nse);
  P=dvector(1,nse);
  PIX=imatrix(1,nse,1,2);
  k=1;
  for(i=1;i<=3*nss;i++)
    for(j=1;j<=3*nen;j++)
      if(fabs(HX[i][j])>1.0e-9){
	PIX[k][1]=i;
	PIX[k][2]=j;
	P[k]=HX[i][j];
	k++;
      }
  free_dmatrix(HX,1,3*nss,1,3*nen);
  PH=dmatrix(1,3*nss,1,3*nss);
  for(i=1;i<=3*nss;i++)
    for(j=i;j<=3*nss;j++)
      PH[i][j]=PH[j][i]=0.0;
  for(i=1;i<=nse;i++){
    r1=PIX[i][1];
    c1=PIX[i][2];
    x=P[i];
    for(j=i;j<=nse;j++){
      r2=PIX[j][1];
      c2=PIX[j][2];
      y=HE[c1][c2]*P[j]*x;
      PH[r1][r2]+=y;
      if(r1==r2 && c1!=c2)
	PH[r1][r2]+=y;
    }
  }

  for(i=1;i<=3*nss;i++)
    for(j=i;j<=3*nss;j++){
      dd=HS[i][j]-PH[i][j];
      if(fabs(dd)>1.0e-10)
	printf("%8d%8d% 25.15e\n",i,j,dd);
    }
  return 0;
}


/* "read_command_line" reads the command line arguments */
void read_command_line(int argc,char *argv[],int *pprm,int *rprm,int *prmd,double *mscl)
{
  int i,ok=0,nopt=0;
  
  *pprm=*rprm=-1;
  *prmd=0;
  *mscl=1.0;

  /* ALLOWED FLAGS: 
     -c Prints the connectivity matrix
     -m prints the mass matrix
     -h Prints the full MW Hessian of the sys-env complex
     -r FILE uses FILE as the inverse MW env-env Hessian to restart
     -s Scale value used in membrane lateral pressure profile
  */

  if(argc>1 && argc<=7){
    ok=1;
    i=1;
    while(i<argc){

      /* Print connectivity matrix (print mode 1) */
      if(strcmp(argv[i],"-c")==0){
	if(*prmd==0) *prmd=1;
	else{
	  fprintf(stderr,
		  "\nMultiple print modes specified.  Aborting run.\n\n");
	  exit(1);}
	nopt++;
      }

      /* Print mass matrix (print mode 2) */
      else if(strcmp(argv[i],"-m")==0){
	if(*prmd==0) *prmd=2;
	else{
	  fprintf(stderr,
		  "\nMultiple print modes specified.  Aborting run.\n\n");
	  exit(1);}
	nopt++;
      }

      /* Print full Hessian (print mode 3) */
      else if(strcmp(argv[i],"-h")==0){
	if(*prmd==0) *prmd=3;
	else{
	  fprintf(stderr,
		  "\nMultiple print modes specified.  Aborting run.\n\n");
	  exit(1);}
	nopt++;
      }

      /* Restart from inverted env-env Hessian */
      else if(strcmp(argv[i],"-r")==0){
	*rprm=++i;
	nopt++;
      }

      /* Scale factor for membrane lateral pressure profile */
      else if(strcmp(argv[i],"-s")==0){
	sscanf(argv[++i],"%lf",mscl);
	nopt++;
      }

      /* Assumed parameter file */
      else{
	if(*pprm==-1) *pprm=i;
	else{
	  fprintf(stderr,
		  "\nMultiple parameter files specified.  Aborting run.\n\n");
	  exit(1);}
	nopt++;
      }
      i++;
    }
    if(nopt>3) ok=0;
  }

  if(ok==0){
    fprintf(stderr,"\nUsage:\n%s file.param [OPTIONS]\n\n",argv[0]);
    fprintf(stderr,"file.param is file containing at least 'syscoords', 'envcoords' and 'massfile' parameters\n\n");
    fprintf(stderr,"OPTIONS:\n");
    fprintf(stderr,"\t-c\tPrint contact matrix and exit\n");
    fprintf(stderr,"\t-m\tPrint mass matrix and exit\n");
    fprintf(stderr,"\t-h\tPrint full Hessian matrix and exit\n");
    fprintf(stderr,"\t-s MSCL\tSpecify membrane lateral pressure profile scaling parameter\n");
    fprintf(stderr,"\t-r FILE\tCalculate sys-sys pseudohessian using FILE as inverse of env-env sub-hessian\n\n");
    fprintf(stderr,"Output:\nPrints the env-env subhessian and exits\n\n");
    exit(1);
  }
  return;
}


/* "read_cutfile2" READS THE FILE OF HALF-CUTOFF DISTANCES 
   FOR ALL CENTROIDS IN BOTH THE SYSTEM AND THE ENVIRONMENT */
double *read_cutfile2(char *file,Centroid *SYS,Centroid *ENV,int nss,int nen)
{
  FILE *data;
  char NAME[200];
  double *X,x,defcut=DEFCUT;
  int nn=nss+nen,i;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_cutfile: unable to open %s\n\n",file);
    exit(1);}
  X=dvector(1,nn);
  while(!feof(data)){
    fscanf(data,"%s%lf",NAME,&x);
    if(!strcmp(NAME,"DEFCUT")){
      defcut=x;
      break;
    }
  }
  rewind(data);
  for(i=1;i<=nn;i++) X[i]=defcut;
  fprintf(stderr,"Default cutoff = %.3f\nOther values:\n",defcut);
  while(!feof(data)){
    fscanf(data,"%s%lf",NAME,&x);
    if(!feof(data)){
      fprintf(stderr,"%s\t%.3f\n",NAME,x);
      for(i=1;i<=nss;i++)
	if(!strcmp(SYS[i].name,NAME))
	  X[i]=x;
      for(i=1;i<=nen;i++)
	if(!strcmp(ENV[i].name,NAME))
	  X[nss+i]=x;
    }
  }
  fclose(data);
  return X;
}


/* "read_springfile_sysenv" READS THE SPRINGFILE FOR 
   BOTH THE 'SYSTEM' AND THE 'ENVIRONMENT' */
Sprngmtx *read_springfile_sysenv(char *file,Centroid *SYS,Centroid *ENV,
			   int nss,int nen,int *ntp)
{
  FILE *data;
  Sprngmtx *foo;
  char **LIST,nup1[NAME_LNG],nup2[NAME_LNG];
  double x;
  int nn=nss+nen,ok,i,j;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_springfile_sysenv: unable to open %s\n\n",file);
    exit(1);}

  /* GET THE LIST OF UNIQUE CENTROID NAMES */
  LIST=cmatrix(1,nn,0,NAME_LNG);
  (*ntp)=0;
  for(i=1;i<=nss;i++){
    ok=1;
    for(j=1;j<=(*ntp);j++){
      if(!strcmp(SYS[i].name,LIST[j])){
	ok=0;
	break;
      }
    }
    if(ok==1)
      strcpy(LIST[++(*ntp)],SYS[i].name);
  }
  for(i=1;i<=nen;i++){
    ok=1;
    for(j=1;j<=(*ntp);j++){
      if(!strcmp(ENV[i].name,LIST[j])){
	ok=0;
	break;
      }
    }
    if(ok==1)
      strcpy(LIST[++(*ntp)],ENV[i].name);
  }

  foo=(Sprngmtx *)malloc((size_t) sizeof(Sprngmtx *));
  foo->name=cmatrix(1,(*ntp),0,NAME_LNG);
  foo->M=dmatrix(1,(*ntp),1,(*ntp));

  for(i=1;i<=(*ntp);i++){
    strcpy(foo->name[i],LIST[i]);
    for(j=i;j<=(*ntp);j++)
      foo->M[i][j]=foo->M[j][i]=DEFGAM;
  }
  free_cmatrix(LIST,1,nn,0,NAME_LNG);

  /* READ THE FILE */
  while(!feof(data)){
    fscanf(data,"%s%s%lf",nup1,nup2,&x);
    for(i=1;i<=(*ntp);i++)
      if(!strcmp(nup1,foo->name[i])){
	for(j=1;j<=(*ntp);j++)
	  if(!strcmp(nup2,foo->name[j])){
	    foo->M[i][j]=foo->M[j][i]=x;
	    break;
	  }
	break;
      }
  }
  fclose(data);
  return foo;
}


/* "num_clusters" RETURNS THE NUMBER OF CLUSTERS IN THE 
   NETWORK DESCRIBED BY THE CONNECTIVITY MATRIX CT */
int num_clusters(int **CT,int n)
{
  int *ID,dun,nc,fp,i,j,ii;

  ID=ivector(1,n);
  for(i=1;i<=n;i++) ID[i]=i;
  nc=0;
  for(ii=1;ii<=n;ii++){
    fp=1;
    dun=0;
    while(dun==0){
      dun=1;
      for(i=1;i<=n;i++){
	if(ID[i]==ii){

	  if(fp==1){
	    fp=0;
	    nc++;
	  }

	  for(j=1;j<=n;j++){
	    if(CT[i][j]==1 && ID[j]!=ii){
	      dun=0;
	      ID[j]=ii;
	    }
	  }
	}
      }
    }
  }
  free_ivector(ID,1,n);
  return nc;
}


/* "radius_contact_sysenv" FINDS THE ARRAY OF CONTACTS 
   BETWEEN 'SYSTEM' AND 'ENVIRONMENT' CENTROIDS BASED 
   ON THE PROVIDED HALF-RADII IN THE ARRAY 'CUT'. */
void radius_contact_sysenv(int **CT,Centroid *SYS,Centroid *ENV,
			   int nss,int nen,double *CUT)
{
  double df,ctsq,dsq;
  int i,j,k,qi,qj;

  for(i=1;i<=nss;i++){
    CT[i][i]=1;
    for(j=i+1;j<=nss;j++){
      dsq=0.0;
      for(k=1;k<=3;k++){
	df = SYS[i].X[k] - SYS[j].X[k];
	dsq+=(df*df);
      }
      ctsq=CUT[i]+CUT[j];
      ctsq*=ctsq;
      if(dsq<ctsq) CT[i][j]=CT[j][i]=1;
      else CT[i][j]=CT[j][i]=0;
    }
    for(j=1;j<=nen;j++){
      dsq=0.0;
      for(k=1;k<=3;k++){
	df = SYS[i].X[k] - ENV[j].X[k];
	dsq+=(df*df);
      }
      ctsq=CUT[i]+CUT[nss+j];
      ctsq*=ctsq;
      if(dsq<ctsq) CT[i][nss+j]=CT[nss+j][i]=1;
      else CT[i][nss+j]=CT[nss+j][i]=0;
    }
  }
  for(i=1;i<=nen;i++){
    qi=nss+i;
    CT[qi][qi]=1;
    for(j=i+1;j<=nen;j++){
      qj=nss+j;
      dsq=0.0;
      for(k=1;k<=3;k++){
	df = ENV[i].X[k] - ENV[j].X[k];
	dsq+=(df*df);
      }
      ctsq=CUT[qi]+CUT[qj];
      ctsq*=ctsq;
      if(dsq<ctsq) CT[qi][qj]=CT[qj][qi]=1;
      else CT[qi][qj]=CT[qj][qi]=0;
    }
  }
}


/* "spring_constants_sysenv" ASSIGNS SPRING CONSTANTS 
   BETWEEN CENTROIDS IN THE 'SYSTEM' AND THE 'ENVIRONMENT'. */
void spring_constants_sysenv(Centroid *SYS,Centroid *ENV,double **GG,int **CT,
			     Sprngmtx *Gamma,int nss,int nen,int ntp)
{
  int i,j,k,p,qi,qj;

  for(i=1;i<=nss;i++)
    for(k=1;k<=ntp;k++)
      if(!strcmp(Gamma->name[k],SYS[i].name)){
	for(j=i+1;j<=nss;j++)
	  if(CT[i][j]!=0)
	    for(p=1;p<=ntp;p++)
	      if(!strcmp(Gamma->name[p],SYS[j].name)){
		GG[i][j]=GG[j][i]=Gamma->M[k][p];
		break;
	      }
	for(j=1;j<=nen;j++){
	  qj=nss+j;
	  if(CT[i][qj]!=0)
	    for(p=1;p<=ntp;p++)
	      if(!strcmp(Gamma->name[p],ENV[j].name)){
		GG[i][qj]=GG[qj][i]=Gamma->M[k][p];
		break;
	      }
	}
	break;
      }
  for(i=1;i<=nen;i++){
    qi=nss+i;
    for(k=1;k<=ntp;k++)
      if(!strcmp(Gamma->name[k],ENV[i].name)){
	for(j=i+1;j<=nen;j++){
	  qj=nss+j;
	  if(CT[qi][qj]!=0)
	    for(p=1;p<=ntp;p++)
	      if(!strcmp(Gamma->name[p],ENV[j].name)){
		GG[qi][qj]=GG[qj][qi]=Gamma->M[k][p];
		break;
	      }
	}
	break;
      }
  }
}



/* "mwhess_sysenv" FINDS THE MASS-WEIGHTED HESSIANS, 'HS', 'HE', AND 'HX' 
   BETWEEN THE 'SYSTEM' AND 'ENVIRONMENT' GIVEN THE ARRAY 'GG' OF SPRING 
   CONSTANTS BETWEEN THE 'nss' SYS CENTROIDS AND THE 'nen' ENV CENTROIDS. */
void mwhess_sysenv(double **HS,double **HE,double **HX,Centroid *SYS,
		   Centroid *ENV,double **GG,int nss,int nen)
{
  double DX[4],*SQTS,*SQTE,df,dsq,scl;
  int i,j,k,ii,jj,qii,qjj;


  /* CLEAR DIAGONAL SUPER-ELEMENTS */
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++){
      for(ii=1;ii<=nss;ii++)
	HS[3*(ii-1)+i][3*(ii-1)+j]=0.0;
      for(ii=1;ii<=nen;ii++)
	HE[3*(ii-1)+i][3*(ii-1)+j]=0.0;
    }

  /* HS AND HX */
  for(ii=1;ii<=nss;ii++){

    /* System-system interactions */
    for(jj=ii+1;jj<=nss;jj++){
      dsq=0.0;
      for(k=1;k<=3;k++){
	DX[k] = SYS[ii].X[k] - SYS[jj].X[k];
	dsq+=(DX[k]*DX[k]);
      }
      for(i=1;i<=3;i++){
	for(j=i;j<=3;j++){
	  df=GG[ii][jj]*DX[i]*DX[j]/dsq;

	  /* OFF-DIAGONAL SUPER-ELEMENTS */
	  HS[3*(ii-1)+i][3*(jj-1)+j]=HS[3*(ii-1)+j][3*(jj-1)+i]=-df;
	  HS[3*(jj-1)+j][3*(ii-1)+i]=HS[3*(jj-1)+i][3*(ii-1)+j]=-df;

	  /* DIAGONAL SUPER-ELEMENTS */
	  HS[3*(ii-1)+i][3*(ii-1)+j]+=df;
	  HS[3*(jj-1)+i][3*(jj-1)+j]+=df;
	  if(i!=j){
	    HS[3*(ii-1)+j][3*(ii-1)+i]+=df;
	    HS[3*(jj-1)+j][3*(jj-1)+i]+=df;
	  }
	}
      }
    }

    /* System-environment interactions */
    for(jj=1;jj<=nen;jj++){

      /* Scaling factor */
      scl=scalefunc1(SYS[ii].X[3],ENV[jj].X[3]);

      qjj=jj+nss;
      dsq=0.0;
      for(k=1;k<=3;k++){
	DX[k] = SYS[ii].X[k] - ENV[jj].X[k];
	dsq+=(DX[k]*DX[k]);
      }
      for(i=1;i<=3;i++){
	for(j=i;j<=3;j++){
	  df=GG[ii][qjj]*DX[i]*DX[j]/dsq;


	  /* ------------------------ Perturbation ------------------------ */
	  /* This block keeps the z-component of SYS-ENV interactions at 
	     unity and scales the lateral components by either GG or 
	     sqrt(GG).  If this is only a lateral component
	     of the submatrix, then nothing is done.  If this is the z-z
	     component, then the force constant is either zero or unity.  If 
	     it is a mixed z/lateral component, then the force constant is 
	     either zero or sqrt(GG[ii][qjj]). 
	  if(i==3 && j==3)
	    df = fabs(df)>1.0e-6 ? df/GG[ii][qjj] : 0.0;
	  else if(i==3 || j==3)
	    df = fabs(df)>1.0e-6 ? df/sqrt(GG[ii][qjj]) : 0.0;
	   */
	  /* ------------------- End of perturbation ----------------------*/



	  /* -------------------------- Scaling -------------------------- */
	  /* In addition to using different force constants for SYS and ENV
	     nodes, the z-component of the force constants involving ENV 
	     is scaled by some factor */
	  if(i!=3) df*=scl;
	  if(j!=3) df*=scl;
	  /* ------------------------ End Scaling ------------------------ */



	  /* OFF-DIAGONAL SUPER-ELEMENTS */
	  HX[3*(ii-1)+i][3*(jj-1)+j]=HX[3*(ii-1)+j][3*(jj-1)+i]=-df;

	  /* DIAGONAL SUPER-ELEMENTS */
	  HS[3*(ii-1)+i][3*(ii-1)+j]+=df;
	  HE[3*(jj-1)+i][3*(jj-1)+j]+=df;
	  if(i!=j){
	    HS[3*(ii-1)+j][3*(ii-1)+i]+=df;
	    HE[3*(jj-1)+j][3*(jj-1)+i]+=df;
	  }
	}
      }
    }
  }

  /* Environment-environment interactions */
  for(ii=1;ii<=nen;ii++){
    qii=ii+nss;
    for(jj=ii+1;jj<=nen;jj++){

      /* Scaling factor */
      scl=scalefunc1(ENV[ii].X[3],ENV[jj].X[3]);

      qjj=jj+nss;
      dsq=0.0;
      for(k=1;k<=3;k++){
	DX[k] = ENV[ii].X[k] - ENV[jj].X[k];
	dsq+=(DX[k]*DX[k]);
      }

      for(i=1;i<=3;i++){
	for(j=i;j<=3;j++){
	  df=GG[qii][qjj]*DX[i]*DX[j]/dsq;

	  /* -------------------------- Scaling -------------------------- */
	  /* In addition to using different force constants for SYS and ENV
	     nodes, the z-component of the force constants involving ENV 
	     is scaled by some factor */
	  if(i!=3) df*=scl;
	  if(j!=3) df*=scl;
	  /* ------------------------ End Scaling ------------------------ */

	  /* OFF-DIAGONAL SUPER-ELEMENTS */
	  HE[3*(ii-1)+i][3*(jj-1)+j]=HE[3*(ii-1)+j][3*(jj-1)+i]=-df;
	  HE[3*(jj-1)+j][3*(ii-1)+i]=HE[3*(jj-1)+i][3*(ii-1)+j]=-df;

	  /* DIAGONAL SUPER-ELEMENTS */
	  HE[3*(ii-1)+i][3*(ii-1)+j]+=df;
	  HE[3*(jj-1)+i][3*(jj-1)+j]+=df;
	  if(i!=j){
	    HE[3*(ii-1)+j][3*(ii-1)+i]+=df;
	    HE[3*(jj-1)+j][3*(jj-1)+i]+=df;
	  }
	}
      }
    }
  }  

  /* ------- MASS-WEIGHTING ------- */
  SQTS=dvector(1,nss);
  SQTE=dvector(1,nen);
  for(i=1;i<=nss;i++) SQTS[i]=sqrt(SYS[i].mass);
  for(i=1;i<=nen;i++) SQTE[i]=sqrt(ENV[i].mass);
  for(ii=1;ii<=nss;ii++)
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++){
	for(jj=1;jj<=nss;jj++)
	  HS[3*(ii-1)+i][3*(jj-1)+j]/=(SQTS[ii]*SQTS[jj]);
	for(jj=1;jj<=nen;jj++)
	  HX[3*(ii-1)+i][3*(jj-1)+j]/=(SQTS[ii]*SQTE[jj]);
      }
  for(ii=1;ii<=nen;ii++)
    for(jj=1;jj<=nen;jj++)
      for(i=1;i<=3;i++)
	for(j=1;j<=3;j++)
	  HE[3*(ii-1)+i][3*(jj-1)+j]/=(SQTE[ii]*SQTE[jj]);
  free_dvector(SQTS,1,nss);
  free_dvector(SQTE,1,nen);
  /**/
}


/* "read_contacts" READS A FILE OF CONTACTS BETWEEN CENTROIDS */
void read_contacts(char *file,int **CT,int nn)
{
  FILE *data;
  int i,j;

  if((data=fopen(file,"r"))==NULL){
    fprintf(stderr,"\nread_contacts: unable to open %s\n\n",file);
    exit(1);}
  for(i=1;i<=nn;i++){
    CT[i][i]=1;
    for(j=i+1;j<=nn;j++)
      CT[i][j]=CT[j][i]=0;
  }
  while(!feof(data)){
    fscanf(data,"%d%d",&i,&j);
    CT[i][j]=CT[j][i]=1;
  }
  fclose(data);
}


/* "scalefunc1" scales the z-component of the Hessian by a depth-dependent 
   factor approximated by two parabolae.  The maximum value that can be 
   returned is sqrt(MSCL), if both nodes are in the region of highest 
   pressure; two nodes at z0 will have a scale factor of EPS*MSCL;
   the scale factor within DELTA from the membrane boundary is less than 1. */
double scalefunc1(double x1,double x2)
{
  double scl=1.0,s1=0.0,s2=0.0,x,eta=0.05;
  static double c1,c2,f0,l0,s0,z0,z1,z2;
  static int first=1;

  /* Scaling function is s(z)=exp(f(z)), where f(z) is parabolic. */

  /* CALCULATE THE CONSTANTS ON FIRST CALL */
  if(first==1){
    first=0;
    z0=(MHI+MLO)/2.0;        /* midpoint of membrane */
    s0=log(sqrt(MSCL));      /* highest scaling value is exp(s0) */
    f0=sqrt(s0);
    l0=sqrt(-log(sqrt(EPS)));
    c1=(l0+f0)/(z0-MLO-DELTA);
    c1*=c1;
    z1=MLO+DELTA+f0*(z0-MLO-DELTA)/(l0+f0);
    c2=(l0+f0)/(MHI-DELTA-z0);
    c2*=c2;
    z2=MHI-DELTA-f0*(MHI-DELTA-z0)/(l0+f0);
  }

  /* ADJUST THE SCALE */
  if(x1>MLO && x1<MHI){
    s1 = s0 - c1*(x1-z1)*(x1-z1);
    s2 = s0 - c2*(x1-z2)*(x1-z2);
    x=0.5*(s1 + s2 + sqrt((s2-s1)*(s2-s1) + 4.0*eta*eta));
    scl*=exp(x);
  }

  if(x2>MLO && x2<MHI){
    s1 = s0 - c1*(x2-z1)*(x2-z1);
    s2 = s0 - c2*(x2-z2)*(x2-z2);
    x=0.5*(s1 + s2 + sqrt((s2-s1)*(s2-s1) + 4.0*eta*eta));
    scl*=exp(x);
  }
  return scl;
}


/* "membounds" finds the z-extrema of the specified structure */
void membounds(Centroid *ENV,int nen,double *max,double *min)
{
  double x;
  int i;

  *max=*min=ENV[1].X[3];

  for(i=2;i<=nen;i++){
    x=ENV[i].X[3];
    if(x>*max) *max=x;
    if(x<*min) *min=x;
  }
}
