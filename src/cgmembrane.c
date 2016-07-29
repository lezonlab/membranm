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

  This code creates a PDB file of a coarse-grained membrane that fits around a 
  protein.             
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "membranmutil.h"
#define DEFHILO 13.0 /* Default lower/upper bound of membrane */
#define DEFSPN 80.0  /* Default radius of membrane segment */
#define DEFRAD 2.5   /* Default radius of membrane centroid */

void read_command_line(int,char*[],double*,double*,double*,double*,int*,int*);
int protclash(double [],PDB_File *,int,double);
void assign_lpvs(double [][3],int);

int main(int argc,char *argv[])
{
  PDB_File PDB;
  double mass;
  double mhi,mlo,rho,rhsq,rr,dd,ds;
  double pxhi,pxlo,pyhi,pylo,pzhi,pzlo;
  double A[3][3],X[3];
  int lat;
  int atm,who,nh,nm,nres,imax,jmax,kmax,i,j,k,p;

  read_command_line(argc,argv,&mlo,&mhi,&rho,&rr,&who,&lat);


  /* Read PDB file */
  if(who!=0){
    fprintf(stderr,"Reading pdb file %s\n",argv[who]);
    pdb_hmr(argv[who],&nh,&nm,&nres);
    pdb_init(&PDB,nh,nres);
    read_pdb3(argv[who],&PDB,nh,nres);

    /* Find the extent of the protein */
    pxhi=pyhi=pzhi=-10000.0;
    pxlo=pylo=pzlo=10000.0;
    for(i=1;i<=nres;i++){
      if(PDB.atom[i].X[0]<pxlo) pxlo=PDB.atom[i].X[0];
      if(PDB.atom[i].X[0]>pxhi) pxhi=PDB.atom[i].X[0];
      if(PDB.atom[i].X[1]<pylo) pylo=PDB.atom[i].X[1];
      if(PDB.atom[i].X[1]>pyhi) pyhi=PDB.atom[i].X[1];
      if(PDB.atom[i].X[2]<pzlo) pzlo=PDB.atom[i].X[2];
      if(PDB.atom[i].X[2]>pzhi) pzhi=PDB.atom[i].X[2];
    }
  }


  /* Calculate the maximum number of points in the volume
  mxpts=(int)(0.75*(mhi-mlo)*rho*rho/rr/rr/rr); */


  /* Assign lattice primitive vectors */
  assign_lpvs(A,lat);

  /* Number of points : imax, jmax and kmax are the expected limits 
     along the first, second and third primitive vectors.  If d is a 
     vector that points from the center of the membrane fragment to the
     point that is maximally distant along the first primitive vector, then
     the projection of d onto A[0] is the required length along A[0].
     The associated number of lattice points is this length divided by the
     centroid radius.  The maximum projection in the x-y plane is rho and 
     the maximum projection in the z-direction is half the membrane thickness. 
  */
  imax = (rho + A[0][2]*(mhi-mlo)/2.0)/rr;
  jmax = (rho + A[1][2]*(mhi-mlo)/2.0)/rr;
  kmax = (rho + A[2][2]*(mhi-mlo)/2.0)/rr;
  
  
  fprintf(stderr,"imax\t%d\njmax\t%d\nkmax\t%d\n",imax,jmax,kmax);


  /* Print PDB header */
  printf("REMARK   7\n");
  printf("REMARK   7 Generated with command:\n");
  printf("REMARK   7 ");
  for(i=0;i<argc;i++) printf("%s ",argv[i]);
  printf("\n");


  /* Identify all candidate lattice sites */
  atm=0;
  rhsq=rho*rho;
  ds=2.0*rr;
  for(i=-imax;i<=imax;i++){
    for(j=-imax;j<=imax;j++){
      for(k=-kmax;k<=kmax;k++){
	for(p=0;p<3;p++){
	  X[p]=ds*((double)i*A[0][p]+(double)j*A[1][p]+(double)k*A[2][p]);
	}

	/* Keep lattice site if it's not too far out */
	dd=0.0;
	for(p=0;p<2;p++) dd+=X[p]*X[p];
	if(dd<=rhsq && X[2]>mlo && X[2]<mhi){

	  /* Keep all points when no PDB file is specified */
	  /* Condition for x-y plane: fabs(X[2])<1.0e-3 */
	  /* Condition for fragment edge: dd>rhsq+4.0*rr*rr-4.0*rr*rho */
	  if(who==0 && dd<rhsq+4.0*rr*rr-4.0*rr*rho){
	    atm++;
	    printf("ATOM% 7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n",
		   atm,atm,X[0],X[1],X[2]);
	  }

	  /* Look for clashes if site is close to protein */
	  else if(X[0]>pxlo && X[0]<pxhi && X[1]>pylo && 
		  X[1]<pyhi && X[2]>pzlo && X[2]<pzhi){
	    if(protclash(X,&PDB,nres,5.0)==0){/* && dd>389.0){ */
	      atm++;
	      printf("ATOM% 7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n",
		     atm,atm,X[0],X[1],X[2]);
	    }
	  }

	  else if(who!=0){
	    atm++;
	    printf("ATOM% 7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n",
		   atm,atm,X[0],X[1],X[2]);
	  }

	}
      }
    }
  }
  printf("END");
  for(i=4;i<=80;i++) printf(" ");
  printf("\n");


  /* Centroid mass = (membrane density)*(membrane volume)/atm 
     The prefactor of 0.005473 converts g/cm^3 to m_0/\AA^3*/
  mass=0.005473*0.9*M_PI*rhsq*(mhi-mlo)/(double)atm;

  fprintf(stderr,"%d centroids\n",atm);
  fprintf(stderr,"Centroid mass = %f*residue mass\n",mass);
  fprintf(stderr,"              = %f Da, assuming residue mass of 110 Da\n",
	  110.0*mass);
  return 0;
}


/* "read_command_line" CHECKS FOR THE INPUT ARGUMENTS */
void read_command_line(int argc,char *argv[],double *mlo,double *mhi,
		       double *rho,double *rad,int *who,int *lat)
{
  double x;
  int i,ok=0,nopt=0;

  /* ASSIGN DEFAULT VALUES */
  *mlo=-DEFHILO;
  *mhi=DEFHILO;
  *rho=DEFSPN;
  *rad=DEFRAD;
  *who=0;
  *lat=0;
  

  /* ALLOWED FLAGS: 
     -b    Specify membrane boundaries
     -r    Specify radius of cylindrical membrane fragment
     -s    Specify size (radius) of membrane particle
     -f    Force a membrane to print with no input PDB file
     -sh,-sc,-fcc  Lattice types: Simple Hexagonal, Simple Cubic, FCC
     :TODO: -c Specify center of membrane fragment 
  */

  if(argc>1 && argc<=10){
    ok=1;
    i=1;
    while(i<argc){

      /* Membrane boundaries */
      if(strcmp(argv[i],"-b")==0){
	sscanf(argv[++i],"%lf",mlo);
	sscanf(argv[++i],"%lf",mhi);
	if(*mlo > *mhi){
	  x=*mhi;
	  *mhi=*mlo;
	  *mlo=x;
	}
	nopt++;
      }

      /* Centroid radius */
      else if(strcmp(argv[i],"-s")==0){
	sscanf(argv[++i],"%lf",rad);
	nopt++;
      }

      /* Membrane fragement radius */
      else if(strcmp(argv[i],"-r")==0){
	sscanf(argv[++i],"%lf",rho);
	nopt++;
      }

      /* Force output, even without protein */
      else if(strcmp(argv[i],"-f")==0){
	nopt++;
      }

      /* Alternative lattice */
      else if(strcmp(argv[i],"-fcc")==0) nopt++;
      else if(strcmp(argv[i],"-sh")==0){
	if(*lat==0) *lat=1;
	else{
	  fprintf(stderr,"\nMultiple lattices specified. Aborting run\n\n");
	  exit(1);}
	nopt++;
      }
      else if(strcmp(argv[i],"-sc")==0){
	if(*lat==0) *lat=2;
	else{
	  fprintf(stderr,"\nMultiple lattices specified. Aborting run\n\n");
	  exit(1);}
	nopt++;
      }

      /* Assumed PDB file */
      else{
	*who=i;
	nopt++;
      }
      i++;
    }
    if(nopt>6) ok=0;
  }

  if(ok==0){
    fprintf(stderr,"\nUsage:\n%s [file.pdb] [OPTIONS]\n\n",argv[0]);
    fprintf(stderr,"file.pdb is an optional PDB file\n\n");
    fprintf(stderr,"OPTIONS:\n");
    fprintf(stderr,"\t-s RAD\tSpecify centroid radius in Angstroms (default %f)\n",DEFRAD);
    fprintf(stderr,"\t-r RHO\tSpecify cylindrical membrane fragment radius in Angstroms (default %f)\n",DEFSPN);
    fprintf(stderr,"\t-b MLO MHI\tSpecify membrane boundaries (default %f,%f)\n",-DEFHILO,DEFHILO);
    fprintf(stderr,"\t-f Force membrane creation without input PDB file\n");
    fprintf(stderr,"\t-sh,-sc SH or SC lattice (default FCC)\n\n");
    fprintf(stderr,"Output:\nPrints PDB file of membrane fragment\n\n");
    exit(1);
  }


  fprintf(stderr,"Membrane limits: %f < z < %f\n",*mlo,*mhi);
  fprintf(stderr,"Fragment radius: %f\n",*rho);
  fprintf(stderr,"Centroid radius: %f\n",*rad);
  if(*lat==0) fprintf(stderr,"FCC lattice\n");
  else if(*lat==1) fprintf(stderr,"SH lattice\n");
  else if(*lat==2) fprintf(stderr,"SC lattice\n");

  return;
}


/* "protclash" Determines whether a point and a 
   protein clash to within the specified distance */
int protclash(double X[],PDB_File *PDB,int nres,double rad)
{
  int i,k;
  double dd,df,rsq;

  rsq=rad*rad;
  for(i=1;i<=nres;i++){
    dd=0.0;
    for(k=0;k<3;k++){
      df=PDB->atom[i].X[k]-X[k];
      dd+=df*df;
    }
    if(dd<=rsq) return 1;
  }
  return 0;
} 


/* "assign_lpvs" assigns lattice primitive vectors. */
void assign_lpvs(double A[][3],int lat)
{
  /* Lattice primitive vectors:
     FCC: 1/sqrt(2)*(0,1,1); 1/sqrt(2)*(1,0,1); 1/sqrt(2)*(1,1,0)
     SC: (1,0,0); (0,1,0); (0,0,1)
     SH/HCP: 0.5*(1,-sqrt(3),0); 0.5*(1,sqrt(3),0); (0,0,1)
  */

  /* FCC: 1/sqrt(2)*(0,1,1); 1/sqrt(2)*(1,0,1); 1/sqrt(2)*(1,1,0) */
  if(lat==0){
    A[0][0]=0.0;
    A[0][1]=M_SQRT1_2;
    A[0][2]=M_SQRT1_2;
    A[1][0]=M_SQRT1_2;
    A[1][1]=0.0;
    A[1][2]=M_SQRT1_2;
    A[2][0]=M_SQRT1_2;
    A[2][1]=M_SQRT1_2;
    A[2][2]=0.0;
  }


  /* SH: 0.5*(1,-sqrt(3),0); 0.5*(1,sqrt(3),0); (0,0,1) */
  else if(lat==1){
    A[0][0]=0.5;
    A[0][1]=-sqrt(3.0)/2.0;
    A[0][2]=0.0;
    A[1][0]=0.5;
    A[1][1]=sqrt(3.0)/2.0;
    A[1][2]=0.0;
    A[2][0]=0.0;
    A[2][1]=0.0;
    A[2][2]=1.0;
  }

  /* SC: (1,0,0); (0,1,0); (0,0,1) */
  else if(lat==2){
    A[0][0]=1.0;
    A[0][1]=0.0;
    A[0][2]=0.0;
    A[1][0]=0.0;
    A[1][1]=1.0;
    A[1][2]=0.0;
    A[2][0]=0.0;
    A[2][1]=0.0;
    A[2][2]=1.0;
  }

  else{
    fprintf(stderr,"\nassign_lpvs: unknown lattice type\n\n");
    exit(1);}
}
