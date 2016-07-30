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

int main(int argc,char *argv[])
{
  PDB_File PDB;
  double mass;
  double mhi,mlo,rho,rhsq,centroid_radius,dd,ds;
  double pxhi,pxlo,pyhi,pylo,pzhi,pzlo;
  double X[3];
  double root3;
  int lat;
  int atm,who,nh,nm,nres,imax,jmax,kmax,i,j,k,p,q;

  read_command_line(argc,argv,&mlo,&mhi,&rho,&centroid_radius,&who,&lat);

  root3 = sqrt(3.0);

  /* Read PDB file */
  if(who!=0){
    fprintf(stderr,"Reading pdb file %s\n",argv[who]);
    pdb_hmr(argv[who],&nh,&nm,&nres);
    pdb_init(&PDB,nh,nres);
    read_pdb3(argv[who],&PDB,nh,nres);

    /* Find the extent of the protein */
    pxhi = pyhi = pzhi = -10000.0;
    pxlo = pylo = pzlo = 10000.0;
    for(i=1; i<=nres; i++){
      if(PDB.atom[i].X[0] < pxlo) pxlo = PDB.atom[i].X[0];
      if(PDB.atom[i].X[0] > pxhi) pxhi = PDB.atom[i].X[0];
      if(PDB.atom[i].X[1] < pylo) pylo = PDB.atom[i].X[1];
      if(PDB.atom[i].X[1] > pyhi) pyhi = PDB.atom[i].X[1];
      if(PDB.atom[i].X[2] < pzlo) pzlo = PDB.atom[i].X[2];
      if(PDB.atom[i].X[2] > pzhi) pzhi = PDB.atom[i].X[2];
    }
  }


  /* For FCC lattice:
     1. A line of points from (-imax,0,0) to (imax,0,0)
     2. The same thing as 1, but repeated from z=-kmax to z=kmax
     3. More stuff
  */

  /* imax, jmax, kmax are x- y- z- limits */
  imax = jmax = rho/centroid_radius/2;
  kmax = (mhi-mlo)/centroid_radius/4;
  
  fprintf(stderr,"imax\t%d\njmax\t%d\nkmax\t%d\n",imax,jmax,kmax);

  /* Print PDB header */
  printf("REMARK   7\n");
  printf("REMARK   7 Generated with command:\n");
  printf("REMARK   7 ");
  for(i=0; i<argc; i++) printf("%s ",argv[i]);
  printf("\n");

  /* Identify all candidate lattice sites */
  /* Construct a hexagonal lattice in the x-y plane (q=0), with one LPV along the x-axis, and one site
     at the origin: Place sites along the x-axis (p=0), and duplicate them for non-zero values of y. 
     Generate sites slighly offset (p=1), and duplicate these along y. Duplicate the entire plane for
     discrete z-values. Similarly, generate an offset plane (q=1) by constructing two rows (p=0,1) and 
     repeating them along y- and z-axes. */
  atm = 0;
  rhsq = rho*rho;
  ds = 2.0*centroid_radius;
  for(i=-imax; i<=imax; i++){
    //jmax = sqrt((double)(imax*imax - i*i));
    for(j=-jmax; j<=jmax; j++){
      for(q=0; q<=1; q++){
	for(p=0; p<=1; p++){
	  X[0] = q==0 ? (double)i*ds + (double)p*centroid_radius : (double)i*ds + (1-p)*centroid_radius;
	  X[1] = (double)j*root3*ds + (double)q*centroid_radius/root3 + (double)p*root3*centroid_radius;

	  dd = X[0]*X[0] + X[1]*X[1];
	  if(dd <= rhsq){
	    for(k=-kmax; k<=kmax; k++){
	      X[2] = (double)k*ds + (double)q*M_SQRT2/root3*centroid_radius;
	  
	      /* Keep all points when no PDB file is specified */
	      if(who==0){
		atm++;
		printf("ATOM% 7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n",atm,atm,X[0],X[1],X[2]);
	      }

	      /* Look for clashes if site is close to protein */
	      else if(X[0]>pxlo && X[0]<pxhi && X[1]>pylo && X[1]<pyhi && X[2]>pzlo && X[2]<pzhi){
		if(protclash(X,&PDB,nres,5.0)==0){
		  atm++;
		  printf("ATOM% 7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n",atm,atm,X[0],X[1],X[2]);
		}
	      }
	      else if(who!=0){
		atm++;
		printf("ATOM% 7d  Q1  NE1 Q%4d% 12.3f% 8.3f% 8.3f\n",atm,atm,X[0],X[1],X[2]);
	      }
	    }
	  }
	}
      }
    }
  }
  printf("END");
  for(i=4; i<=80; i++) printf(" ");
  printf("\n");


  /* Centroid mass = (membrane density)*(membrane volume)/atm 
     The prefactor of 0.005473 converts g/cm^3 to m_0/\AA^3*/
  mass = 0.005473*0.9*M_PI*rhsq*(mhi-mlo)/(double)atm;

  fprintf(stderr,"%d centroids\n",atm);
  fprintf(stderr,"Centroid mass = %f*residue mass\n",mass);
  fprintf(stderr,"              = %f Da, assuming residue mass of 110 Da\n",110.0*mass);
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
