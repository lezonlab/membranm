/*
  Licensed under the MIT License (MIT)

  Copyright (c) 2008-2015 Timothy Lezon

  This is the main file of the imANM software package

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

  RTB ANM for membrane proteins                      

*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include "membranmutil.h"

#define DEFGAM 1.0   /* Default force constant */
#define DEFCUT 11.0  /* Default cutoff distance in angstroms */
#define DEFMSCL 16.0  /* Default membrane spring scale factor */
#define DEFMLO -13.4 /* Default lower membrane boundary */
#define DEFMHI 13.4  /* Default upper membrane boundary */

void read_command_line(int,char *[],int *,double *,double *,double *,double *,int *,double *);
int calc_blessian_mem(PDB_File *,dSparse_Matrix *,int,int,int,double **,double);
void hess_superrow_mem(double **,int **,PDB_File *,int,int,double,double);
void bless_to_hess(double **,double **,dSparse_Matrix *,int,int);
double scalefunc0(double,double);
void backproject(double **,dSparse_Matrix *,int,int,int,int);
int dblock_projections2(dSparse_Matrix *,PDB_File *,int,int,int);

double MHI,MLO,CUT,MSCL;

int main(int argc,char *argv[])
{
  Rigid_Block *BLX;
  PDB_File PDB1;
  dSparse_Matrix PP,HH;
  char *pdbfile;
  double **HB;
  double gam;
  int all,prm,nh1,nm1,nres,nrg,nblx,bmx,elm,imx,i,j;


  /* Formalities */
  read_command_line(argc,argv,&prm,&MLO,&MHI,&CUT,&MSCL,&all,&gam);
  read_blockfile(argv[1],&pdbfile,&BLX,&nrg);


  /* Read PDB File */
  pdb_hmr(pdbfile,&nh1,&nm1,&nres);
  pdb_init(&PDB1,nh1,nres);
  read_pdb3(pdbfile,&PDB1,nh1,nres);


  /* Assign each residue to a block */
  if(all==0){
    nblx=assign_rigid_blocks(&PDB1,BLX,nres,nrg,&bmx);
  }
  else{
    nblx=nres;
    bmx=1;
    for(i=1;i<=nres;i++)
      PDB1.atom[i].model=i;
  }
  fprintf(stderr,"\n%s: %d block ranges specified\n",argv[1],nrg);
  fprintf(stderr,"%s: %d residues, %d asymmetric units\n",pdbfile,nres,nm1);
  fprintf(stderr,"%d total rigid blocks in structure\n",nblx);
  fprintf(stderr,"largest block contains %d residues\n",bmx);
  free(pdbfile);
  free(BLX);


  /* Find the projection matrix */
  HH.IDX=imatrix(1,18*bmx*nblx,1,2);
  HH.X=dvector(1,18*bmx*nblx);
  elm=dblock_projections2(&HH,&PDB1,nres,nblx,bmx);
  PP.IDX=imatrix(1,elm,1,2);
  PP.X=dvector(1,elm);
  for(i=1;i<=elm;i++){
    PP.IDX[i][1]=HH.IDX[i][1];
    PP.IDX[i][2]=HH.IDX[i][2];
    PP.X[i]=HH.X[i];
  }
  free_imatrix(HH.IDX,1,18*bmx*nblx,1,2);
  free_dvector(HH.X,1,18*bmx*nblx);
  dsort_PP2(&PP,elm,1);


  /* Print the projection matrix */
  if(prm==1){
    print_prj_ofst(&PP,elm);
    return 0;
  }
  /**/


  /* Calculate the block Hessian */
  HB=dmatrix(1,6*nblx,1,6*nblx);
  imx=calc_blessian_mem(&PDB1,&PP,nres,nblx,elm,HB,gam);




  /* Print the block Hessian */
  for(i=1;i<=imx;i++)
    for(j=i;j<=imx;j++)
      if(fabs(HB[i][j])>1.0e-10)
	printf("%8d%8d% 20.10e\n",i,j,HB[i][j]);
  return 0;
  /**/  

}


/* "read_command_line" checks for the input arguments. */
void read_command_line(int argc,char *argv[],int *prm,double *mlo,double *mhi,
		       double *cut,double *mscl,int *all,double *gam)
{
  char *param;
  double x;
  int i,ok=0,nopt=0;

  /* Assign default values */
  *prm=-1;
  *cut=DEFCUT;
  *mscl=DEFMSCL;
  *mlo=DEFMLO;
  *mhi=DEFMHI;
  *all=0;
  *gam=DEFGAM;
  

  /* Allowed flags:
     -p    Print the projection matrix
     -a    Use all residues instead of the specified blocks
     -s    Specify membrane scaling factor
     -c    Specify cutoff distance 
     -b    Specify membrane boundaries
     -k    Specify force constant
  */

  if(argc>1 && argc<=13){


    /* Assign values from blockfile, if any */
    param=(char *)malloc((size_t) (99*sizeof(char)));
    if((param=get_param(argv[1],"mlo"))!=NULL)
      sscanf(param,"%lf",mlo);
    if((param=get_param(argv[1],"mhi"))!=NULL)
      sscanf(param,"%lf",mhi);
    if((param=get_param(argv[1],"cut"))!=NULL)
      sscanf(param,"%lf",cut);
    if((param=get_param(argv[1],"mscl"))!=NULL)
      sscanf(param,"%lf",mscl);

    ok=1;
    i=2;
    while(i<argc){
      if(strcmp(argv[i],"-p")==0){       /* Print mode */
	*prm=1;
	nopt++;
      }
      else if(strcmp(argv[i],"-a")==0){  /* Forget blocks; use all residues */
	*all=1;
	nopt++;
      }
      else if(strcmp(argv[i],"-s")==0){  /* Membrane scaling factor */
	sscanf(argv[++i],"%lf",mscl);
	nopt++;
      }
      else if(strcmp(argv[i],"-c")==0){  /* Cutoff distance */
	sscanf(argv[++i],"%lf",cut);
	nopt++;
      }
      else if(strcmp(argv[i],"-b")==0){  /* Membrane boundaries */
	sscanf(argv[++i],"%lf",mlo);
	sscanf(argv[++i],"%lf",mhi);
	if(*mlo > *mhi){
	  x=*mhi;
	  *mhi=*mlo;
	  *mlo=x;
	}
	nopt++;
      }
      else if(strcmp(argv[i],"-k")==0){  /* Force constant */
	sscanf(argv[++i],"%lf",gam);
	nopt++;
      }
      else{
	fprintf(stderr,"\n%s: Unknown argument: %s\n\n",argv[0],argv[i]);
	exit(1);}
      i++;
    }
    if(nopt>6) ok=0;
  }
  if(ok==0){
    fprintf(stderr,"\nUsage:\n%s blockfile.blk [OPTIONS]\n\n",argv[0]);
    fprintf(stderr,"OPTIONS:\n%s%s%s%s%s%s",
	    "\t-p\tPrint the projection matrix instead of Hessian\n",
	    "\t-a\tDisregard blocks and use all residues (standard ANM)\n",
	    "\t-s MSCL\tSpecify z-scaling factor in membrane\n",
	    "\t-c CUT\tSpecify cut-off distance\n",
	    "\t-k GAM\tSpecify force constant\n",
	    "\t-b MLO MHI\tSpecify membrane boundaries\n\n");
    fprintf(stderr,
	    "Output:\nPrints Hessian (default) or RTB projection matrix\n\n");
    exit(1);
  }
  fprintf(stderr,"%s%d\n%s%d\n%s%f\n%s%f\n%s%f\n%s%f\n%s%d\n%s%f\n",
	  "Print: ",*prm,"All: ",*all,"mscl: ",*mscl,"mlo: ",*mlo,
	  "mhi: ",*mhi,"cut: ",*cut,"nopt: ",nopt,"gam: ",*gam);
  return;
}
  

/* "calc_blessian_mem" is the membrane version of dblock_hessian5: 
   in it, 'hess_superrow1' is replaced with 'hess_superrow_mem'. */
int calc_blessian_mem(PDB_File *PDB,dSparse_Matrix *PP1,int nres,int nblx,int elm,double **HB,double gam)
{
  dSparse_Matrix *PP2;
  double **HR,***HT;
  int **CT,*BST1,*BST2;
  int ii,i,j,k,p,q,q1,q2,ti,tj,bi,bj,sb,nc,out;


  /* ------------------- INITIALIZE LOCAL VARIABLES ------------------- */
  fprintf(stderr,"calc_blessian_mem: initializing local variables...\n");

  /* HR holde three rows (corresponding to 1 residue) of the full Hessian */
  HR=zero_dmatrix(1,3*nres,1,3);

  /* CT is an array of contacts between blocks */
  CT=unit_imatrix(0,nblx);

  /* Copy PP1 to PP2 and sort by second element */
  PP2=(dSparse_Matrix *)malloc((size_t)sizeof(dSparse_Matrix));
  PP2->IDX=imatrix(1,elm,1,2);
  PP2->X=dvector(1,elm);
  copy_dsparse(PP1,PP2,1,elm);
  dsort_PP2(PP2,elm,2);

  /* BST1: for all j: BST1[i]<=j<BST[i+1], PP1->IDX[j][1]=i */
  /* BST2: for all j: BST2[i]<=j<BST2[i+1], PP2->IDX[j][2]=i */
  BST1=ivector(1,3*nres+1);
  BST2=ivector(1,6*nblx+1);
  init_bst(BST1,PP1,elm,3*nres+1,1);
  init_bst(BST2,PP2,elm,6*nblx+1,2);
  /* ------------------- LOCAL VARIABLES INITIALIZED ------------------ */



  /* ------------- FIND WHICH BLOCKS ARE IN CONTACT --------------- */
  fprintf(stderr,"calc_blessian_mem: finding block contacts...\n");
  nc=find_contacts1(CT,PDB,nres,nblx,CUT);


  /* Allocate a tensor for the block Hessian */
  HT=zero_d3tensor(1,nc,1,6,1,6);


  /* Calculate each super-row of the full Hessian */
  fprintf(stderr,"calc_blessian_mem: calculating full hessian...\n");
  for(ii=1;ii<=nres;ii++){

    if(PDB->atom[ii].model!=0){

      /* ----------------- FIND SUPER-ROW OF FULL HESSIAN --------------- */
      hess_superrow_mem(HR,CT,PDB,nres,ii,CUT,gam);

      /* Update elements of block hessian */
      q1=BST1[3*(ii-1)+2];
      q2=BST1[3*(ii-1)+3];
      /* Sum over elements of projection matrix corresponding to residue ii:
	 for each k in the following loop, PP1->IDX[k][1]==3*ii + 0,1,2 */
      for(k=BST1[3*ii-2];k<BST1[3*ii+1];k++){
	if(k<q1) q=1;
	else if(k<q2) q=2;
	else q=3;
	i=PP1->IDX[k][2];
	bi=(i-1)/6+1;
	ti=i-6*(bi-1);
	/* Sum over all elements of projection matrix with column j>=i */
	for(p=BST2[i];p<=elm;p++){
	  j=PP2->IDX[p][2];
	  bj=(j-1)/6+1;
	  sb=CT[bi][bj];
	  if(i<=j && sb!=0){  /* the first condition should ALWAYS hold */
	    tj=j-6*(bj-1);
	    HT[sb][ti][tj]+=(PP1->X[k]*PP2->X[p]*HR[PP2->IDX[p][1]][q]);
	  }
	}
      }
    }
  }


  /* Print the block Hessian in sparse format */
  fprintf(stderr,"calc_blessian_mem: projecting into block space...\n");
  out=bless_from_tensor(HB,HT,CT,nblx);

  /* Free up memory */
  free_dmatrix(HR,1,3*nres,1,3);
  free_d3tensor(HT,1,nc,1,6,1,6);
  free_imatrix(CT,0,nblx,0,nblx);
  free_ivector(BST1,1,3*nres+1);
  free_ivector(BST2,1,6*nblx+1);
  free_imatrix(PP2->IDX,1,elm,1,2);
  free_dvector(PP2->X,1,elm);
  return out;
}


/* "hess_superrow_mem" calculates the 'who'-th super-row 
   of the Hessian, using 'cut' as the cutoff and 'gam' as the 
   spring constant for all interactions. */
void hess_superrow_mem(double **HR,int **CT,PDB_File *PDB,int nres,int who,double cut,double gam)
{
  int i,j,k,jj;
  double DX[3],csq=cut*cut,dsq,df,scl=1.0;

  /* Clear the diagonal super-element */
  for(i=1;i<=3;i++)
    for(j=1;j<=3;j++)
      HR[3*(who-1)+i][j]=0.0;

  /* Calculate the submatrices */
  for(jj=1;jj<=nres;jj++){

    /* scalefunc0: Uniform z-scaling inside membrane */
    scl=scalefunc0(PDB->atom[who].X[2],PDB->atom[jj].X[2]);

    if(jj!=who && PDB->atom[jj].model!=0 && 
       CT[PDB->atom[who].model][PDB->atom[jj].model]!=0){
      dsq=0.0;
      for(k=0;k<3;k++){
	DX[k] = (double)PDB->atom[who].X[k] - PDB->atom[jj].X[k];
	dsq+=(DX[k]*DX[k]);
      }


      if(dsq<csq){
	for(i=1;i<=3;i++){
	  for(j=i;j<=3;j++){

	    df=gam*DX[i-1]*DX[j-1]/dsq;

	    /* Strong backbone bonds */
	    if((int)fabs(PDB->atom[who].resnum-PDB->atom[jj].resnum)==1 && 
	       PDB->atom[who].chain==PDB->atom[jj].chain)
	      df*=100.0;


	    /* -------- MEMBRANE RULES -------- */
	    /* Scale lateral components */
	    if(i!=3) df*=scl;
	    if(j!=3) df*=scl;

	    /* Off-diagonal super-elements */
	    HR[3*(jj-1)+i][j]=HR[3*(jj-1)+j][i]=-df;
	      
	    /* Diagonal super-elements */
	    HR[3*(who-1)+i][j]+=df;
	    if(i!=j)
	      HR[3*(who-1)+j][i]+=df;
	  }
	}
      } /* <----- if(dsq<csq) */
      else
	for(i=1;i<=3;i++)
	  for(j=1;j<=3;j++)
	    HR[3*(jj-1)+i][j]=HR[3*(jj-1)+j][i]=0.0;
    } /* <---- if(jj!=who &&...) */


    /* **** THIS IS PROBABLY NOT NECESSARY, AS THESE ELEMENTS ARE LIKELY
       NOT USED IN PROJECTING INTO THE BLOCK SPACE **** */
    else if(jj!=who && PDB->atom[jj].model!=0)
      for(i=1;i<=3;i++)
	for(j=1;j<=3;j++)
	  HR[3*(jj-1)+i][j]=HR[3*(jj-1)+j][i]=0.0;
  }
}


/* "bless_to_hess" CONVERTS A BLOCK HESSIAN BACK INTO 
   A FULL HESSIAN USING THE PROJECTION MATRIX PP */
void bless_to_hess(double **HB,double **HF,dSparse_Matrix *PP,int nres,int elm)
{
  int *I1,*I2,ii,jj,i,j,a,b,max=0;

  /* Get a list of indices */
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


  /* Zero the matrix */
  for(i=1;i<=3*nres;i++)
    for(j=i;j<=3*nres;j++)
      HF[i][j]=HF[j][i]=0.0;


  /* ----- Calculate full Hessian ----- */
  for(ii=1;ii<=elm;ii++)
    for(jj=1;jj<=elm;jj++){

      i=PP->IDX[ii][1];
      a=I2[PP->IDX[ii][2]];
      j=PP->IDX[jj][1];
      b=I2[PP->IDX[jj][2]];
      HF[i][j]+=(PP->X[ii]*PP->X[jj]*HB[a][b]);
      HF[j][i]=HF[i][j];
    }
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}


/* "scalefunc0" returns the factor by which the x- and y- rows of 
   the Hessian superelement are to be scaled.  It returns 1 if neither
   residue is in the membrane, sqrt(MSCL) if both residues are in the 
   membrane, and MSCL^(1/4) if only one residue is in the membrane.
 */
double scalefunc0(double z1,double z2)
{
  double scl=1.0;
  static double s0;
  static int first=1;

  /* CALCULATE THE CONSTANTS ON FIRST CALL */
  if(first==1){
    first=0;
    s0=sqrt(MSCL);
    s0=sqrt(s0);
  }

  if(z1<MHI && z1>MLO) scl*=s0;
  if(z2<MHI && z2>MLO) scl*=s0;
  
  return scl;
}


/* "backproject" projects block modes back into the full-residue space.  
   'VEC' is a 'dim' by 'nev' matrix that holds the block vectors in its 
   first 'imx' rows.  The matrix 'PP' has 'elm' elements and projects 
   between block and full spaces. */ 
void backproject(double **VEC,dSparse_Matrix *PP,
		 int dim,int nev,int imx,int elm)
{
  double **VT,dd;
  int *I1,*I2,max,i,j,k,q;
  int imax=0,kmax=0;


  if(imx>dim){
    fprintf(stderr,
	    "\nbackproject: block dimension %d exceeds full dimension %d\n\n",
	    imx,dim);
    exit(1);}


  /* I1[i] = 0 if column 'i' contains only zeros; otherwise I1[i] = i
     I2[i] = rank of 'i' among non-zero columns */
  max=0;
  for(i=1;i<=elm;i++)
    if(PP->IDX[i][2]>max)
      max=PP->IDX[i][2];
  I1=ivector(1,max);
  I2=ivector(1,max);
  for(i=1;i<=max;i++) I1[i]=0;
  for(i=1;i<=elm;i++) I1[PP->IDX[i][2]]=PP->IDX[i][2];
  j=0;
  for(i=1;i<=max;i++){
    if(I1[i]!=0) j++;
    I2[i]=j;
  }


  /* Copy the first 'imx' rows of 'VEC' to a temporary array 
     and clear all of the elements of 'VEC' */
  VT=dmatrix(1,imx,1,nev);
  for(j=1;j<=nev;j++){
    for(i=1;i<=imx;i++){
      VT[i][j]=VEC[i][j];
      VEC[i][j]=0.0;
    }
    for(i=imx+1;i<=dim;i++) VEC[i][j]=0.0;
  }
  

  /* Calculate: VEC_{ij} = \sum_k PP_{ik}*VT_{kj} */
  for(q=1;q<=elm;q++){
    i=PP->IDX[q][1];
    k=I2[PP->IDX[q][2]];
    if(i>imax) imax=i;
    if(k>kmax) kmax=k;
    dd=PP->X[q];
    for(j=1;j<=nev;j++) VEC[i][j]+=dd*VT[k][j];
  }

  free_dmatrix(VT,1,imx,1,nev);
  free_ivector(I1,1,max);
  free_ivector(I2,1,max);
}



/* "dblock_projections2" CALCULATES THE PROJECTION 
   FROM FULL RESIDUE SPACE TO RIGID BLOCK SPACE */
int dblock_projections2(dSparse_Matrix *PP,PDB_File *PDB,int nres,int nblx,int bmx)
{
  double **X,**I,**IC,*CM,*W,**A,**ISQT;
  double x,tr,dd,df;
  int *IDX,nbp,b,i,j,k,ii,jj,aa,bb,elm;


  /* INITIALIZE BLOCK ARRAYS */
  elm=0;
  X=dmatrix(1,bmx,1,3);
  IDX=ivector(1,bmx);
  CM=dvector(1,3);
  I=dmatrix(1,3,1,3);
  IC=dmatrix(1,3,1,3);
  W=dvector(1,3);
  A=dmatrix(1,3,1,3);
  ISQT=dmatrix(1,3,1,3);

  /* CYCLE THROUGH BLOCKS */
  for(b=1;b<=nblx;b++){

    /* CLEAR MATRICES */
    for(j=1;j<=3;j++){
      CM[j]=0.0;
      for(i=1;i<=3;i++) I[i][j]=0.0;
      for(i=1;i<=bmx;i++) X[i][j]=0.0;
    }

    /* STORE VALUES FOR CURRENT BLOCK */
    nbp=0;
    for(i=1;i<=nres;i++){
      if(PDB->atom[i].model==b){
	IDX[++nbp]=i;
	for(j=1;j<=3;j++){
	  x=(double)PDB->atom[i].X[j-1];
	  X[nbp][j]=x;
	  CM[j]+=x;
	}   
      }
    }

    /* TRANSLATE BLOCK CENTER OF MASS TO ORIGIN */
    for(j=1;j<=3;j++) CM[j]/=(double)nbp;
    for(i=1;i<=nbp;i++)
      for(j=1;j<=3;j++)
	X[i][j]-=CM[j];

    /* CALCULATE INERTIA TENSOR */
    for(k=1;k<=nbp;k++){
      dd=0.0;
      for(j=1;j<=3;j++){
	df=X[k][j];
	dd+=df*df;
      }
      for(i=1;i<=3;i++){
	I[i][i]+=(dd-X[k][i]*X[k][i]);
	for(j=i+1;j<=3;j++){
	  I[i][j]-=X[k][i]*X[k][j];
	  I[j][i]=I[i][j];
	}
      }
    }

    /* DIAGONALIZE INERTIA TENSOR */
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++)
	IC[i][j]=I[i][j];
    dsvdcmp(IC,3,3,W,A);
    deigsrt(W,A,3);
    righthand2(W,A,3);
    /* must this be a right-handed coordinate system? */

    /* FIND ITS SQUARE ROOT */
    for(i=1;i<=3;i++)
      for(j=1;j<=3;j++){
	dd=0.0;
	for(k=1;k<=3;k++)
	  dd+=A[i][k]*A[j][k]/sqrt(W[k]);
	ISQT[i][j]=dd;
      }

    /* UPDATE PP WITH THE RIGID MOTIONS OF THE BLOCK */
    tr=1.0/sqrt((float)nbp);
    for(i=1;i<=nbp;i++){

      /* TRANSLATIONS: 3*(IDX[i]-1)+1 = x-COORDINATE OF RESIDUE IDX[i];
	 6*(b-1)+1 = x-COORDINATE OF BLOCK b */
      for(j=1;j<=3;j++){
	elm++;
	PP->IDX[elm][1] = 3*(IDX[i]-1)+j;
	PP->IDX[elm][2] = 6*(b-1)+j;
	PP->X[elm] = tr;
      }

      /* ROTATIONS */
      if(nbp>1){
	for(ii=1;ii<=3;ii++){
	  for(jj=1;jj<=3;jj++){
	    if(jj==1) {aa=2; bb=3;}
	    else if(jj==2) {aa=3; bb=1;}
	    else {aa=1; bb=2;}
	    dd=ISQT[ii][aa]*X[i][bb]-ISQT[ii][bb]*X[i][aa];
	    elm++;
	    PP->IDX[elm][1] = 3*(IDX[i]-1)+jj;
	    PP->IDX[elm][2] = 6*(b-1)+3+ii;
	    PP->X[elm] = dd;
	  }
	}
      }
    }
  }
  free_dmatrix(X,1,bmx,1,3);
  free_ivector(IDX,1,bmx);
  free_dvector(CM,1,3);
  free_dmatrix(I,1,3,1,3);
  free_dmatrix(IC,1,3,1,3);
  free_dvector(W,1,3);
  free_dmatrix(A,1,3,1,3);
  free_dmatrix(ISQT,1,3,1,3);

  return elm;
}

