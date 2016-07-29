/* Header file for membrane ANM utilities */
#ifndef _MEMBRANMUTIL_H_
#define _MEMBRANMUTIL_H_

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#ifndef M_SQRT1_2
# define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif


#define FREE_ARG char*

#define PDB_MAX_LINE 90
#define NAME_LNG 30  /* Maximum length of a protein name */


/* PDB file-related structures */
typedef struct {int num;char chain;} Resid;
typedef struct {char HEAD[7];int atmnum;char ATOM[5];char RES[4];char chain;int resnum;float X[3];float beta;char ELEMENT[3];int model;} Atom_Line;
typedef struct {char **HEADER;Atom_Line *atom;} PDB_File;


/* Rigid block-related structures */
typedef struct {char RES[4];char chain;int resnum;float disp;} Disp_Iso;
typedef struct {char RES[4];char chain;int resnum;float X[3];} Disp_Aniso;
typedef struct {int blknum;char LORES[4];char lochain;int lonum;char HIRES[4];char hichain;int hinum;} Rigid_Block;


/* Centroid-related structures */
typedef struct {double mass;double radius;double *X;int num;char name[NAME_LNG];double beta;} Centroid;
typedef struct {char name1[NAME_LNG];char name2[NAME_LNG];double mass;} Map;


/* Linear algebra-related structure */
typedef struct {int **IDX;double *X;} dSparse_Matrix;


/* ------------------ Numerical Recipes Macros ------------------ */

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


/* Membrane ANM Functions */
int assign_rigid_blocks(PDB_File *,Rigid_Block *,int,int,int *);
int bless_from_tensor(double **,double ***,int **,int);
Centroid *centroid_vector(int lo,int hi);
void copy_dsparse(dSparse_Matrix *,dSparse_Matrix *,int,int);
char **cmatrix(long nrl,long nrh,long ncl,long nch);
void cross(double [],double [],double []);
void dblock_hessian5(PDB_File *,dSparse_Matrix *,int,int,int,double);
void dsort_PP2(dSparse_Matrix *,int,int);
int filerow(char *file);
int find_contacts1(int **,PDB_File *,int,int,double);
void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch);
char *get_param(char *file,char *param);
void hess_superrow1(double **,int **,PDB_File *,int,int,double,double);
void init_bst(int *,dSparse_Matrix *,int,int,int);
int pdb_field_count(char *file,char *str);
void pdb_hmr(char *file,int *nhed,int *nmod,int *nca);
void pdb_init(PDB_File *PDB,int nhed,int nres);
void print_hess_tensor(double ***,int **,int);
void print_prj_ofst(dSparse_Matrix *,int);
void read_blockfile(char *,char **,Rigid_Block **,int *);
int read_centfile(char *file,Centroid *PTS,Map *map,int np,int nmp);
Centroid *read_centroids1(char *pdbfile,char *mapfile,int *np);
int read_map(char *file,Map *map,int np);
int read_pdb3(char *file,PDB_File *PDB,int nhed,int nres);
void read_sparsemtx(char *file,double **MM,int row,int col);
void righthand2(double *,double **,int);
int **unit_imatrix(long,long);
double ***zero_d3tensor(long,long,long,long,long,long);
double **zero_dmatrix(long,long,long,long);

/* ------------------ Numerical Recipes Routines ------------------ */
double dpythag(double a, double b);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void deigsrt(double d[], double **v, int n);
void nrerror(char error_text[]);
int *ivector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector(int *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);

#endif
