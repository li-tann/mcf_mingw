/*
  Copyright 2008 Gamma Remote Sensing AG Switzerland clw/uw
  v2.9 28-Feb-2005 added orbit propagation and state vector calculation routines
  v3.0 11-Nov-2005 phg_est correction to avoid null data values
  v3.1 16-May-2008 added support for LARGEFILES on Win32 (MSDOS)
  v3.2 11-Jun-2008 increased polynomial interpolator order to 8
*/
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef FFTW			/* FFTW v2.1.5 subroutines */
#ifdef SFFTW
#include <sfftw.h>
#include <srfftw.h>
#endif
#ifndef SFFTW
#include <fftw.h>
#include <rfftw.h>
#endif
#endif

#define FFT_F		1	/* forward FFT */
#define FFT_I		-1	/* inverse FFT */

#ifndef DATA_TYPES
#define FCMPLX_DATA 	0
#define SCMPLX_DATA	1
#define FLOAT_DATA	2
#define INT_DATA	3
#define SHORT_DATA	4
#define BYTE_DATA	5
#define RASTER_DATA	6
#define POINT_DATA	7
#define COORD_DATA	8
#define DATA_TYPES      1
#endif

#define PI	3.1415926535897932384
#define TWO_PI  6.2831853071795864769
#define FOUR_PI 12.566370614359172953
#define SQRT2   1.41421356237 		/* square root of 2 */
#define RTD	57.2957795131		/* radians to degrees */
#define DTR	.0174532925199		/* degrees to radians */
#define C	2.99792458e8
#define INERTIAL_TO_EARTH     1		/* inertial to rotating earth-fixed reference frame */
#define EARTH_TO_INERTIAL    -1		/* rotating earth-fixed to inertial reference frame */
/* International Earth Rotation Service http://hpiers.obspm.fr */
#define SIDEREAL_OMEGA	7.29211514672e-5	/* siderial rotation rate */
#define SOLAR_SIDEREAL	1.002737909350795	/* ratio of conventional solar day to sidereal day */
#define G_EARTH		3.986004418e14		/* geocentric grav. const. = G*Mass_earth (m**3/s**2) */
#define SPACE_MIN_ALT	1.5e5			/* minimum altitude to qualify as in orbit */

#define J2		.00108263 		/* spherical harmonics of earth gravitational field */
#define J3              -2.5323078e-06
#define J4              -1.6204300e-06
#define WGS84_MAJOR 	6378137.0		/* Earth semi-major axis */

#define NR_END 1
#define FREE_ARG char*
#define MAX_INTERP_TIME		90.0	/* range of orbit propagation +/- 90.0 sec */
#define PPTS         8			/* number of points to use in the position and velocity interpolation */

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define SQR(a)    ( (a)*(a) )
#define CUBE(a)   ( (a)*(a)*(a) )
#define MAX_STR_LEN	511

static float maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1, minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define MAX(a,b)  ( ( (a) > (b) ) ? (a) : (b) )
#define MIN(a,b)  ( ( (a) < (b) ) ? (a) : (b) )

#ifndef STRUCTURES
typedef struct {
  double re, im;
} dcomplex;	/* double precision complex data type */
typedef struct {
  float re, im;
} fcomplex; /* single precision complex data type */
typedef struct {
  short re, im;
} scomplex;	/* short complex data type */
typedef struct {
  unsigned char re, im;
} bcomplex;  /* byte data complex structure */

typedef struct {
  double x, y, z;
} VEC;	/* vector structure cartesian (X,Y,Z)*/
typedef struct {
  double t, c, n;
} VEC_TCN;	/* vector structure (Track, Cross-Track, Normal) */
typedef struct {
  double s, c, h;
} VEC_SCH;	/* vector structure (Track, Cross-Track, Height) */
typedef struct {
  double n, e, d;
} VEC_NED;	/* vector structure (North, East, Down) */
typedef struct {
  double v, c, p;
} VEC_VCP;	/* vector structure (Velocity, Cross-Track, Perpendicular) */
typedef struct {
  double v, p, q;
} VEC_VPQ;  /* vector structure (V,P,Q) */
typedef struct {
  float az, cr;
} Slope; 	/* azimuth,range slope */

typedef struct {			/* state vector data structure */
  VEC pos;			/* position vector */
  VEC vel;			/* velocity vector */
} STATE;

typedef double **MAT;		/* 2D matrix data structure */

typedef struct {		/* Keplerian Orbital Elemens */
  double a;		/* semimajor axis  (m) */
  double e;		/* eccentricity */
  double inc;		/* inclination: angle between orbit normal and the z-axis
			of the Earth Centered Inertial coordinate system  (rad) */
  double rasc;		/* right ascension of the ascending node (rad) */
  double omega;		/* argument of perigee (rad) */
  double theta;		/* true anomoly (rad) */
  double E;		/* Eccentric anomaly (rad) */
  double M;		/* mean anomaly (rad) */
} KOE;

typedef struct {
  int x, y;
}COORD;

#define STRUCTURES 	1
#endif

char *sstr(char *a, int na, char *s, int ns);			/* find a specific set of characters in a string */
int rd_str(char *sarpar, int len, char *key, char *str);	/* read different data types */
int rd_int(char *sarpar, int len, char *key, int *);
int rd_dbl(char *sarpar, int len, char *key, double *);
int rd_flt(char *sarpar, int len, char *key, float *);
int rd_vec(char *sarpar, int len, char *key, double *dbl, int nv);

double a2_dbl(char *s1, int offset, int len); 			/* convert string extracted from array to double */
int a2_int(char *s1, int offset, int len);			/* convert string extracted from array to int */
int read_int4(char *a, int j);					/* extract a 4 byte short integer from a byte array */
int read_short2(char *a, int j);				/* extract a 2 byte short integer from a byte array */
float read_float(char *a, int j);
double read_double(char *a, int j);

void ***calloc_3d(size_t ncols, size_t nlines, size_t nlayers, size_t size); /* 3D memory array allocation and initializaton */
void **calloc_2d(size_t ncols, size_t nlines, size_t size);	/* 2D memory array allocation and initializaton */
void *calloc_1d(size_t nv, size_t size);			/* 1D memory array allocation and initializaton */

void zero_3d(void ***a, int width, int nlines, int nlayers, size_t size);/* 3D memory initialization */
void zero_2d(void **a, int width, int nlines, size_t size);	/* 2D memory initialization */
void zero_1d(void *a, int width, size_t size);			/* 1D memory initialization */
/* read a block of fcomplex or scomplex data from a file */
int rd_clist(char *fn, COORD **pt);				/* read a list of x,y coordinates in text format */
int rdp_slc(fcomplex **a, int width, int nls, int roff, int azoff, int nr, int naz, int dtype, double sc, FILE *fltf); /* read a block of scomplex or fcomplex data */
int rdp_flt(float **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *fltf); /* read a block of float data */

int julday(int id, int mm, int iyyy);				/* compute Julian day since 1/1/1950 */
void caldat(int julian, int *id, int *mm, int *iyyy);		/* calculate date given julian day */
void polint(double xa[], double ya[], int n, double x, double *y, double *dy); /* polynomial interpolation */
int doy(int id, int mm, int iyyy);				/* calculate day of year */
int doy2jd(int doy1, int iyyy);					/* convert day of year to julian day */

void four1(float data[], unsigned int n, int isign);		/* complex FFT */
void fourn(float data[], unsigned int nn[], int ndim, int isign); /* multi-dimensional complex FFT */
void realft(float data[], unsigned int n, int isign);		/* FFT of real valued data */
double bessi0(double x);					/* modified Bessel function of order 0 */

void lfit1d(double x1[], double y[], double sig[], int ndat, double ac[], int ia[],
            int ma, double **covar, double *chisq, void (*funcs)(double, double *));
void lfit(double x1[], double x2[], double y[], double sig[], int ndat, double ac[], int ia[],
          int ma, double **covar, double *chisq, void (*funcs)(double, double, double *));
void covsrt(double **covar, int ma, int ia[], int mfit);
void gaussj(double **a, int n, double **b, int m);
void hpsort(int n, float ra[]);					/* heap sort */
void svdfit(double x1[], double x2[], double y[], double sig[], int ndata, double a[], int ma,
            double **u, double **v, double w[], double *chisq, void (*funcs)(double, double, double *));
void svdvar(double **v, int ma, double w[], double **cvm);	/* matrix SVD */
void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
void svdcmp(double **a, int m, int n, double w[], double **v);
double pythag(double a, double b);				/* Pythagorean theorem */
void nrerror(char error_text[]);

unsigned char *uvector(int nl, int nh);				/* allocate space for 1-D data vectors */
int *ivector(int nl, int nh);
short *i2vector(int nl, int nh);
float *vector(int nl, int nh);
double *dvector(int nl, int nh);
fcomplex *cvector(int nl, int nh);
scomplex *scvector(int nl, int nh);

void free_uvector(unsigned char *v, int nl, int nh);		/* free vector arrays */
void free_ivector(int *v, int nl, int nh);
void free_i2vector(short *v, int nl, int nh);
void free_vector(float *v, int nl, int nh);
void free_dvector(double *v, int nl, int nh );
void free_cvector(fcomplex *v, int nl, int nh);
void free_scvector(scomplex *v, int nl, int nh );

unsigned char **umatrix(int nrl, int nrh, int ncl, int nch);	/* allocate 2-D matrices */
int **imatrix(int nrl, int nrh, int ncl, int nch);
short **i2matrix(int nrl, int nrh, int ncl, int nch);
float **matrix(int nrl, int nrh, int ncl, int nch);
extern double **dmatrix(int nrl, int nrh, int ncl, int nch);
fcomplex **cmatrix(int nrl, int nrh, int ncl, int nch);
scomplex **scmatrix(int nrl, int nrh, int ncl, int nch);

void free_umatrix(unsigned char **m, int nrl, int nrh, int ncl, int nch);
void free_i2matrix(short **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void free_matrix(float **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void free_cmatrix(fcomplex **m, int nrl, int nrh, int ncl, int nch);
void free_scmatrix(scomplex **m, int nrl, int nrh, int ncl, int nch);

/* vector routines */
void vprn(char *title, VEC *a);			/* print vector */
double dot(VEC *a, VEC *b);			/* dot product a . b */
void cross(VEC *a, VEC *b, VEC *c);		/* cross product c = a X B */
void unit(VEC *a, VEC *u);			/* u is a unit vector made from a */
void add(VEC *a, VEC *b, VEC *c);		/* add vectors c = a + b */
void sub(VEC *a, VEC *b, VEC *c);		/* subtract vectors c = a - b */
void smul(double sm, VEC *a, VEC *b);		/* scalar multiply b = sm * a */
void unit(VEC *a, VEC *n);			/* create unit vector */
void lc(double s1, VEC *a, double s2, VEC *b, VEC *r);	/* linear combination r = s1*a + s2*b */
double norm(VEC *a);				/* norm of a vector */

/* matrix routines */
void vec_mat(VEC *a, VEC *b, VEC *c, MAT z);	/* copy three 3x1 vectors, create 3x3 matrix, vectors are the rows */
void mat_vec(MAT z, VEC *a, VEC *b, VEC *c);	/* copy the rows of a 3x3 matrix into 3  3x1 vectors */
void trans(MAT a, MAT b);			/* transpose matrix a, result is b */
void matmat(MAT a, MAT b, MAT c);		/* multiply two 3x3 matrices, c = a x b */
void matvm(MAT a, VEC *b, VEC *c);		/* multiply a 3x1 vector by a 3x3 matrix to get a 3x1 vector result */
void matprn(char *title, MAT a);		/* print 3x3 matrix */

/* geometry and orbit state vector routines */
double compute_hour_angle(int iyear, int imonth, int iday, double sec);	/* calculate Greenwich hour angle */
void ref_frame(STATE *sv1, STATE *sv2, double gha, int cflg);	/* change between earth fixed and inertial reference frame */
void sv_calc2(STATE *sv1, STATE *sv2, double delta); 		/* prop. of state vectors in an inertial ref. frame using integration of equations of motion */
void sv2kep(VEC *r, VEC *v, KOE *koe);				/* convert state vector to Keplerian orbital elements */
void kep2sv(VEC *r, VEC *v, KOE *koe);				/* convert Keplerian orbitial elements to a state vector */
double radius(double ra, double e2, double hdg, double lat);
void geo_loc(double ra, double e2, double rho, double fd, double fc, double alt, int s_flag, VEC *s, VEC *v, double *lat, double *lon);
void ll_xyz(double ra, double e2, double lat, double lon, double alt, VEC *s);
void xyz_ll(double ra, double e2, VEC *s, double *lat, double *lon, double *alt);
void c_tcn(VEC *pc, VEC *vc, MAT tcn);			/* construct TCN basis vectors */
void euler(double yaw, double pitch, double roll, MAT e);
void ned_calc(double lat, double lon, VEC *nv, VEC *ev, VEC *dv);
void pos_sv(double tpos, STATE *sv, double t_state, double tis, int nstate, VEC *pos);
void vel_sv(double tpos, STATE *sv, double t_state, double tis, int nstate, VEC *vel);

double dop_cen(double *dcf, double r1, double r);						/* 1-D Doppler polynomial evaluation */
double dop_cen_2d(double *fdp, double *fdp_dot, double *fdp_ddot, double dt,  double dr);	/* 2-D Doppler polynomial evaluation */
void bp_filter(double bw, double wc, int nps, int nfft, double beta, fcomplex *bpf);		/* band-pass filter coefficients */
void fp_peak(double x1, double x2, double p[]);						/* model function for determining peak of 2-D polynomial */
void wrp_cpx(fcomplex **tcpx, int ip, int nr, int naz, FILE *cpxf);			/* write out a block of complex array */
int rd_cpx(fcomplex **a, int width, int nls, int roff, int azoff, int nr, int naz, int dat_type, FILE *datf);
void phg_est(fcomplex **cpx, int roff, int azoff, int nr, int naz, double *rpg, double *azpg);	/*estimate a phase gradient */
void phg_apply(fcomplex **cpx, int nr, int naz, double rpg, double azpg);		/* apply a phase gradient to an fcomplex array */
double peak4(fcomplex **cc, int rwin2, int azwin2, int pwin, double *rpl, double *azpl, int iflg); /* estimate peak of a 2D function */
void cpx_interp(fcomplex **d1, fcomplex **d2, int nx1, int ny1, int n_ovr, int iflg);	/* oversample a 2D fcomplex array by n_ovr */

void start_timing();						/* timing routines */
void stop_timing();

#ifndef FILE_DEF
#ifdef MSDOS
#define FOPEN_RDONLY_TEXT   "r"
#define FOPEN_RDONLY_BINARY "rb"
#define FOPEN_RDWR_TEXT     "r+"
#define FOPEN_RDWR_BINARY   "r+b"
#define FOPEN_WRONLY_TEXT   "w"
#define FOPEN_WRONLY_BINARY "wb"
#define SCAN_LONG_INT       "%ld"
#else
#define FOPEN_RDONLY_TEXT   "r"
#define FOPEN_RDONLY_BINARY "r"
#define FOPEN_RDWR_TEXT     "r+"
#define FOPEN_RDWR_BINARY   "r+"
#define FOPEN_WRONLY_TEXT   "w"
#define FOPEN_WRONLY_BINARY "w"
#define O_BINARY            0
#define O_TEXT              0
#define SCAN_LONG_INT       "%d"
#endif
#define FILE_DEF   1
#endif

#ifndef ENDIAN_DEF
/* for CPU with native LITTLE ENDIAN native byte order CPU_LITTLE_END is defined */
/* to select LITTLE ENDIAN byte order for Gamma Software binary files LITTLE_END is defined */
#ifdef CPU_LITTLE_END
#ifdef LITTLE_END
/* Native byte order: LITTLE ENDIAN, selected byte order: LITTLE ENDIAN */
/* native byte order == selected byte order: no swapping for fread and fwrite commands */
#undef fread
#undef fwrite
#else
/* Native byte order: LITTLE ENDIAN, selected byte order: BIG ENDIAN */
/* native byte order != selected byte order: redefinition of fread and fwrite for swapping */
#define fread fr_swap
#define fwrite fw_swap
int fr_swap(char *, int, int, FILE *);
int fw_swap(char *, int, int, FILE *);
#endif
#else
#ifdef LITTLE_END
/* Native byte order: BIG ENDIAN, selected byte order: LITTLE ENDIAN */
/* native byte order != selected byte order: redefinition of fread and fwrite for swapping */
#define fread fr_swap
#define fwrite fw_swap
int fr_swap(char *, int, int, FILE *);
int fw_swap(char *, int, int, FILE *);
#else
/* Native byte order: BIG ENDIAN, selected byte order: BIG ENDIAN */
/* native byte order == selected byte order: redefinition of fread and fwrite for swapping */
#undef fread
#undef fwrite
#endif
#endif
#define ENDIAN_DEF    1
#endif

#ifndef LARGEFILES_DEF
#ifdef LARGEFILES
#if defined SGI || defined MSDOS

#ifdef SGI
#define fseek(a,b,c) fseek64((a),(off_t)(b),(c))
#define ftell(a) ftell64((a))
#endif

#ifdef MSDOS
#define fseek(a,b,c) fseeko64((a),(off64_t)(b),(c))
#define ftell(a) ftello64((a))
#define off_t off64_t
#define fstat(a, b) _fstati64((a), (__stat64*)(b));
#define stat(a, b) _stati64((a), (__stat64*)(b));
#endif

#else
#define fseek(a,b,c) fseeko((a),(off_t)(b),(c))
#define ftell(a) ftello((a))
#endif

#endif
#define LARGEFILES_DEF   1
#endif

#ifndef SUN_MACRO
static double nintarg;
#define nint(a) ( ((nintarg=(a)) >= 0.0 )?(int)(nintarg+0.5):(int)(nintarg-0.5) )
#define aint(a) ((double)(int)(a))
#define log2(a) (log(a)/.693147181)
#define exp2(a) (pow(2.0,a))
#define SUN_MACRO	1
#endif

#ifndef STRUCTURES
typedef struct {
  double re, im;
} dcomplex;	/* double precision complex data type */
typedef struct {
  float re, im;
} fcomplex; /* single precision complex data type */
typedef struct {
  short re, im;
} scomplex;	/* short complex data type */
typedef struct {
  unsigned char re, im;
} bcomplex;  /* byte data complex structure */

#define STRUCTURES 	1
#endif

/*** block with SWAP_io routines ***/
#undef REMEMBER_fread
#ifdef fread
#define REMEMBER_fread
#undef fread
#endif

int fr_swap(char *data, int nbyte, int size, FILE *outf) {
  static int i;
  static char tmp[8];

  i = fread((char *)data, nbyte, size, outf);

  switch (nbyte) {
  case 2:	/* swap pair of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*2];
      tmp[1] = data[i*2+1];
      data[i*2] = tmp[1];
      data[i*2+1] = tmp[0];
    }
    break;
  case 4:	/* swap quadruple of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*4];
      tmp[1] = data[i*4+1];
      tmp[2] = data[i*4+2];
      tmp[3] = data[i*4+3];
      data[i*4]  = tmp[3];
      data[i*4+1] = tmp[2];
      data[i*4+2] = tmp[1];
      data[i*4+3] = tmp[0];
    }
    break;
  case 8:	/* swap octuple of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*8];
      tmp[1] = data[i*8+1];
      tmp[2] = data[i*8+2];
      tmp[3] = data[i*8+3];
      tmp[4] = data[i*8+4];
      tmp[5] = data[i*8+5];
      tmp[6] = data[i*8+6];
      tmp[7] = data[i*8+7];
      data[i*8]  = tmp[7];
      data[i*8+1] = tmp[6];
      data[i*8+2] = tmp[5];
      data[i*8+3] = tmp[4];
      data[i*8+4] = tmp[3];
      data[i*8+5] = tmp[2];
      data[i*8+6] = tmp[1];
      data[i*8+7] = tmp[0];
    }
    break;
  default: /* no swapping for other nbyte values */
    break;
  }

  return(i);
}

#ifdef REMEMBER_fread
#define fread fr_swap
#undef REMEMBER_fread
#endif

#undef REMEMBER_fwrite
#ifdef fwrite
#define REMEMBER_fwrite
#undef fwrite
#endif

int fw_swap(char *data, int nbyte, int size, FILE *outf) {
  static int i;
  static char tmp[8];

  switch (nbyte) {
  case 2:	/* swap pair of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*2];
      tmp[1] = data[i*2+1];
      data[i*2] = tmp[1];
      data[i*2+1] = tmp[0];
    }
    break;
  case 4:	/* swap quadruple of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*4];
      tmp[1] = data[i*4+1];
      tmp[2] = data[i*4+2];
      tmp[3] = data[i*4+3];
      data[i*4]  = tmp[3];
      data[i*4+1] = tmp[2];
      data[i*4+2] = tmp[1];
      data[i*4+3] = tmp[0];
    }
    break;
  case 8:	/* swap octuple of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*8];
      tmp[1] = data[i*8+1];
      tmp[2] = data[i*8+2];
      tmp[3] = data[i*8+3];
      tmp[4] = data[i*8+4];
      tmp[5] = data[i*8+5];
      tmp[6] = data[i*8+6];
      tmp[7] = data[i*8+7];
      data[i*8]  = tmp[7];
      data[i*8+1] = tmp[6];
      data[i*8+2] = tmp[5];
      data[i*8+3] = tmp[4];
      data[i*8+4] = tmp[3];
      data[i*8+5] = tmp[2];
      data[i*8+6] = tmp[1];
      data[i*8+7] = tmp[0];
    }
    break;
  default: /* no swapping for other nbyte values */
    break;
  }

  i = fwrite((char *)data, nbyte, size, outf);

  return(i);
}

#ifdef REMEMBER_fwrite
#define fwrite fw_swap
#undef REMEMBER_fwrite
#endif

int swap(char *data, int nbyte, int size) {
  static int i;
  static char tmp[8];

  switch (nbyte) {
  case 2:	/* swap pair of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*2];
      tmp[1] = data[i*2+1];
      data[i*2] = tmp[1];
      data[i*2+1] = tmp[0];
    }
    break;
  case 4:	/* swap quadruple of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*4];
      tmp[1] = data[i*4+1];
      tmp[2] = data[i*4+2];
      tmp[3] = data[i*4+3];
      data[i*4]  = tmp[3];
      data[i*4+1] = tmp[2];
      data[i*4+2] = tmp[1];
      data[i*4+3] = tmp[0];
    }
    break;
  case 8:	/* swap octuple of bytes */
    for (i = 0;i < size;i++) {
      tmp[0] = data[i*8];
      tmp[1] = data[i*8+1];
      tmp[2] = data[i*8+2];
      tmp[3] = data[i*8+3];
      tmp[4] = data[i*8+4];
      tmp[5] = data[i*8+5];
      tmp[6] = data[i*8+6];
      tmp[7] = data[i*8+7];
      data[i*8]  = tmp[7];
      data[i*8+1] = tmp[6];
      data[i*8+2] = tmp[5];
      data[i*8+3] = tmp[4];
      data[i*8+4] = tmp[3];
      data[i*8+5] = tmp[2];
      data[i*8+6] = tmp[1];
      data[i*8+7] = tmp[0];
    }
    break;
  default: /* no swapping for other nbyte values */
    break;
  }

  return(i);
}

int swap_int(int value_in) {
  static union {
    char char_array[sizeof(int)];
    int out;
  } value;

  value.out = value_in;
  swap( (char *)(&value.char_array), sizeof(int), 1);
  return(value.out);
}

float swap_float(float value_in) {
  static union {
    char char_array[sizeof(float)];
    float out;
  } value;

  value.out = value_in;
  swap( (char *)(&value.char_array), sizeof(float), 1);
  return(value.out);
}

double swap_double(double value_in) {
  static union {
    char char_array[sizeof(double)];
    double out;
  } value;

  value.out = value_in;
  swap( (char *)(&value.char_array), sizeof(double), 1);
  return(value.out);
}

short swap_short(short value_in) {
  static union {
    char char_array[sizeof(short)];
    short out;
  } value;

  value.out = value_in;
  swap( (char *)(&value.char_array), sizeof(short), 1);
  return(value.out);
}

fcomplex swap_fcomplex(fcomplex value_in) {
  static union {
    char char_array[sizeof(fcomplex)];
    fcomplex out;
  } value;

  value.out = value_in;
  swap( (char *)(&value.char_array), sizeof(float), 2);
  return(value.out);
}

scomplex swap_scomplex(scomplex value_in) {
  static union {
    char char_array[sizeof(scomplex)];
    scomplex out;
  } value;

  value.out = value_in;
  swap( (char *)(&value.char_array), sizeof(short), 2);
  return(value.out);
}

/*** block with sstr and rd_*** routines ***/

#define MAX(a,b)  ( ( (a) > (b) ) ? (a) : (b) )
#define MIN(a,b)  ( ( (a) < (b) ) ? (a) : (b) )

char *sstr(char *a, int na, char *s, int ns);		/* find a specific set of characters in a string */
int rd_str(char *sarpar, int len, char *key, char *str);
int rd_int(char *sarpar, int len, char *key, int *);
int rd_dbl(char *sarpar, int len, char *key, double *);
int rd_flt(char *sarpar, int len, char *key, float *);
int rd_vec(char *sarpar, int len, char *key, double *dbl, int nv);

int rd_str(char *sarpar, int len, char *key, char *str1) {
  /*
    subroutine to read a specific string from parameter files
    Gamma Remote Sensing v2.1 1-Feb-2000 clw

     0 = successful
    -1 = string not found

  */
  static int ns, nc;
  static char *s1, *s2, eol[3];

  eol[0] = '\n';			/* new line */
  eol[1] = ';';			/* semicolon */
  eol[2] = '\0';

  ns = strlen(key);
  s1 = sstr(sarpar, len, key, ns);

  if (s1 == NULL) {
    printf("WARNING: subroutine rd_str, keyword not found: %s\n", key);
    str1[0] = '\0';
    return -1;
  }
  s1 += ns;

  while ((*s1 == ' ') || (*s1 == '\t') || (*s1 == ':') || (*s1 == '"'))s1++;	/* remove leading spaces, tabs, and colons */
  s2 = strpbrk(s1, eol) - 1;					/* last location of the valid string */
  while ((*s2 == ' ') || (*s2 == '\t') || (*s2 == '\r') || (*s2 == '"'))s2--;	/* remove trailing spaces, tabs, and returns */
  nc = MAX(s2 - s1 + 1, 0);				/* number of characters of the string, check for whitespace or NULL string */
  strncpy(str1, s1, nc);				/* destination, source, number of characters */
  str1[nc] = '\0';				/* null terminate copied string */
  return 0;
}

int rd_int(char *sarpar, int len, char *key, int *val) {
  /*
    subroutine to read an integer from parameter files
    Gamma Remote Sensing v2.1 1-Feb-2000 clw/uw

     0 = successful
    -1 = string not found

  */
  static int ns, nc;
  static char *s1, *s2, eol[3];
  static char s3[80];

  eol[0] = '\n';			/* new line */
  eol[1] = ';';			/* semicolon */
  eol[2] = '\0';
  *val = 0;

  ns = strlen(key);
  s1 = sstr(sarpar, len, key, ns);
  if (s1 == NULL) {
    printf("WARNING: subroutine rd_int, keyword not found: %s\n", key);
    return -1;
  }
  s1 += ns;

  while ((*s1 == ' ') || (*s1 == '\t') || (*s1 == ':') || (*s1 == '"'))s1++;	/* remove leading spaces, tabs, and colons */
  s2 = strpbrk(s1, eol) - 1;					/* last location of the valid string */
  while ((*s2 == ' ') || (*s2 == '\t') || (*s2 == '\r') || (*s2 == '"'))s2--;	/* remove trailing spaces, tabs, and returns */
  nc = MAX(s2 - s1 + 1, 0);				/* number of characters of the string, check for null string */

  if (nc > 0) {
    strncpy(s3, s1, nc);				/* destination, source , number of characters */
    s3[nc] = '\0';					/* null terminate copied string */
    if (sscanf(s3, "%d", val) != 1) {
      printf("WARNING rd_int: corrupted data value with key: %s  value: %s\n", key, s3);
      return -1;
    }
    return 0;
  } else {
    printf("WARNING rd_int: no value given for data with key: %s\n", key);
    return -1;
  }
}

int rd_dbl(char *sarpar, int len, char *key, double *dbl_val) {
  /*
    subroutine to read a specific parameter of type double from parameter files
    Gamma Remote Sensing v2.1 1-Feb-2000 clw/uw

     0 = successful
    -1 = string not found

  */
  static int ns, nc;
  static char *s1, *s2, eol[3];
  static char s3[80];

  eol[0] = '\n';			/* new line */
  eol[1] = ';';			/* semicolon */
  eol[2] = '\0';
  *dbl_val = 0.0;

  ns = strlen(key);
  s1 = sstr(sarpar, len, key, ns);
  if (s1 == NULL) {
    printf("WARNING: subroutine rd_dbl, keyword not found: %s\n", key);
    return -1;
  }
  s1 += ns;

  while ((*s1 == ' ') || (*s1 == '\t') || (*s1 == ':') || (*s1 == '"'))s1++;	/* remove leading spaces, tabs, and colons */
  s2 = strpbrk(s1, eol) - 1;					/* last location of the valid string */
  while ((*s2 == ' ') || (*s2 == '\t') || (*s2 == '\r') || (*s2 == '"'))s2--;	/* remove trailing spaces, tabs, and returns */
  nc = MAX(s2 - s1 + 1, 0);				/* number of characters of the string, check for null string */

  if (nc > 0) {
    strncpy(s3, s1, nc);				/* destination, source , number of characters */
    s3[nc] = '\0';				/* null terminate copied string */
    if (sscanf(s3, "%le", dbl_val) != 1) {
      printf("WARNING rd_dbl: corrupted data value with key: %s  value: %s\n", key, s3);
      return -1;
    }
    return 0;
  } else {
    printf("WARNING rd_dbl: no value given for data with key: %s\n", key);
    return -1;
  }
}

int rd_flt(char *sarpar, int len, char *key, float *flt_val) {
  /*
    subroutine to read a specific parameter of type float from parameter files
    Gamma Remote Sensing v2.1 1-Feb-2000 clw/uw

     0 = successful
    -1 = string not found

  */
  static int ns, nc;
  static char *s1, *s2, eol[3];
  static char s3[80];

  eol[0] = '\n';			/* new line */
  eol[1] = ';';			/* semicolon */
  eol[2] = '\0';
  *flt_val = 0.0;

  ns = strlen(key);
  s1 = sstr(sarpar, len, key, ns);
  if (s1 == NULL) {
    printf("WARNING: subroutine rd_flt, keyword not found: %s\n", key);
    return -1;
  }
  s1 += ns;

  while ((*s1 == ' ') || (*s1 == '\t') || (*s1 == ':') || (*s1 == '"'))s1++;	/* remove leading spaces, tabs, and colons */
  s2 = strpbrk(s1, eol) - 1;					/* last location of the valid string */
  while ((*s2 == ' ') || (*s2 == '\t') || (*s2 == '\r') || (*s2 == '"'))s2--;	/* remove trailing spaces, tabs, and returns */
  nc = MAX(s2 - s1 + 1, 0);				/* number of characters of the string, check for null string */

  if (nc > 0) {
    strncpy(s3, s1, nc);				/* destination, source , number of characters */
    s3[nc] = '\0';					/* null terminate copied string */
    if (sscanf(s3, "%e", flt_val) != 1) {
      printf("WARNING rd_flt: corrupted data value with key: %s  value: %s\n", key, s3);
      return -1;
    }
    return 0;
  } else {
    printf("WARNING rd_flt: no value given for data with key: %s\n", key);
    return -1;
  }
}

#define MAX_STR  256

int rd_vec(char *sarpar, int len, char *key, double *dbl, int nv) {
  /*
    subroutine to read an array of doubles (type VEC) from the SAR parameter file
    Gamma Remote Sensing v2.1 1-Feb-2000 clw/uw

     0 = successful
    -1 = string not found

    sarpar	parameter data as a character array
    len		length of parameter array in bytes
    key		search key
    dbl		array of doubles
    nv		number of doubles to read
  */

  int ns, nc, j;
  char *s1, *s2, *s3, *s4;
  char eol[3], sep[4];

  s3 = (char *)malloc(MAX_STR);
  if (s3 == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error in subroutine rd_vec\n\n");
    exit(-1);
  }

  eol[0] = '\n';			/* new line */
  eol[1] = ';';			/* semicolon */
  eol[2] = '\0';
  sep[0] = ' ';			/* valid separators */
  sep[1] = ',';
  sep[2] = '\t';			/* tab */
  sep[3] = '\0';
  for (j = 0; j < nv; j++)dbl[j] = 0.0;		/* initialize data */

  ns = strlen(key);
  s1 = sstr(sarpar, len, key, ns);
  if (s1 == NULL) {
    printf("WARNING: subroutine rd_vec, keyword not found: %s\n", key);
    return -1;
  }
  s1 += ns;

  while ((*s1 == ' ') || (*s1 == '\t') || (*s1 == ':') || (*s1 == '"'))s1++;	/* remove leading spaces, tabs, and colons */
  s2 = strpbrk(s1, eol) - 1;					/* last location of the valid string */
  while ((*s2 == ' ') || (*s2 == '\t') || (*s2 == '\r') || (*s2 == '"'))s2--;	/* remove trailing spaces, tabs, and returns */
  nc = s2 - s1 + 1;					/* number of characters of the string */
  strncpy(s3, s1, nc);				/* destination, source, number of characters */
  s3[nc] = '\0';				/* null terminate copied string */

  for (j = 0; j < nv; j++) {
    s4 = strtok(s3, sep);
    s3 = (char *)NULL;				/* continue scanning for next token */
    if (s4 == NULL) {
      printf("WARNING rd_vec: insufficient number of values: %d    key: %s\n", nv, key);
      return -1;
    }
    if (sscanf(s4, "%le", dbl + j) != 1) {
      printf("WARNING rd_vec: invalid entry for data with key: %s    value: %s\n", key, s4);
      return -1;
    }
  }
  free(s3);
  return 0;
}

char *sstr(char *a, int na, char *s, int ns)
/*
	subroutine to find the occurance of the character array s
        with length ns in the character array a with length na
	clw 21-Oct-93
*/
{
  int i, j;
  for (i = 0; i <= (na - ns); i++) {
    if (a[i] == s[0]) {		/* check first character */
      j = 1;
      while ((j < ns) && (a[i+j] == s[j])) j++;
      if (j == ns) return &a[i];
    }
  }
  return NULL;			/* not found! */
}

void ***calloc_3d(size_t ncols, size_t nlines, size_t nlayers, size_t size) {
  /*
    allocate a 3D memory array as a single contiguous block and initialized to 0

    ncols		number of columns in the 2-D arrays
    nlines	number of lines in each 2-D array
    nlayers       number of 2D layers
    size		size of each array element in bytes

    return value is the address of an array of pointers referencing the
    2D layers of the the data

    24-Jan-2002 clw
  */

  void ***a, **a1, *a2;
  int i, j;

  if (ncols <= 0) {
    fprintf(stderr, "\nERROR: calloc_3d: number of columns <= 0: %d\n", ncols);
    exit(-1);
  }
  if (nlines <= 0) {
    fprintf(stderr, "\nERROR: calloc_3d: number of lines <= 0: %d\n", nlines);
    exit(-1);
  }
  if (nlayers <= 0) {
    fprintf(stderr, "\nERROR: calloc_3d: number of layers <= 0: %d\n", nlayers);
    exit(-1);
  }
  if (size <= 0) {
    fprintf(stderr, "\nERROR: calloc_3d: element size <= 0: %d\n", size);
    exit(-1);
  }

  a2 = (void *)calloc(ncols * nlines * nlayers, size);
  if (a2 == NULL) {
    fprintf(stderr, "\nERROR calloc_3d: memory allocation error for data values\n\n");
    exit(-1);
  }

  a1 = (void **)calloc(nlines * nlayers, sizeof(void *));
  if (a2 == NULL) {
    fprintf(stderr, "\nERROR calloc_3d: memory allocation error for line pointers\n\n");
    exit(-1);
  }

  a = (void ***)malloc(nlayers * sizeof(void **));
  if (a == NULL) {
    fprintf(stderr, "\nERROR calloc_3d: memory allocation error for layer pointers\n\n");
    exit(-1);
  }

  for (i = 0; i < (int)nlayers; i++) {
    a[i] = a1 + i * nlines;		/* pointers to the individual layers */
    for (j = 0; j < (int)nlines; j++)a1[j+i*nlines] = (unsigned char *)a2 + ((i * nlines + j) * ncols) * size; /* pointers to the lines of the 2D array */
  }

  return a;
}

/*** memory allocation and clear memory routines ***/

void **calloc_2d(size_t ncols, size_t nlines, size_t size) {
  /*
    allocate a 2D memory as a single contiguous block and initialized to 0

    ncols		number of columns in the array
    nlines	number of lines
    size		size of each element in bytes

    return value is the address of an array of pointers referencing the
    lines of the 2D array

    24-Jan-2002 clw
  */

  void **a, *a1;
  int i;

  if (ncols <= 0) {
    fprintf(stderr, "\nERROR: calloc_2d: number of columns <= 0: %d\n", ncols);
    exit(-1);
  }
  if (nlines <= 0) {
    fprintf(stderr, "\nERROR: calloc_2d: number of lines <= 0: %d\n", nlines);
    exit(-1);
  }
  if (size <= 0) {
    fprintf(stderr, "\nERROR: calloc_2d: element size <= 0: %d\n", size);
    exit(-1);
  }

  a1 = (void *)calloc(ncols * nlines, size);
  if (a1 == NULL) {
    fprintf(stderr, "\nERROR calloc_2d: memory allocation error for data values\n\n");
    exit(-1);
  }

  a = (void **)malloc(nlines * sizeof(void *));
  if (a == NULL) {
    fprintf(stderr, "\nERROR calloc_2d: memory allocation error for line pointers\n\n");
    exit(-1);
  }

  for (i = 0; i < (int)nlines; i++)a[i] = (unsigned char *)a1 + i * ncols * size;

  return a;
}

void *calloc_1d(size_t nv, size_t size) {
  /*
    allocate a 1D memory array as a single contiguous block and initialized to 0

    nv		number of values
    size		size of each element in bytes

    return value is the address of the array

    24-Jan-2002 clw
  */

  void *a;

  if (nv <= 0) {
    fprintf(stderr, "\nERROR: calloc_1d: number of elements <= 0: %d\n", nv);
    exit(-1);
  }
  if (size <= 0) {
    fprintf(stderr, "\nERROR: calloc_1d: element size <= 0: %d\n", size);
    exit(-1);
  }

  a = (void *)calloc(nv, size);
  if (a == NULL) {
    fprintf(stderr, "\nERROR calloc_1d: memory allocation error for data values\n\n");
    exit(-1);
  }

  return a;
}

void zero_3d(void ***a, int width, int nlines, int nlayers, size_t size) {
  int i, j;

  for (i = 0; i < nlayers; i++) {
    for (j = 0; j < nlines; j++)memset(a[i][j], 0, size*width);
  }
}

void zero_2d(void **a, int width, int nlines, size_t size) {
  int i;
  for (i = 0; i < nlines; i++)memset(a[i], 0, size*width);
}

void zero_1d(void *a, int width, size_t size) {
  memset(a, 0, size*width);
}

int rd_clist(char *fn, COORD **pt) {
  /*
    read list of coordinated entries and return in a data array of coordinates x,y pairs
    26-Jan-2005 clw
  */
  char str[MAX_STR_LEN];
  int npt, ix, iy;
  FILE *cf;

  cf = fopen(fn, FOPEN_RDONLY_TEXT);
  if (cf == NULL) {
    fprintf(stderr, "\nERROR: cannot open clist file: %s\n\n", fn);
    exit( -1);
  }
  npt = 0;

  while (1) { 					/* determine number of records */
    fgets(str, MAX_STR_LEN, cf);
    if (feof(cf) != 0) break;
    if (str[0] == '#')continue;
    if (sscanf(str, "%d %d", &ix, &iy) != 2)continue;
    npt++;
  }

  if (npt == 0) {
    fprintf(stderr, "ERROR: no valid points in the clist: %s\n", fn);
    exit(-1);
  }
  printf("reading clist file: %s   number of points: %d\n", fn, npt);
  fflush(stdout);
  *pt = calloc_1d(npt, sizeof(COORD));
  rewind(cf);
  npt = 0;

  while (1) { 					/* read ITAB records */
    fgets(str, MAX_STR_LEN, cf);
    if (feof(cf) != 0)break;
    if (str[0] == '#')continue;			/* ignore lines starting with # */
    if (sscanf(str, "%d %d", &ix, &iy) != 2)continue;

    (*pt)[npt].x = ix;
    (*pt)[npt].y = iy;
#ifdef DEBUG
    fprintf(stderr, "%7d  %d  %d\n", npt, (*pt)[npt].x, (*pt)[npt].y);
#endif
    npt++;
  }

  fclose(cf);
  return npt;
}

int rdp_slc(fcomplex **a, int width, int nls, int roff, int azoff, int nr, int naz, int dtype, double sc, FILE *fltf) {
  /*
     read fcomplex or scomplex data file patch from a region starting at (line:azoff1, sample: roff1).
     nominal patch size is naz lines x nr samples. The return value is lines actually read from the file
     and is set to 0 if no lines read. All entries in the fcomplex data array are set to (0.0,0.0)
     if no data read for that location
   
     30-Sep-2003 clw
   
     a		output 2D fcomplex point phase derived from complex data
     width	width of input data file, number of samples/line
     nls	number of lines in the input fcomplex data
     roff       range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz        number of azimuth lines in the patch
     dtype      SLC type 0=fcomplex 1=scomplex
     sc         SLC scale factor
     fltf	complex floating point format FILE data structure
   
  */
  int i, j, jj, naz1, nr1, err;
  fcomplex *xtmp;
  scomplex *stmp;

  /*  printf("rdp_slc width:%d  nls:%d  roff:%d  azoff:%d  nr:%d  naz:%d  dtype:%d\n",width,nls,roff,azoff,nr,naz,dtype); */

  if (roff >= width || azoff >= nls)
    return (0);
  if ((roff < 0) || (azoff < 0))
    return (0);

  xtmp = (fcomplex *)calloc_1d(width, sizeof(fcomplex));
  stmp = (scomplex *)calloc_1d(width, sizeof(scomplex));

  if ((naz + azoff) > nls)	/* check if blank lines need to be written at the end of the array */
    naz1 = nls - azoff;
  else
    naz1 = naz;
  if ((nr + roff) > width)	/* check if blank samples need to be written at the end of the lines */
    nr1 = width - roff;
  else
    nr1 = nr;

  if (dtype == 0)		/* seek to start of line */
    err = fseek(fltf, (off_t)sizeof(fcomplex) * azoff * width, SEEK_SET);
  else
    err = fseek(fltf, (off_t)sizeof(scomplex) * azoff * width, SEEK_SET);

  for (i = 0; i < naz1; i++) {
    if (dtype == 0) {
      fread((char *)xtmp, sizeof(float), 2*width, fltf);
      if (feof(fltf)) {
        fprintf(stderr, "\nERROR rdp_slc: unexpected end of data file at line: %d\n", azoff + i);
        exit( -1);
      }
      for (j = 0; j < nr1; j++) {
        jj = j + roff;			/* offset in range */
        a[i][j].re = xtmp[jj].re * sc;
        a[i][j].im = xtmp[jj].im * sc;
      }
    } else {
      fread((char *)stmp, sizeof(short), 2*width, fltf);
      if (feof(fltf)) {
        fprintf(stderr, "\nERROR rdp_slc: unexpected end of data file at line: %d\n", azoff + i);
        exit( -1);
      }
      for (j = 0; j < nr1; j++) {
        jj = j + roff;			/* offset in range */
        a[i][j].re = stmp[jj].re * sc;
        a[i][j].im = stmp[jj].im * sc;
      }
    }
    if (nr1 < nr) {			/* fill in missing range samples */
      for (j = nr1; j < nr; j++) {
        a[i][j].re = 0.0;
        a[i][j].im = 0.0;
      }
    }
  }

  if (naz1 < naz) {			/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++) {
        a[i][j].re = 0.0;
        a[i][j].im = 0.0;
      }
    }
  }

  free(xtmp);
  free(stmp);
  return (naz1);
}

int rdp_flt(float **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *fltf) {
  /*
     read float data file patch from a region starting at (line:azoff1, sample: roff1).
     nominal patch size is naz lines x nr samples. The return value is lines actually read from the file
     and is set to 0 if no lines read. All entries in the fcomplex data array are set to (0.0)
     if no data read for that location
   
     10-Mar-2003 clw
   
     a		output 2D fcomplex point phase derived from complex data
     width	width of input data file, number of samples/line
     nls	number of lines in the input float data
     roff       range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz        number of azimuth lines in the patch
     fltf	complex floating point format FILE data structure
   
  */
  int i, j, jj, naz1, nr1, err;
  float *ftmp;

  if (roff >= width || azoff >= nls)return (0);
  if ((roff < 0) || (azoff < 0))return (0);

  ftmp = (float *)calloc_1d(width, sizeof(float));

  if ((naz + azoff) > nls)	/* check if blank lines need to be written at the end of the array */
    naz1 = nls - azoff;
  else
    naz1 = naz;
  if ((nr + roff) > width)	/* check if blank samples need to be written at the end of the lines */
    nr1 = width - roff;
  else
    nr1 = nr;

  err = fseek(fltf, (off_t)sizeof(float) * azoff * width, SEEK_SET);

  for (i = 0; i < naz1; i++) {
    fread((char *)ftmp, sizeof(float), width, fltf);
    if (feof(fltf)) {
      fprintf(stderr, "\nERROR rdp_flt: unexpected end of data file at line: %d\n", azoff + i);
      exit( -1);
    }
    for (j = 0; j < nr1; j++) {
      jj = j + roff;			/* offset in range */
      a[i][j] = ftmp[jj];
    }

    if (nr1 < nr) {			/* fill in missing range samples */
      for (j = nr1; j < nr; j++) {
        a[i][j] = 0.0;
      }
    }
  }

  if (naz1 < naz) {			/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++) {
        a[i][j] = 0.0;
      }
    }
  }

  free(ftmp);
  return (naz1);
}

void cross(VEC *a, VEC *b, VEC *c) {
  /**** cross product c = a X b ****/
  c->x = a->y * b->z - a->z * b->y;
  c->y = a->z * b->x - a->x * b->z;
  c->z = a->x * b->y - a->y * b->x;
}

void add(VEC *a, VEC *b, VEC *c) {
  /**** vector sum c = a + b ****/
  c->x = a->x + b->x;
  c->y = a->y + b->y;
  c->z = a->z + b->z;
}

void sub(VEC *a, VEC *b, VEC *c) {
  /**** vector difference c = a - b ****/
  c->x = a->x - b->x;
  c->y = a->y - b->y;
  c->z = a->z - b->z;
}

void smul(double sm, VEC *a, VEC *b) {
  /**** vector scalar multiply b = a * sm ****/
  b->x = a->x * sm;
  b->y = a->y * sm;
  b->z = a->z * sm;
}

void unit(VEC *a, VEC *n) {
  /**** unit vector is constructed from vector a ****/
  double mag;
  mag = sqrt(a->x * a->x + a->y * a->y + a->z * a->z);

  n->x = a->x / mag;
  n->y = a->y / mag;
  n->z = a->z / mag;
}

double norm(VEC *a) {
  /**** magnitude of vector a ****/
  double mag;
  mag = sqrt(a->x * a->x + a->y * a->y + a->z * a->z);
  return (mag);
}

double dot(VEC *a, VEC *b) {
  /**** dot product <a,b> ****/
  double sum;
  sum = a->x * b->x;
  sum += a->y * b->y;
  sum += a->z * b->z;
  return (sum);
}

void lc(double s1, VEC *a, double s2, VEC *b, VEC *r) {
  /**** linear combination of two VEC, r = s1*a + s2*b ****/
  r->x = s1 * a->x + s2 * b->x;
  r->y = s1 * a->y + s2 * b->y;
  r->z = s1 * a->z + s2 * b->z;
}

void vprn(char *title, VEC *a) {
  /* print out a vector */
  printf("%s %13.6e %13.6e %13.6e\n", title, a->x, a->y, a->z);
}

void vec_mat( VEC *a, VEC *b, VEC *c, MAT z) {
  /* pack 3 vectors a,b,c as the rows of z */
  z[0][0] = a->x;
  z[0][1] = a->y;
  z[0][2] = a->z;
  z[1][0] = b->x;
  z[1][1] = b->y;
  z[1][2] = b->z;
  z[2][0] = c->x;
  z[2][1] = c->y;
  z[2][2] = c->z;
}

void mat_vec(MAT z, VEC *a, VEC *b, VEC *c) {
  /* unpack rows of z into vectors a,b,c */
  a->x = z[0][0];
  a->y = z[0][1];
  a->z = z[0][2];
  b->x = z[1][0];
  b->y = z[1][1];
  b->z = z[1][2];
  c->x = z[2][0];
  c->y = z[2][1];
  c->z = z[2][2];
}

void trans(MAT a, MAT b) {
  /* transpose 3x3 matrix a, result is in b */
  double tmp;

  b[0][0] = a[0][0];
  b[1][1] = a[1][1];
  b[2][2] = a[2][2];

  tmp = a[1][0]; /* allow b to overwrite a */
  b[1][0] = a[0][1];
  b[0][1] = tmp;

  tmp = a[2][0];
  b[2][0] = a[0][2];
  b[0][2] = tmp;

  tmp = a[1][2];
  b[1][2] = a[2][1];
  b[2][1] = tmp;
}

void matmat(MAT a, MAT b, MAT c) {
  /* multiply c = a x b where a,b,c are 3x3 matrices */
  int i, j;

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
    }
  }
}

void matvm(MAT a, VEC *b, VEC *c) {
  /* multiply a 3x1 vector by a 3x3 matrix to get a 3x1 vector result */

  c->x = a[0][0] * b->x + a[0][1] * b->y + a[0][2] * b->z;
  c->y = a[1][0] * b->x + a[1][1] * b->y + a[1][2] * b->z;
  c->z = a[2][0] * b->x + a[2][1] * b->y + a[2][2] * b->z;
}

void matprn(char *title, MAT a) {
  /* print out a 3x3 matrix */
  printf("\n%s\n", title);
  printf("%10.3e %10.3e %10.3e\n", a[0][0], a[0][1], a[0][2]);
  printf("%10.3e %10.3e %10.3e\n", a[1][0], a[1][1], a[1][2]);
  printf("%10.3e %10.3e %10.3e\n", a[2][0], a[2][1], a[2][2]);
}

double dop_cen(double *dcf, double r1, double r) {
  /*
  subroutine to evaluate Doppler polynomial clw 23-aug-93

    dcf  cubic Doppler polynomial coefficients
    r1   range to center of the raw data swath (meters)
    r    slant range (meters)
  */
  double fd;			/* Doppler centroid */
  double dr;

  dr = r - r1;
  fd = dcf[0] + dr * (dcf[1] + dr * (dcf[2] + dr * dcf[3]));
  return fd;
}

double dop_cen_2d(double *fdp, double *fdp_dot, double *fdp_ddot, double dt,  double dr)
/*
  evaluate Doppler polynomial including along_track clw 11-Jun-2003

  fdp		cubic Doppler polynomial containg slant range info
  fdp_dot	first derivative of the coefficents in fdp
  fdp_ddot      second derivative of the coefficents in fpd
  dt		time relative to the center of the image (along-track)
  dr		slant range relative to the center of the swath
*/
{
  double fd;			/* Doppler centroid */

  fd = fdp[0] + dt * fdp_dot[0] + + dt * dt * fdp_ddot[0] +
       dr * (fdp[1] + fdp_dot[1] * dt + dr * (fdp[2] + dr * fdp[3]));
  return fd;
}

int doy(int id, int mm, int iyyy) {
  /*
    calculate day of year from day, month, year

    7-Jan-1999 clw

  */
  int jd1, jd2, doy1;

  if ((id < 1) || (id > 31)) {
    fprintf(stderr, "\nERROR subroutine doy: error in day of month: %d\n", id);
    exit(-1);
  }
  if ((mm < 1) || (mm > 12)) {
    fprintf(stderr, "\nERROR subroutine doy: error in month: %d\n", mm);
    exit(-1);
  }

  jd1 = julday(1, 1, iyyy);		/* julian day of Jan 1 */
  jd2 = julday(id, mm, iyyy);		/* julian day of date */
  doy1 = jd2 - jd1 + 1;
  if ((doy1 < 1) || (doy1 > 366)) {
    fprintf(stderr, "\nERROR subroutine doy: invalid day of year %d\n", doy1);
    exit(-1);
  }
  return doy1;
}

int doy2jd(int doy1, int iyyy) {
  /*
    calculate julian day from day of year and year

    7-Jan-1999 clw

  */
  int jd1, jd2;

  jd1 = julday(1, 1, iyyy);		/* julian day of Jan 1 */
  jd2 = jd1 + doy1 - 1;			/* julian day of date */
  return jd2;
}

#define WE 	4.1666666666666667e-03	/* degrees/sec 360/86400 */
#define JD2000  2451545.0		/* Julian day of 1-Jan-2000  12h UT */

double compute_hour_angle(int iyear, int imonth, int iday, double sec)
/*
Gamma Remote Sensing AG clw 11-Jun-2004.

inputs: year, month (1-12), and day (int)
       sec: UT seconds of day  (double)
output: Greewich hour angle (0 - 360.0 degrees)

This function is for getting the mean Hour Angle of Greenwich at
input date/time.  The time must be specified in Universal Time,
which is offset from  civil time (UTC) by a daily varying amount
(UT1 - UTC).  Tables of UT1-UTC (Bull. A) are available from  the
Earth rotation service http://hpiers.obspm.fr/eop-pc. or from the
US Naval Observatory http://maia.usno.navy.mil/search/search.html.

The expression for GST is described in: "The New Definition of
Universal Time", S. Aoki, et al. Astronomy and Astrophysics vol 105,
pages 359-361 (1982).  Dates are expressed relative to 1-Jan-2000.
The expression for the Greenwich mean sidereal time  at 0h UT
(GMST1) is given in equation 13:

GMST1 = 24110.54841 + 8640184.812866*Tu + .093104*Tu*Tu -6.210E-06 * Tu*Tu*Tu

Adding a term for the time offset since 0h UT1:

GMST = 24110.54841 + 8640184.812866*Tu + .093104*Tu*Tu -6.210E-06 * Tu*Tu*Tu + sec*1.002737909350795

where Tu is the time in Julian centuries  (days - JD2000)/3625 and JD2000 is the Julian day of 1-Jan-2000
is the number of days that have elapsed since JD 2451545.0 UT1 (2000 Jan 1 12 hours UT1).
The expression here is converted to degrees by scaling by 360./86400.

*/
{
  double day, t;
  double hour_angle, gst;

  day = julday(iday, imonth, iyear) - .5;	/* day has values +/-.5, +/-1.5... */
  t = (day - JD2000) / 36525.0;			/* elapsed time in Julian centuries */
  gst = 24110.54841 + 8640184.812866 * t + .093104 * t * t - 6.210e-6 * t * t * t + sec * SOLAR_SIDEREAL;
  gst = fmod(gst, 86400.0);

  if (gst < 0.0)gst += 86400.;
  hour_angle = gst * WE;

#ifdef DEBUG
  printf("\ndate (day, month, year): %d  %d  %d  sec: %f \n", iday, imonth, iyear, sec);
  printf("JD since 1-Jan-2000, 12h UT1: %f   Julian centuries: %f\n", day - JD2000, t);
  printf("GST: %f  hour angle: %f\n", gst, hour_angle);
#endif

  return (hour_angle);
}

void ref_frame(STATE *sv1, STATE *sv2, double gha, int cflg) {
  /*
    convert state vectors between earth fixed and inertial reference frame
    apply earth rotation angle to get coordiates relative to Greenwich

    sv1	  input state vector
    sv2     output state vector
    gha     Greenwich hour angle (degrees)
    cflg    conversion flag 1: inertial to earth-fixed rotating  -1 earth-fixed rotating to inertial

    10-Apr-2003  clw
  */

  double sin_gha, cos_gha;
  VEC spin, c1;
  VEC p1, p2, v1, v2, v3;

  p1 = sv1->pos;
  v1 = sv1->vel;
  spin.x = 0.0;
  spin.y = 0.0;
  spin.z = SIDEREAL_OMEGA;  	/* rotation rate of Earth (rad/sec) */
  gha *= DTR;						/* convert gha to radians */
  sin_gha = sin(gha);
  cos_gha = cos(gha);

  switch (cflg) {
  case INERTIAL_TO_EARTH:
    p2.z = p1.z;					/* rotate by gha */
    p2.x = p1.x * cos_gha + p1.y * sin_gha;
    p2.y = p1.y * cos_gha - p1.x * sin_gha;
    v2.z = v1.z;					/* rotate by gha */
    v2.x = v1.x * cos_gha + v1.y * sin_gha;
    v2.y = v1.y * cos_gha - v1.x * sin_gha;
    cross(&spin, &p2, &c1);				/* coriolis velocity (spin vector x position) */
    sub(&v2, &c1, &v3);					/* subtract Coriolis velocity */

#ifdef DEBUG
    printf("\n*** inertial to earth-fixed ***\n");
    printf("inertial position x,y,z:         %12.3f %12.3f %12.3f\n", p1.x, p1.y, p1.z);
    printf("inertial velocity x,y,z:         %12.4f %12.4f %12.4f\n", v1.x, v1.y, v1.z);
    printf("inertial longitude:    %12.6f\n", atan2(p1.y, p1.x)*RTD);
    printf("earth-fixed longitude: %12.6f\n", atan2(p2.y, p2.x)*RTD);
    printf("earth-fixed position x,y,z:      %12.3f %12.3f %12.3f\n", p2.x, p2.y, p2.z);
    printf("coriolis velocity vector x,y,z:  %12.4f %12.4f %12.4f\n", c1.x, c1.y, c1.z);
    printf("earth-fixed velocity x,y,z:      %12.4f %12.4f %12.4f\n", v3.x, v3.y, v3.z);
#endif

    break;

  case EARTH_TO_INERTIAL:
    cross(&spin, &p1, &c1);				/* coriolis velocity (spin vector x position) */
    add(&v1, &c1, &v2);					/* add coriolis */

    p2.z = p1.z;					/* rotate by gha */
    p2.x = p1.x * cos_gha - p1.y * sin_gha;
    p2.y = p1.y * cos_gha + p1.x * sin_gha;
    v3.z = v2.z;					/* rotate by gha */
    v3.x = v2.x * cos_gha - v2.y * sin_gha;
    v3.y = v2.y * cos_gha + v2.x * sin_gha;

#ifdef DEBUG
    printf("\n*** earth-fixed to inertial ***\n");
    printf("earth-fixed position x,y,z:      %12.3f %12.3f %12.3f\n", p1.x, p1.y, p1.z);
    printf("earth-fixed velocity x,y,z:      %12.4f %12.4f %12.4f\n", v1.x, v1.y, v1.z);
    printf("earth-fixed longitude: %12.6f\n", atan2(p1.y, p1.x)*RTD);
    printf("inertial longitude:    %12.6f\n", atan2(p2.y, p2.x)*RTD);
    printf("inertial position x,y,z:         %12.3f %12.3f %12.3f\n", p2.x, p2.y, p2.z);
    printf("coriolis velocity vector  x,y,z: %12.4f %12.4f %12.4f\n", c1.x, c1.y, c1.z);
    printf("inertial velocity x,y,z:         %12.5f %12.5f %12.5f\n", v3.x, v3.y, v3.z);
#endif
    break;
  default:
    fprintf(stderr, "\nERROR ref_frame: invalid state vector conversion flag: %d\n", cflg);
  }
  sv2->pos = p2;
  sv2->vel = v3;
  return;
}

#define T_STEP  0.1
#define DX    1.0

void sv_calc2(STATE *sv1, STATE *sv2, double delta) {
  /*
     propagate an inertial input state vector sv1 by time delta using integration
     of equations of motion

     sv1		input starting state vector
     sv2		output state vector
     delta	time from starting state vector

     clw 7-Jun-2004
  */
  double u, u0;
  double v[3], r[3], f[3];
  double rm, rs, s1, s2, s3, s4, uj2, uj3, uj4, rdotv, as, ac;
  double dt = T_STEP;
  double dx;
  int ns, i, j;

  if (delta < 0)dt = -T_STEP;		/* step backwards in time if delta < 0 */
  ns = nint(delta / dt);
  as = WGS84_MAJOR * WGS84_MAJOR;
  ac = as * WGS84_MAJOR;
  dx = DX;

  r[0] = sv1->pos.x;
  r[1] = sv1->pos.y;
  r[2] = sv1->pos.z;
  v[0] = sv1->vel.x;
  v[1] = sv1->vel.y;
  v[2] = sv1->vel.z;

  for (j = 0; j < ns; j++) {
    rs = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    rm = sqrt(rs);
    s1 = r[2] / rm;	 /* sin(latitude) */
    s2 = r[2] * r[2] / rs;	 /* sin^2(latitude) */
    s3 = s2 * r[2] / rm; /* sin^3(latitude) */
    s4 = s2 * s2;        /* sin^4(latitude) */
    uj2 = as * J2 * (3.*s2 - 1.) / (2.*rs);
    uj3 = ac * J3 * 5.*(s3 - (3.*s1) / 5.) / (2.*rs * rm);
    uj4 = as * as * J4 * 35. * (s4 - (6.*s2) / 7. + 3. / 35.) / (8.*rs * rs);
    u0 = -G_EARTH / rm * (-1.0 + uj2 + uj3 + uj4);/* gravitational pototential */

#ifdef DEBUG
    if (j == 0)printf("u0: %f  uj2: %e  uj3: %e  uj4: %e\n", u0, uj2, uj3, uj4);
#endif

    for (i = 0; i < 3; i++) {			/* calculate acceleration vector grad(u) */
      r[i] += DX;
      rs = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
      rm = sqrt(rs);
      s1 = r[2] / rm;	 	/* sin(latitude) */
      s2 = r[2] * r[2] / rs;	/* sin^2(latitude) */
      s3 = s2 * r[2] / rm;		/* sin^3(latitude) */
      s4 = s2 * s2;		/* sin^4(latitude) */
      uj2 = as * J2 * (3.*s2 - 1.) / (2.*rs);
      uj3 = ac * J3 * 5.*(s3 - (3.*s1) / 5.) / (2.*rs * rm);
      uj4 = as * as * J4 * 35. * (s4 - (6.*s2) / 7. + 3. / 35.) / (8.*rs * rs);
      u = -G_EARTH / rm * (-1.0 + uj2 + uj3 + uj4);
      f[i] = (u - u0) / DX;			/* acceleration */
      r[i] -= DX;
    }
    rdotv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

    for (i = 0; i < 3; i++) {
      r[i] = r[i] + v[i] * dt + f[i] * dt * dt / 2.;
      v[i] = v[i] + f[i] * dt - (G_EARTH * v[i] / (rs * rm) + 3.*f[i] * rdotv / rs) * dt * dt / 2.;
    }
  }

  dt = delta - ns * dt;
#ifdef DEBUG
  printf("number of steps: %d\n", ns);
  printf("dt: %f\n", dt);
#endif
  if (dt != 0.0) {				/* check if there is a small extra increment */
    rs = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    rm = sqrt(rs);
    s1 = r[2] / rm;
    s2 = r[2] * r[2] / rs;
    s3 = s2 * r[2] / rm;
    s4 = s2 * s2;
    uj2 = as * J2 * (3.*s2 - 1.) / (2.*rs);
    uj3 = ac * J3 * 5.*(s3 - (3.*s1) / 5.) / (2.*rs * rm);
    uj4 = as * as * J4 * 35. * (s4 - (6.*s2) / 7. + 3. / 35.) / (8.*rs * rs);
    u0 = -G_EARTH / rm * (-1.0 + uj2 + uj3 + uj4);

    for (i = 0; i < 3; i++) {
      r[i] += DX;
      rs = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
      rm = sqrt(rs);
      s1 = r[2] / rm;
      s2 = r[2] * r[2] / rs;
      s3 = s2 * r[2] / rm;
      s4 = s2 * s2;
      uj2 = as * J2 * (3.*s2 - 1.) / (2.*rs);
      uj3 = ac * J3 * 5.*(s3 - (3.*s1) / 5.) / (2.*rs * rm);
      uj4 = as * as * J4 * 35. * (s4 - (6.*s2) / 7. + 3. / 35.) / (8.*rs * rs);
      u = -G_EARTH / rm * (-1.0 + uj2 + uj3 + uj4);
      f[i] = (u - u0) / DX;
      r[i] -= DX;
    }

    rdotv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2];

    for (i = 0; i < 3; i++) {
      r[i] = r[i] + v[i] * dt + f[i] * dt * dt / 2.;
      v[i] = v[i] + f[i] * dt - (G_EARTH * v[i] / (rs * rm) + 3.*f[i] * rdotv / rs) * dt * dt / 2.;
    }
  }

  sv2->pos.x = r[0];
  sv2->pos.y = r[1];
  sv2->pos.z = r[2];
  sv2->vel.x = v[0];
  sv2->vel.y = v[1];
  sv2->vel.z = v[2];
  return;
}

void sv2kep(VEC *r, VEC *v, KOE *koe) {
  /*
    convert orbital state vector to Keplerian orbital elements for an earth orbiting satellite
    clw 23-Jun-2004
  */
  VEC h;	/* angular momentum */
  VEC e;	/* eccentricity vector */
  VEC r1;	/* unit vector in the range direction */
  VEC t1, t2;
  VEC k1;	/* unit vector in the z direction */
  VEC node;	/* node vector */

  double h2, hm, nodem;

  k1.x = 0.0;
  k1.y = 0.0;
  k1.z = 1.0;

  unit(r, &r1);			/* unit vector in the range direction */
  cross(r, v, &h);		/* angular momentum vector h */
  hm = norm(&h);		/* magnitude of angular momemtum */
  h2 = hm * hm;			/* magnitude squared of h */

  cross(v, &h, &t1);		/* eccentricity vector e */
  smul(1. / G_EARTH, &t1, &t2);
  sub(&t2, &r1, &e);
  koe->e = norm(&e);

  koe->a = h2 / (G_EARTH * (1. - (koe->e * koe->e)));/* semi-major axis */
  koe->inc = acos(h.z / hm);

  cross(&k1, &h, &node);		/* node vector from earth center to ascending node */
  nodem = norm(&node);
  koe->rasc = acos(node.x / nodem);	/* right ascension of the ascending node */
  if (node.y < 0.0)koe->rasc = TWO_PI - koe->rasc;

  koe->omega = acos(dot(&node, &e) / (nodem * koe->e));	/* argument of perigee */
  if (e.z < 0.0) koe->omega = TWO_PI - koe->omega;

  koe->theta = acos(dot(&e, &r1) / koe->e);		/* true anomoly */
  if (dot(r, v) < 0.0) koe->theta = TWO_PI - koe->theta;

  koe->E = acos((koe->e + cos(koe->theta)) / (1.0 + koe->e * cos(koe->theta)));	/* eccentric anomoly */
  if ((koe->theta > PI) && (koe->theta < TWO_PI)) koe->E = TWO_PI - koe->E;

  koe->M = koe->E - koe->e * sin(koe->E);		/* mean anomoly */

#ifdef DEBUG
  {
    double n;		/* orbit angular frequency (rad/s) */
    double tau;		/* time since perigee (s) */
    double vperi, vapo, vcirc, tper;

    vperi = sqrt(G_EARTH / koe->a) * sqrt((1 + koe->e) / (1 - koe->e));
    vapo = sqrt(G_EARTH / koe->a) * sqrt((1 - koe->e) / (1 + koe->e));
    vcirc = sqrt(G_EARTH / koe->a);

    n = sqrt(G_EARTH) / pow(koe->a, 1.5);	/* orbit angular frequency (rad/s) */
    tper = TWO_PI / n; 			/* orbit period (s) */
    tau = koe->M / n;			/* time since perigee */

    printf("\n*** Cartesian to Keplerian Elements ***\n");
    printf("r (x,y,z):  %15.5lf  %15.5lf  %15.5lf  r:   %15.5lf\n", r->x, r->y, r->z, norm(r));
    printf("v (x,y,z):  %15.5lf  %15.5lf  %15.5lf  vel: %15.5lf\n", v->x, v->y, v->z, norm(v));
    /*
        printf("velocity (m/s) perigee:%15.5lf   apogee:%15.5lf   circular:%15.5lf\n",vperi,vapo,vcirc);
    */
    printf("semi-major axis (m):        %15.4lf\n", koe->a);
    printf("eccentricity:               %15.8lf\n", koe->e);
    printf("inclination (deg.):         %15.8lf     deg: %15.8lf\n", koe->inc, koe->inc*RTD);
    printf("right asc. asc.node (deg.): %15.8lf     deg: %15.8lf\n", koe->rasc, koe->rasc*RTD);
    printf("argument of perigee (deg.): %15.8lf     deg: %15.8lf\n", koe->omega, koe->omega*RTD);
    printf("mean anomoly (rad.):        %15.8lf     deg: %15.8lf\n", koe->M, koe->M*RTD);
    printf("eccentric anomoly (rad.):   %15.8lf     deg: %15.8lf\n", koe->E, koe->E*RTD);
    printf("true anomoly (rad.):        %15.8lf     deg: %15.8lf\n", koe->theta, koe->theta*RTD);
    printf("orbit period (s):           %15.8lf\n", tper);
    printf("time since perigee (s):     %15.8lf\n\n", tau);
  }
#endif
}

#define K_EPS  1.0e-15
void kep2sv(VEC *r, VEC *v, KOE *koe) {
  /*
    convert Keplerian orbital elements to state vectors for an earth orbiting satellite
    calculation proceeds with the assumption that the mean anomoly is the input.

    To propagate an orbit state vector using koe:

    KOE koe, koe1;				#allocate KOE data structures
    VEC r,v,r1,v1;				#allocate vectors r=position, v=velocity;

    dt = 30.0;					#propagate 30 s.
    sv2kep(&r, &v, &koe);				#generate keplerian elements koe
    koe1 = koe;					#copy orbital elements to koe1
    n = sqrt(G_EARTH/(koe.a * koe.a * koe.a));	#calculate orbit angular frequency (rad/s)
    koe1.M = koe.M + n*dt;			#update mean anomoly
    koe1.M = fmod(koe1.M, TWO_PI);		#make sure within range 0->2PI
    kep2sv(&r1, &v1, &koe1);			#convert back to get propagated state vector

    clw 23-Jun-2004
  */

  double e1, e2;
  double cap, sap;	/* cos and sin argument of perigee  (omega) */
  double ci, si;	/* cos and sin of inclination angle (i) */
  double cr, sr;	/* cos and sin of right ascension of the ascending node  (Omega) */
  double a1, a2, a3, a4;
  double n;		/* radians/sec orbit angular frequency */
  double Edot;		/* rate of change of eccentric anomoly */

  VEC p, q;
  VEC p1, q1;

  e1 = koe->M;		/* initial estimate of eccentric anomoly */

  while (1) {		/* Newton's method iteration */
    e2 = e1 - ((e1 - koe->e * sin(e1) - koe->M) / (1. - koe->e * cos(e1)));
    if (fabs(e2 - e1) < K_EPS) break;
    e1 = e2;
  }
  koe->E = e2; 		/* save solution */
  /* construct p and q vectors defining orbital plane */
  cap = cos(koe->omega);
  sap = sin(koe->omega);	/* omega: argument of perigee */
  ci = cos(koe->inc);
  si = sin(koe->inc);	/* inc: inclination */
  cr = cos(koe->rasc);
  sr = sin(koe->rasc);	/* rasc: right ascension of the ascending node */

  p.x = cap * cr - sap * ci * sr;
  p.y = cap * sr + sap * ci * cr ;
  p.z = sap * si;

  q.x = -sap * cr - cap * ci * sr;
  q.y = -sap * sr + cap * ci * cr;
  q.z = cap * si;

  a1 = koe->a * (cos(koe->E) - koe->e);		/* calculate position vector r */
  a2 = koe->a * sqrt(1.0 - koe->e * koe->e) * sin(koe->E);

  smul(a1, &p, &p1);
  smul(a2, &q, &q1);
  add(&p1, &q1, r);

  n = sqrt(G_EARTH / (koe->a * koe->a * koe->a));	/* angular frequency of orbit */
  Edot = n / (1.0 - koe->e * cos(koe->E));

  a3 = -koe->a * sin(koe->E) * Edot;
  a4 = koe->a * sqrt(1.0 - koe->e * koe->e) * cos(koe->E) * Edot;

  smul(a3, &p, &p1);
  smul(a4, &q, &q1);
  add(&p1, &q1, v);

  koe->theta = acos((cos(koe->E) - koe->e) / (1.0 - koe->e * cos(koe->E)));
  if ((koe->E > PI) && (koe->E < TWO_PI)) koe->theta = TWO_PI - koe->theta;

#ifdef DEBUG
  {
    printf("\n*** Keplerian Elements to Cartesian ***\n");
    printf("r (x,y,z):  %15.5lf  %15.5lf  %15.5lf  r:   %15.5lf\n", r->x, r->y, r->z, norm(r));
    printf("v (x,y,z):  %15.5lf  %15.5lf  %15.5lf  vel: %15.5lf\n", v->x, v->y, v->z, norm(v));
    printf("semi-major axis (m):        %15.4lf\n", koe->a);
    printf("eccentricity:               %15.8lf\n", koe->e);
    printf("inclination:                %15.8lf     deg: %15.8lf\n", koe->inc, koe->inc*RTD);
    printf("right asc. asc.node (deg.): %15.8lf     deg: %15.8lf\n", koe->rasc, koe->rasc*RTD);
    printf("argument of perigee (deg.): %15.8lf     deg: %15.8lf\n", koe->omega, koe->omega*RTD);
    printf("mean anomoly (rad.):        %15.8lf     deg: %15.8lf\n", koe->M, koe->M*RTD);
    printf("eccentric anomoly (rad.):   %15.8lf     deg: %15.8lf\n", koe->E, koe->E*RTD);
    printf("true anomoly (rad.):        %15.8lf     deg: %15.8lf\n\n", koe->theta, koe->theta*RTD);
  }
#endif
}

double radius(double ra, double e2, double hdg, double lat) {
  /*
      subroutine to calculate the radius of curvature for a particular
      heading (0 degrees = North, +90 = East, -90 = West)

      ra		semi-major axis of the earth for the reference ellipsoid (m)
      e2		square of the excentricity ra**2 - rb**2/ra**2
      hdg   	heading (degrees)
      lat		latitude (degrees)

      clw 23-Jan-97
  */

  double s_lat;
  double rad_n, rad_e, rad_hdg;

  hdg *= DTR;
  lat *= DTR;		/* convert to radians */
  s_lat = sin(lat);

  rad_n = ra * (1. - e2) / pow((1. - e2 * SQR(s_lat)), 1.5);
  rad_e = ra / sqrt(1.0 - e2 * SQR(s_lat));
  rad_hdg = (rad_n * rad_e) / (rad_e * SQR(cos(hdg)) + rad_n * SQR(sin(hdg)));
  printf("north radius (m): %12.3f  east radius (m):    %12.3f\n", rad_n, rad_e);
  printf("heading (deg.):   %12.3f  heading radius (m): %12.3f\n", RTD*hdg, rad_hdg);
  return rad_hdg;
}

void geo_loc(double ra, double e2, double rho, double fd, double fc, double alt, int s_flag, VEC *s, VEC *v, double *lat, double *lon) {
  /*
  	subroutine to calculate the lat/lon/altitude of the intersection of the
  	radar look vector with the geoid for a specfied average terrain height.

          ra	ellipsoid semi-major axis
          e2	ellipsoid square of eccentricity
      	rho 	slant range (m)
  	fd	Doppler frequency (Hz)
  	fc	radar carrier frequency (Hz)
  	alt	nominal terrain altitude (m)
          s_flag  flag to determine if right or left looking (left looking=1, right looking= -1)
  	s	position vector of the radar, (earth fixed coordinates)(m)
  	v	velocity vector of the radar (earth fixed coordinates)(m/s)
  	lat	output latitude (decimal degrees)
  	lon	output longitude(decimal degrees)

  */
  double a, b, c, r, r1, lam;
  double s2, v2, t2, sv, det;
  double lat1, lon1, alt1, c1, c2, c3;
  double r_new;
  VEC q, t;
  int iter;

  lam = C / fc;
  s2 = dot(s, s);
  v2 = dot(v, v);
  sv = dot(v, s);
  xyz_ll(ra, e2, s, &lat1, &lon1, &alt1);
  r = sqrt(s2) - alt1;
  cross(s, v, &t);
  t2 = dot(&t, &t);

#ifdef DEBUG_GEO
  printf("T vector:        %13.6e %13.6e %13.6e\n", t.x, t.y, t.z);
  printf("position vector: %13.6e %13.6e %13.6e\n", s->x, s->y, s->z);
  printf("velocity vector: %13.6e %13.6e %13.6e\n", v->x, v->y, v->z);
  printf("SAR lat,lon,alt: %12.5f %12.5f %12.3f\n", lat1, lon1, alt1);
  printf("initial estimate of earth radius: %12.3f\n", r);
#endif

  iter = 0;
  while (iter < 8 ) {
    iter++;
    a = (s2 + SQR(r + alt) - SQR(rho)) / 2.0;
    b = (lam * rho * fd + 2.0 * sv) / 2.0;
    c = SQR(r + alt);
#ifdef DEBUG_GEO
    printf("a,b,c: %13.6e %13.6e %13.6e\n", a, b, c);
#endif

    det = s2 * v2 - SQR(sv);

    c1 = (a * v2 - b * sv) / det;
    c2 = (b * s2 - a * sv) / det;
    c3 = (double)s_flag * sqrt((c - SQR(c1) * s2 - SQR(c2) * v2 - 2.0 * c1 * c2 * sv) / t2);
#ifdef DEBUG_GEO
    printf("\niteration: %2d c1,c2,c3:        %13.6e %13.6e %13.6e\n", iter, c1, c2, c3);
#endif

    /*  radius vector to image point is q = c1*s + c2*v + c3*t */

    q.x = c1 * s->x + c2 * v->x + c3 * t.x;
    q.y = c1 * s->y + c2 * v->y + c3 * t.y;
    q.z = c1 * s->z + c2 * v->z + c3 * t.z;
#ifdef DEBUG_GEO
    printf("iteration: %2d position vector: %13.6e %13.6e %13.6e\n", iter, q.x, q.y, q.z);
#endif
    xyz_ll(ra, e2, &q, &lat1, &lon1, &alt1);
#ifdef DEBUG_GEO
    printf("latitude:  %12.6f  longitude: %12.6f   altitude: %12.3f\n", lat1, lon1, alt1);
#endif

    r1 = norm(&q);			/* includes surface altitude of alt1, earth radius=r1-alt1 */
    r_new = r1 - alt1;
    if (fabs(r_new - r) < .001)break;
    r = r_new;
  }

  *lat = lat1;				/* return solution */
  *lon = lon1;

#ifdef DEBUG_GEO
  printf("re*s:  %13.6e\n", dot(&q, s));
  printf("re*v:  %13.6e\n", dot(&q, v));
  printf("re*re: %13.6e\n", dot(&q, &q));
#endif
}

void ll_xyz(double ra, double e2, double lat, double lon, double alt, VEC *s) {
  /*
      convert from latitude, longitude and altitude above the reference geoid to
      a cartesian vector.

      ra    input semi major axis
      e2    input squared eccentricity
      lat	  latitude (decimal degrees)
      lon	  longitude (decimal degrees)
      alt	  altitude (meters)
      s	  radial cartesian vector (x,y,z meters)

    clw 23-Jan-97
  */
  static double nu;

  lat *= DTR;
  lon *= DTR;		/* convert to radians */
  nu = ra / sqrt(1. - e2 * SQR(sin(lat)));	/* curvature prime vertical (east-west) */
  s->x = (nu + alt) * cos(lon) * cos(lat);
  s->y = (nu + alt) * sin(lon) * cos(lat);
  s->z = (nu * (1.0 - e2) + alt) * sin(lat);
}

void xyz_ll(double ra, double e2, VEC *s, double *lat, double *lon, double *alt) {
  /*
     convert from a catesian vector to latitude, longitude and altitude above
     the reference geoid

      ra    input semi major axis
      e2    input squared eccentricity
      s	  input radial cartesian vector (x,y,z meters)
      lat	  output latitude (decimal degrees) (pointer to the latitude)
      lon	  output longitude (decimal degrees) (pointer to the longitude)
      alt	  output altitude (meters) (pointer to the altitude)

      clw 23-Mar-97
  */
  static double e1, p, nu;
  static double phi, theta, lam;
  static double rb;
  static double st, ct;

  rb = ra * sqrt(1. - e2);
  e1 = e2 * SQR(ra) / SQR(rb);

  p = sqrt(SQR(s->x) + SQR(s->y));
  theta = atan((s->z * ra) / (p * rb));
#ifdef DEBUG_GEO
  printf("p:  %13.6e\n", p);
  printf("rb:  %13.6e\n", rb);
  printf("theta: %13.6e\n", theta);
#endif

  st = sin(theta);
  ct = cos(theta);
  phi = atan((s->z + e1 * rb * st * st * st) / (p - e2 * ra * ct * ct * ct));
  /*phi = atan2((s->z +e1*rb*st*st*st),(p - e2*ra * ct*ct*ct)); */
  lam = atan2(s->y, s->x);

  *lat = phi * RTD;
  *lon = lam * RTD;
  nu = ra / sqrt(1. - e2 * SQR(sin(phi)));	/* radius of curvature in the prime vertical (east-west) */
  *alt = p / cos(phi) - nu;
}

void euler(double yaw, double pitch, double roll, MAT e) {
  /*
    Calculate Euler matrix from yaw, pitch, roll, ,nput angles are in radians.

    The first rotation is is the yaw angle about the z-axis, the second is the pitch angle
    about an intermediary y-axis, and the third is a bank or roll angle about the final x-axis.
    This is the so called XYZ or 321 sequence convention for rotations.

    reprojection of a vector in the XYZ basis into the new basis (TCN into body frame (IJK) coordinates)
    The transpose of the matrix can be used to transform vectors (e.g. IMU lever arms)
    expressed in IJK body coordinates into the TCN reference frame.

    If this routine is used to calculate a transformation matrix (e.g. calculate the look vector given
    a yaw,pitch and roll of the T vector), then the angles must be negated since rotation of
    the vector is the opposite of rotating the reference basis set.

    Reference: Classical Mechanics, 2nd Ed., author: H. Goldstein, Appendix B, pp. 606-610.
               Addison-Wesley, 1980.

    16-aug-94 clw

  */
  double cy, sy;		/* cosine and sin of yaw angle */
  double cp, sp;		/* cosine and sin of pitch angle */
  double cr, sr;		/* cosine and sin of roll angle */


  cy = cos(yaw);
  sy = sin(yaw);
  cp = cos(pitch);
  sp = sin(pitch);
  cr = cos(roll);
  sr = sin(roll);

  e[0][0] = cp * cy;
  e[0][1] = cp * sy;
  e[0][2] = (-sp);

  e[1][0] = sr * sp * cy - cr * sy;
  e[1][1] = sr * sp * sy + cr * cy;
  e[1][2] = cp * sr;

  e[2][0] = cr * sp * cy + sr * sy;
  e[2][1] = cr * sp * sy - sr * cy;
  e[2][2] = cp * cr;
}

void c_tcn(VEC *pc, VEC *vc, MAT tcn) {
  /*
     subroutine to construct the TCN orthonormal basis vector set given the
     position vector pc, and velocity vector vc. These two
     vectors should be supplied in the XYZ, or East,North,Up (ENU)
     orthogonal basis sets.

  	T:	Track vector
  	C:	Cross-track vector
  	N:	Normal (down pointing )

  	clw 6-sep-93
  */
  VEC t, c, n, tmp;

  tmp.x = (-pc->x);	/* construct N vector from platform to origin */
  tmp.y = (-pc->y);	/* pointing from platform to the origin */
  tmp.z = (-pc->z);
  unit(&tmp, &n);	/* make unit normal vector N */

  cross(&n, vc, &tmp);	/* construct cross track vector C */
  unit(&tmp, &c); 	/* unit vector C */

  cross(&c, &n, &t);	/* construct unit  T vector */
  vec_mat(&t, &c, &n, tcn); /* construct tcn matrix from vectors */
}

void ned_calc(double lat, double lon, VEC *nv, VEC *ev, VEC *dv) {
  /* Calc. north, east and down vectors given latitude and longitude */

  lat *= DTR;
  lon *= DTR;

  nv->x = -sin(lat) * cos(lon);
  nv->y = -sin(lat) * sin(lon);
  nv->z =  cos(lat);

  ev->x = -sin(lon);
  ev->y =  cos(lon);
  ev->z = 0.0;

  cross(nv, ev, dv);
}

double a2_dbl(char *s1, int offset, int len) {
  /*
    subroutine to convert string extracted from array to double
    29-Aug-98 clw
  */

  static char str[80];
  double v;

  strncpy(str, &s1[offset], len);
  str[len] = '\0';

  if (sscanf(str, "%lf", &v) != 1) {
    fprintf(stderr, "\nERROR: subroutine a2_dbl, cannot convert string to double floating point: %s\n\n", str);
    exit(-2);
  }
  return v;
}

int a2_int(char *s1, int offset, int len) {
  /*
    subroutine to convert string extracted from array to int
    29-Aug-98 clw
  */

  static char str[80];
  int v;

  strncpy(str, &s1[offset], len);
  str[len] = '\0';

  if (sscanf(str, "%d", &v) != 1) {
    fprintf(stderr, "\nERROR: subroutine a2_int, cannot convert string to integer: %s\n\n", str);
    exit(-2);
  }
  return v;
}

int read_int4(char *a, int j)
/* subroutine to extract a 4 byte integer from a byte array.
   uw/clw 2-Feb-2000
*/
{
  unsigned char b[4];
  int *x;
  b[0] = a[j];
  b[1] = a[j+1];
  b[2] = a[j+2];
  b[3] = a[j+3];
  x = (int *) & b[0];
  return(*x);
}

float read_float(char *a, int j)
/* subroutine to extract a 4 byte float from a byte array.
   clw 4-Feb-2003
*/
{
  unsigned char b[4];
  float *x;
  b[0] = a[j];
  b[1] = a[j+1];
  b[2] = a[j+2];
  b[3] = a[j+3];
  x = (float *) & b[0];
  return(*x);
}

double read_double(char *a, int j)
/* subroutine to extract an 8-byte double from a byte array.
   clw 26-May-2004
*/
{
  unsigned char b[8];
  double *x;
  b[0] = a[j];
  b[1] = a[j+1];
  b[2] = a[j+2];
  b[3] = a[j+3];
  b[4] = a[j+4];
  b[5] = a[j+5];
  b[6] = a[j+6];
  b[7] = a[j+7];
  x = (double *) & b[0];
  return(*x);
}

int read_short2(char *a, int j)
/* subroutine to extract a 2 byte short integer from a byte array.
   uw/clw 2-Feb-2000
*/
{
  unsigned char b[2];
  short *x;
  b[0] = a[j];
  b[1] = a[j+1];
  x = (short *) & b[0];
  return(*x);
}

unsigned char *uvector(int nl, int nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
  unsigned char *v;

  v = (unsigned char *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
  if (!v) nrerror("allocation failure in uvector()");
  return v -nl + NR_END;
}

void free_uvector(unsigned char *v, int nl, int nh )
/* free an unsigned char vector allocated with uvector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}

int *ivector(int nl, int nh) {
  /* allocate an int vector with subscript range v[nl..nh] */
  int *v;

  v = (int *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(int)));
  if (!v) nrerror("allocation failure in ivector()");
  return v -nl + NR_END;
}

void free_ivector(int *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}

short *i2vector(int nl, int nh) {
  /* allocate an int vector with subscript range v[nl..nh] */
  short *v;

  v = (short *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(short)));
  if (!v) nrerror("allocation failure in ivector()");
  return v -nl + NR_END;
}

void free_i2vector(short *v, int nl, int nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}


float *vector(int nl, int nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;

  v = (float *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(float)));
  if (!v) nrerror("ERROR: memory allocation failure in vector()");
  return v -nl + NR_END;
}

void free_vector(float *v, int nl, int nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}


double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v = (double *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(double)));
  if (!v) nrerror("allocation failure in dvector()");
  return v -nl + NR_END;
}

void free_dvector(double *v, int nl, int nh )
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}

fcomplex *cvector(int nl, int nh)
/* allocate a fcomplex vector with subscript range v[nl..nh] */
{
  fcomplex *v;

  v = (fcomplex *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(fcomplex)));
  if (!v) nrerror("ERROR: allocation failure in cvector()");
  return v -nl + NR_END;
}

void free_cvector(fcomplex *v, int nl, int nh )
/* free a fcomplex vector allocated with cvector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}

scomplex *scvector(int nl, int nh)
/* allocate a scomplex vector with subscript range v[nl..nh] */
{
  scomplex *v;

  v = (scomplex *)malloc((size_t) ((nh - nl + 1 + NR_END) * sizeof(scomplex)));
  if (!v) nrerror("ERROR: allocation failure in scvector()");
  return v -nl + NR_END;
}

void free_scvector(scomplex *v, int nl, int nh )
/* free a scomplex vector allocated with scvector() */
{
  free((FREE_ARG) (v + nl - NR_END));
}

short **i2matrix(int nrl, int nrh, int ncl, int nch)
/* allocate a short int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  short **m;

  /* allocate pointers to rows */
  m = (short **) malloc((size_t)((nrow + NR_END) * sizeof(short*)));
  if (!m) nrerror("ERROR: allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (short *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(short)));
  if (!m[nrl]) nrerror("ERROR: allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_i2matrix(short **m, int nrl, int nrh, int ncl, int nch)
/* free a short matrix allocated by i2matrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

int **imatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc((size_t)((nrow + NR_END) * sizeof(int*)));
  if (!m) nrerror("ERROR: allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
  if (!m[nrl]) nrerror("ERROR: allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
/* free a int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

float **matrix(int nrl, int nrh, int ncl, int nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  float **m;

  /* allocate pointers to rows */
  m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
  if (!m) nrerror("ERROR: allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (float *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
  if (!m[nrl]) nrerror("ERROR: allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_matrix(float **m, int nrl, int nrh, int ncl, int nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc((size_t)((nrow + NR_END) * sizeof(double*)));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (double *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
  if (!m[nrl]) nrerror("allocation failure 2 in dmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

fcomplex **cmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a fcomplex matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  fcomplex **m;

  /* allocate pointers to rows */
  m = (fcomplex **)malloc((size_t)((nrow + NR_END) * sizeof(fcomplex*)));
  if (!m) nrerror("ERROR: memory allocation failure 1 in cmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (fcomplex *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(fcomplex)));
  if (!m[nrl]) nrerror("ERROR: memory allocation failure 2 in cmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_cmatrix(fcomplex **m, int nrl, int nrh, int ncl, int nch)
/* free a fcomplex matrix allocated by cmatrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

scomplex **scmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a scomplex matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  scomplex **m;

  /* allocate pointers to rows */
  m = (scomplex **)malloc((size_t)((nrow + NR_END) * sizeof(scomplex*)));
  if (!m) nrerror("ERROR: memory allocation failure 1 in scmatrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] = (scomplex *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(scomplex)));
  if (!m[nrl]) nrerror("ERROR: memory allocation failure 2 in scmatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_scmatrix(scomplex **m, int nrl, int nrh, int ncl, int nch)
/* free a scomplex matrix allocated by scmatrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

unsigned char **umatrix(int nrl, int nrh, int ncl, int nch)
/* allocate an unsigned char matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  unsigned char **m;

  /* allocate pointers to rows */
  m = (unsigned char **) malloc((size_t)((nrow + NR_END) * sizeof(unsigned char*)));
  if (!m) nrerror("allocation failure 1 in umatrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl] = (unsigned char *) malloc((size_t)((nrow * ncol + NR_END) * sizeof(unsigned char)));
  if (!m[nrl]) nrerror("allocation failure 2 in umatrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1;i <= nrh;i++) m[i] = m[i-1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

void free_umatrix(unsigned char **m, int nrl, int nrh, int ncl, int nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl] + ncl - NR_END));
  free((FREE_ARG) (m + nrl - NR_END));
}

#define IGREG (15+31L*(10+12L*1582))
int julday( int id, int mm, int iyyy) {
  int jul;
  int ja, jy = iyyy, jm;

  if (jy == 0) {
    fprintf(stderr, "\nERROR in subroutine julday: there is no year zero!\n");
    exit(-1);
  }

  if (jy < 0) ++jy;
  if (mm > 2) {
    jm = mm + 1;
  } else {
    --jy;
    jm = mm + 13;
  }
  jul = (int) (floor(365.25 * jy) + floor(30.6001 * jm) + id + 1720995);
  if (id + 31L*(mm + 12L*iyyy) >= IGREG) {
    ja = (int)(0.01 * jy);
    jul += 2 - ja + (int) (0.25 * ja);
  }
  return jul;
}
#undef IGREG

#define IGREG 2299161

void caldat(int julian, int *id, int *mm, int *iyyy) {
  int ja, jalpha, jb, jc, jd, je;

  if (julian >= IGREG) {
    jalpha = (int)(((float) (julian - 1867216) - 0.25) / 36524.25);
    ja = julian + 1 + jalpha - (int) (0.25 * jalpha);
  } else
    ja = julian;
  jb = ja + 1524;
  jc = (int)(6680.0 + ((float) (jb - 2439870) - 122.1) / 365.25);
  jd = (int)(365 * jc + (0.25 * jc));
  je = (int)((jb - jd) / 30.6001);
  *id = jb - jd - (int) (30.6001 * je);
  *mm = je - 1;
  if (*mm > 12) *mm -= 12;
  *iyyy = jc - 4715;
  if (*mm > 2) --(*iyyy);
  if (*iyyy <= 0) --(*iyyy);
}
#undef IGREG

void polint(double xa[], double ya[], int n, double x, double *y, double *dy) {
  int i, m, ns = 1;
  double den, dif, dift, ho, hp, w;
  double *c, *d;

  dif = fabs(x - xa[1]);
  c = dvector(1, n);
  d = dvector(1, n);
  for (i = 1;i <= n;i++) {
    if ( (dift = fabs(x - xa[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  *y = ya[ns--];
  for (m = 1;m < n;m++) {
    for (i = 1;i <= n - m;i++) {
      ho = xa[i] - x;
      hp = xa[i+m] - x;
      w = c[i+1] - d[i];
      if ( (den = ho - hp) == 0.0) nrerror("ERROR: subroutine polint, infinite slope!");
      den = w / den;
      d[i] = hp * den;
      c[i] = ho * den;
    }
    *y += (*dy = (2 * ns < (n - m) ? c[ns+1] : d[ns--]));
  }
  free_dvector(d, 1, n);
  free_dvector(c, 1, n);
}


void pos_sv(double tpos, STATE *sv, double t_state, double tis, int nstate, VEC *pos) {
  /*
   v1.2 clw 4-Apr-2003 interpolate position state vectors
   v1.3 clw 11-Jun-2008 increased polynomial order to 8

      tpos	time to estimate position (s)
      sv	state vectors (position and velocity vectors)
      tstate	time of first state vector (s)
      tis	time interval between state vectors (s)
      nstate	number of state vectors
      pos	interpolated position vector
  */

  double t[PPTS], px[PPTS], py[PPTS], pz[PPTS];
  double px_err, py_err, pz_err;
  double tmin, tmax;
  int i, ix, i1;
  int ppts;

  ppts = IMIN(PPTS, nstate);
  tmin = t_state - MAX_INTERP_TIME;
  tmax = t_state + (nstate - 1) * tis + MAX_INTERP_TIME;

  if (tpos < tmin) {
    fprintf(stderr, "\nERROR pos_sv: time outside of range: %12.5f  min:%12.5f  max:%12.5f\n", tpos, tmin, tmax);
    exit(-1);
  }

  if (tpos > tmax) {
    fprintf(stderr, "\nERROR pos_sv: time outside of range: %12.5f  min:%12.5f  max:%12.5f\n", tpos, tmin, tmax);
    exit(-1);
  }

  ix = nint((tpos - t_state) / tis - (ppts - 1.0) / 2.0);
  ix = MAX(ix, 0);
  ix = MIN(ix, nstate - ppts);

  for (i = 0; i < ppts; i++) {
    i1 = ix + i;
    px[i] = sv[i1].pos.x;
    py[i] = sv[i1].pos.y;
    pz[i] = sv[i1].pos.z;
    t[i] = t_state + i1 * tis;
  }
  polint(t - 1, px - 1, ppts, tpos, &(pos->x), &px_err);
  polint(t - 1, py - 1, ppts, tpos, &(pos->y), &py_err);
  polint(t - 1, pz - 1, ppts, tpos, &(pos->z), &pz_err);
}

void vel_sv(double tpos, STATE *sv, double t_state, double tis, int nstate, VEC *vel) {
  /*
   v1.2 clw 4-Apr-2003 interpolate velocity state vectors
   v1.3 clw 11-Jun-2008 increased polynomial order to 8

      tpos	time to estimate position (s)
      sv	state vectors (position and velocity vectors)
      tstate	time of first state vector (s)
      tis	time interval between state vectors (s)
      nstate	number of state vectors
      vel	interpolated velocity vector
  */
  double t[PPTS], vx[PPTS], vy[PPTS], vz[PPTS];
  double vx_err, vy_err, vz_err;
  double tmin, tmax;
  int i, ix, i1;
  int ppts;

  ppts = IMIN(PPTS, nstate);
  tmin = t_state - MAX_INTERP_TIME;
  tmax = t_state + (nstate - 1) * tis + MAX_INTERP_TIME;

  if (tpos < tmin) {
    fprintf(stderr, "\nERROR vel_sv: time outside of range: %12.5f  min:%12.5f  max:%12.5f\n", tpos, tmin, tmax);
    exit(-1);
  }

  if (tpos > tmax) {
    fprintf(stderr, "\nERROR vel_sv: time outside of range: %12.5f  min:%12.5f  max:%12.5f\n", tpos, tmin, tmax);
    exit(-1);
  }

  ix = nint((tpos - t_state) / tis - (ppts - 1.0) / 2.0);
  ix = MAX(ix, 0);
  ix = MIN(ix, nstate - ppts);

  for (i = 0; i < ppts; i++) {
    i1 = ix + i;
    vx[i] = sv[i1].vel.x;
    vy[i] = sv[i1].vel.y;
    vz[i] = sv[i1].vel.z;
    t[i] = t_state + i1 * tis;
  }
  polint(t - 1, vx - 1, ppts, tpos, &(vel->x), &vx_err);
  polint(t - 1, vy - 1, ppts, tpos, &(vel->y), &vy_err);
  polint(t - 1, vz - 1, ppts, tpos, &(vel->z), &vz_err);
}

double bessi0(double x) {
  double ax, ans;
  double y;

  if ((ax = fabs(x)) < 3.75) {
    y = x / 3.75;
    y *= y;
    ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                      + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
  } else {
    y = 3.75 / ax;
    ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                                  + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                                                            + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                                                                                        + y * 0.392377e-2))))))));
  }
  return ans;
}

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void four1(float data[], unsigned int nn, int isign) {
  unsigned int n, mmax, m, j, istep, i;
  double wtemp, wr, wpr, wpi, wi, theta;
  float tempr, tempi, temp;

  n = nn << 1;
  j = 1;
  for (i = 1;i < n;i += 2) {
    if (j > i) {
      SWAP(data[j], data[i]);
      SWAP(data[j+1], data[i+1]);
    }
    m = n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax = 2;
  while (n > mmax) {
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1;m < mmax;m += 2) {
      for (i = m;i <= n;i += istep) {
        j = i + mmax;
        tempr = wr * data[j] - wi * data[j+1];
        tempi = wr * data[j+1] + wi * data[j];
        data[j] = data[i] - tempr;
        data[j+1] = data[i+1] - tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr = (wtemp = wr) * wpr - wi * wpi + wr;
      wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
  }
}
#undef SWAP

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
void fourn(float data[], unsigned int nn[], int ndim, int isign) {
  int idim;
  unsigned int i1, i2, i3, i2rev, i3rev, ip1, ip2, ip3, ifp1, ifp2;
  unsigned int ibit, k1, k2, n, nprev, nrem, ntot;
  float tempi, tempr;
  double theta, wi, wpi, wpr, wr, wtemp;

  for (ntot = 1, idim = 1;idim <= ndim;idim++)
    ntot *= nn[idim];
  nprev = 1;
  for (idim = ndim;idim >= 1;idim--) {
    n = nn[idim];
    nrem = ntot / (n * nprev);
    ip1 = nprev << 1;
    ip2 = ip1 * n;
    ip3 = ip2 * nrem;
    i2rev = 1;
    for (i2 = 1;i2 <= ip2;i2 += ip1) {
      if (i2 < i2rev) {
        for (i1 = i2;i1 <= i2 + ip1 - 2;i1 += 2) {
          for (i3 = i1;i3 <= ip3;i3 += ip2) {
            i3rev = i2rev + i3 - i2;
            SWAP(data[i3], data[i3rev]);
            SWAP(data[i3+1], data[i3rev+1]);
          }
        }
      }
      ibit = ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
        i2rev -= ibit;
        ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1 = ip1;
    while (ifp1 < ip2) {
      ifp2 = ifp1 << 1;
      theta = isign * 6.28318530717959 / (ifp2 / ip1);
      wtemp = sin(0.5 * theta);
      wpr = -2.0 * wtemp * wtemp;
      wpi = sin(theta);
      wr = 1.0;
      wi = 0.0;
      for (i3 = 1;i3 <= ifp1;i3 += ip1) {
        for (i1 = i3;i1 <= i3 + ip1 - 2;i1 += 2) {
          for (i2 = i1;i2 <= ip3;i2 += ifp2) {
            k1 = i2;
            k2 = k1 + ifp1;
            tempr = (float)wr * data[k2] - (float)wi * data[k2+1];
            tempi = (float)wr * data[k2+1] + (float)wi * data[k2];
            data[k2] = data[k1] - tempr;
            data[k2+1] = data[k1+1] - tempi;
            data[k1] += tempr;
            data[k1+1] += tempi;
          }
        }
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
      }
      ifp1 = ifp2;
    }
    nprev *= n;
  }
}
#undef SWAP

void realft(float data[], unsigned int n, int isign) {
  void four1(float data[], unsigned int nn, int isign);
  unsigned int i, i1, i2, i3, i4, np3;
  float c1 = 0.5, c2, h1r, h1i, h2r, h2i;
  double wr, wi, wpr, wpi, wtemp, theta;

  theta = 3.141592653589793 / (double) (n >> 1);
  if (isign == 1) {
    c2 = -0.5;
    four1(data, n >> 1, 1);
  } else {
    c2 = 0.5;
    theta = -theta;
  }
  wtemp = sin(0.5 * theta);
  wpr = -2.0 * wtemp * wtemp;
  wpi = sin(theta);
  wr = 1.0 + wpr;
  wi = wpi;
  np3 = n + 3;
  for (i = 2;i <= (n >> 2);i++) {
    i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
    h1r = c1 * (data[i1] + data[i3]);
    h1i = c1 * (data[i2] - data[i4]);
    h2r = -c2 * (data[i2] + data[i4]);
    h2i = c2 * (data[i1] - data[i3]);
    data[i1] = h1r + wr * h2r - wi * h2i;
    data[i2] = h1i + wr * h2i + wi * h2r;
    data[i3] = h1r - wr * h2r + wi * h2i;
    data[i4] = -h1i + wr * h2i + wi * h2r;
    wr = (wtemp = wr) * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
  }
  if (isign == 1) {
    data[1] = (h1r = data[1]) + data[2];
    data[2] = h1r - data[2];
  } else {
    data[1] = c1 * ((h1r = data[1]) + data[2]);
    data[2] = c1 * (h1r - data[2]);
    four1(data, n >> 1, -1);
  }
}

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void lfit1d(double x1[], double y[], double sig[], int ndat, double a[], int ia[],
            int ma, double **covar, double *chisq, void (*funcs)(double, double *)) {
  void covsrt(double **covar, int ma, int ia[], int mfit);
  void gaussj(double **a, int n, double **b, int m);

  int i, j, k, l, m, mfit = 0;
  double ym, wt, sum, sig2i, **beta, *afunc;

  beta = dmatrix(1, ma, 1, 1);
  afunc = dvector(1, ma);
  for (j = 1;j <= ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) nrerror("ERROR lfit: no parameters to be fitted");
  for (j = 1;j <= mfit;j++) {
    for (k = 1;k <= mfit;k++) covar[j][k] = 0.0;
    beta[j][1] = 0.0;
  }
  for (i = 1;i <= ndat;i++) {
    (*funcs)(x1[i], afunc);
    ym = y[i];
    if (mfit < ma) {
      for (j = 1;j <= ma;j++)
        if (!ia[j]) ym -= a[j] * afunc[j];
    }
    sig2i = 1.0 / SQR(sig[i]);
    for (j = 0, l = 1;l <= ma;l++) {
      if (ia[l]) {
        wt = afunc[l] * sig2i;
        for (j++, k = 0, m = 1;m <= l;m++)
          if (ia[m]) covar[j][++k] += wt * afunc[m];
        beta[j][1] += ym * wt;
      }
    }
  }
  for (j = 2;j <= mfit;j++)
    for (k = 1;k < j;k++)
      covar[k][j] = covar[j][k];
  gaussj(covar, mfit, beta, 1);
  for (j = 0, l = 1;l <= ma;l++)
    if (ia[l]) a[l] = beta[++j][1];
  *chisq = 0.0;
  for (i = 1;i <= ndat;i++) {
    (*funcs)(x1[i], afunc);
    for (sum = 0.0, j = 1;j <= ma;j++) sum += a[j] * afunc[j];
    *chisq += SQR((y[i] - sum) / sig[i]);
  }
  covsrt(covar, ma, ia, mfit);
  free_dvector(afunc, 1, ma);
  free_dmatrix(beta, 1, ma, 1, 1);
}
#undef SWAP

#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
void lfit(double x1[], double x2[], double y[], double sig[], int ndat, double a[], int ia[],
          int ma, double **covar, double *chisq, void (*funcs)(double, double, double *)) {
  void covsrt(double **covar, int ma, int ia[], int mfit);
  void gaussj(double **a, int n, double **b, int m);

  int i, j, k, l, m, mfit = 0;
  double ym, wt, sum, sig2i, **beta, *afunc;

  beta = dmatrix(1, ma, 1, 1);
  afunc = dvector(1, ma);
  for (j = 1;j <= ma;j++)
    if (ia[j]) mfit++;
  if (mfit == 0) nrerror("ERROR lfit: no parameters to be fitted");
  for (j = 1;j <= mfit;j++) {
    for (k = 1;k <= mfit;k++) covar[j][k] = 0.0;
    beta[j][1] = 0.0;
  }
  for (i = 1;i <= ndat;i++) {
    (*funcs)(x1[i], x2[i], afunc);
    ym = y[i];
    if (mfit < ma) {
      for (j = 1;j <= ma;j++)
        if (!ia[j]) ym -= a[j] * afunc[j];
    }
    sig2i = 1.0 / SQR(sig[i]);
    for (j = 0, l = 1;l <= ma;l++) {
      if (ia[l]) {
        wt = afunc[l] * sig2i;
        for (j++, k = 0, m = 1;m <= l;m++)
          if (ia[m]) covar[j][++k] += wt * afunc[m];
        beta[j][1] += ym * wt;
      }
    }
  }
  for (j = 2;j <= mfit;j++)
    for (k = 1;k < j;k++)
      covar[k][j] = covar[j][k];
  gaussj(covar, mfit, beta, 1);

  for (j = 0, l = 1;l <= ma;l++)
    if (ia[l]) a[l] = beta[++j][1];
  *chisq = 0.0;
  for (i = 1;i <= ndat;i++) {
    (*funcs)(x1[i], x2[i], afunc);
    for (sum = 0.0, j = 1;j <= ma;j++) sum += a[j] * afunc[j];
    *chisq += SQR((y[i] - sum) / sig[i]);
  }
  covsrt(covar, ma, ia, mfit);
  free_dvector(afunc, 1, ma);
  free_dmatrix(beta, 1, ma, 1, 1);
}
#undef SWAP
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void gaussj(double **a, int n, double **b, int m) {
  int *indxc, *indxr, *ipiv;
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, temp;

  indxc = ivector(1, n);
  indxr = ivector(1, n);
  ipiv = ivector(1, n);
  for (j = 1;j <= n;j++) ipiv[j] = 0;
  for (i = 1;i <= n;i++) {
    big = 0.0;
    for (j = 1;j <= n;j++)
      if (ipiv[j] != 1)
        for (k = 1;k <= n;k++) {
          if (ipiv[k] == 0) {
            if (fabs(a[j][k]) >= big) {
              big = fabs(a[j][k]);
              irow = j;
              icol = k;
            }
          } else if (ipiv[k] > 1) {
            fprintf(stderr, "\nERROR: gaussj: singular matrix 1\n");
            exit(-1);
          }
        }
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l = 1;l <= n;l++) SWAP(a[irow][l], a[icol][l])
        for (l = 1;l <= m;l++) SWAP(b[irow][l], b[icol][l])
        }
    indxr[i] = irow;
    indxc[i] = icol;
    if (a[icol][icol] == 0.0) {
      fprintf(stderr, "\nERROR: gaussj: singular matrix 2\n");
      exit(-1);
    }
    pivinv = 1.0 / a[icol][icol];
    a[icol][icol] = 1.0;
    for (l = 1;l <= n;l++) a[icol][l] *= pivinv;
    for (l = 1;l <= m;l++) b[icol][l] *= pivinv;
    for (ll = 1;ll <= n;ll++)
      if (ll != icol) {
        dum = a[ll][icol];
        a[ll][icol] = 0.0;
        for (l = 1;l <= n;l++) a[ll][l] -= a[icol][l] * dum;
        for (l = 1;l <= m;l++) b[ll][l] -= b[icol][l] * dum;
      }
  }
  for (l = n;l >= 1;l--) {
    if (indxr[l] != indxc[l])
      for (k = 1;k <= n;k++)
        SWAP(a[k][indxr[l]], a[k][indxc[l]]);
  }
  free_ivector(ipiv, 1, n);
  free_ivector(indxr, 1, n);
  free_ivector(indxc, 1, n);
}
#undef SWAP

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
void covsrt(double **covar, int ma, int ia[], int mfit) {
  int i, j, k;
  double swap;

  for (i = mfit + 1;i <= ma;i++)
    for (j = 1;j <= i;j++) covar[i][j] = covar[j][i] = 0.0;
  k = mfit;
  for (j = ma;j >= 1;j--) {
    if (ia[j]) {
      for (i = 1;i <= ma;i++) SWAP(covar[i][k], covar[i][j])
        for (i = 1;i <= ma;i++) SWAP(covar[k][i], covar[j][i])
          k--;
    }
  }
}
#undef SWAP
#define TOL 1.0e-20
void svdfit(double x1[], double x2[], double y[], double sig[], int ndata, double a[], int ma,
            double **u, double **v, double w[], double *chisq, void (*funcs)(double, double, double *)) {
  void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]);
  void svdcmp(double **a, int m, int n, double w[], double **v);

  int j, i;
  double wmax, tmp, thresh, sum, *b, *afunc;

  b = dvector(1, ndata);
  afunc = dvector(1, ma);
  for (i = 1;i <= ndata;i++) {
    (*funcs)(x1[i], x2[i], afunc);
    tmp = 1.0 / sig[i];
    for (j = 1;j <= ma;j++) u[i][j] = afunc[j] * tmp;
    b[i] = y[i] * tmp;
  }
  svdcmp(u, ndata, ma, w, v);
  wmax = 0.0;
  for (j = 1;j <= ma;j++)
    if (w[j] > wmax) wmax = w[j];
  thresh = TOL * wmax;
  for (j = 1;j <= ma;j++)
    if (w[j] < thresh) w[j] = 0.0;
  svbksb(u, w, v, ndata, ma, b, a);
  *chisq = 0.0;
  for (i = 1;i <= ndata;i++) {
    (*funcs)(x1[i], x2[i], afunc);
    for (sum = 0.0, j = 1;j <= ma;j++) sum += a[j] * afunc[j];
    *chisq += (tmp = (y[i] - sum) / sig[i], tmp * tmp);
  }
  free_dvector(afunc, 1, ma);
  free_dvector(b, 1, ndata);
}
#undef TOL

void svdvar(double **v, int ma, double w[], double **cvm) {
  int k, j, i;
  double sum, *wti;

  wti = dvector(1, ma);
  for (i = 1;i <= ma;i++) {
    wti[i] = 0.0;
    if (w[i]) wti[i] = 1.0 / (w[i] * w[i]);
  }
  for (i = 1;i <= ma;i++) {
    for (j = 1;j <= i;j++) {
      for (sum = 0.0, k = 1;k <= ma;k++) sum += v[i][k] * v[j][k] * wti[k];
      cvm[j][i] = cvm[i][j] = sum;
    }
  }
  free_dvector(wti, 1, ma);
}

void svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[]) {
  int jj, j, i;
  double s, *tmp;

  tmp = dvector(1, n);
  for (j = 1;j <= n;j++) {
    s = 0.0;
    if (w[j]) {
      for (i = 1;i <= m;i++) s += u[i][j] * b[i];
      s /= w[j];
    }
    tmp[j] = s;
  }
  for (j = 1;j <= n;j++) {
    s = 0.0;
    for (jj = 1;jj <= n;jj++) s += v[j][jj] * tmp[jj];
    x[j] = s;
  }
  free_dvector(tmp, 1, n);
}

void svdcmp(double **a, int m, int n, double w[], double **v) {
  double pythag(double a, double b);
  int flag, i, its, j, jj, k, l = 0, nm = 0;
  double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

  rv1 = dvector(1, n);
  g = scale = anorm = 0.0;
  for (i = 1;i <= n;i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if (i <= m) {
      for (k = i;k <= m;k++) scale += fabs(a[k][i]);
      if (scale) {
        for (k = i;k <= m;k++) {
          a[k][i] /= scale;
          s += a[k][i] * a[k][i];
        }
        f = a[i][i];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i][i] = f - g;
        for (j = l;j <= n;j++) {
          for (s = 0.0, k = i;k <= m;k++) s += a[k][i] * a[k][j];
          f = s / h;
          for (k = i;k <= m;k++) a[k][j] += f * a[k][i];
        }
        for (k = i;k <= m;k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if (i <= m && i != n) {
      for (k = l;k <= n;k++) scale += fabs(a[i][k]);
      if (scale) {
        for (k = l;k <= n;k++) {
          a[i][k] /= scale;
          s += a[i][k] * a[i][k];
        }
        f = a[i][l];
        g = -SIGN(sqrt(s), f);
        h = f * g - s;
        a[i][l] = f - g;
        for (k = l;k <= n;k++) rv1[k] = a[i][k] / h;
        for (j = l;j <= m;j++) {
          for (s = 0.0, k = l;k <= n;k++) s += a[j][k] * a[i][k];
          for (k = l;k <= n;k++) a[j][k] += s * rv1[k];
        }
        for (k = l;k <= n;k++) a[i][k] *= scale;
      }
    }
    anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }
  for (i = n;i >= 1;i--) {
    if (i < n) {
      if (g) {
        for (j = l;j <= n;j++)
          v[j][i] = (a[i][j] / a[i][l]) / g;
        for (j = l;j <= n;j++) {
          for (s = 0.0, k = l;k <= n;k++) s += a[i][k] * v[k][j];
          for (k = l;k <= n;k++) v[k][j] += s * v[k][i];
        }
      }
      for (j = l;j <= n;j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  for (i = IMIN(m, n);i >= 1;i--) {
    l = i + 1;
    g = w[i];
    for (j = l;j <= n;j++) a[i][j] = 0.0;
    if (g) {
      g = 1.0 / g;
      for (j = l;j <= n;j++) {
        for (s = 0.0, k = l;k <= m;k++) s += a[k][i] * a[k][j];
        f = (s / a[i][i]) * g;
        for (k = i;k <= m;k++) a[k][j] += f * a[k][i];
      }
      for (j = i;j <= m;j++) a[j][i] *= g;
    } else for (j = i;j <= m;j++) a[j][i] = 0.0;
    ++a[i][i];
  }
  for (k = n;k >= 1;k--) {
    for (its = 1;its <= 30;its++) {
      flag = 1;
      for (l = k;l >= 1;l--) {
        nm = l - 1;
        if ((double)(fabs(rv1[l]) + anorm) == anorm) {
          flag = 0;
          break;
        }
        if ((double)(fabs(w[nm]) + anorm) == anorm) break;
      }
      if (flag) {
        c = 0.0;
        s = 1.0;
        for (i = l;i <= k;i++) {
          f = s * rv1[i];
          rv1[i] = c * rv1[i];
          if ((double)(fabs(f) + anorm) == anorm) break;
          g = w[i];
          h = pythag(f, g);
          w[i] = h;
          h = 1.0 / h;
          c = g * h;
          s = -f * h;
          for (j = 1;j <= m;j++) {
            y = a[j][nm];
            z = a[j][i];
            a[j][nm] = y * c + z * s;
            a[j][i] = z * c - y * s;
          }
        }
      }
      z = w[k];
      if (l == k) {
        if (z < 0.0) {
          w[k] = -z;
          for (j = 1;j <= n;j++) v[j][k] = -v[j][k];
        }
        break;
      }
      if (its == 30) nrerror("no convergence in 30 svdcmp iterations");
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
      g = pythag(f, 1.0);
      f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
      c = s = 1.0;
      for (j = l;j <= nm;j++) {
        i = j + 1;
        g = rv1[i];
        y = w[i];
        h = s * g;
        g = c * g;
        z = pythag(f, h);
        rv1[j] = z;
        c = f / z;
        s = h / z;
        f = x * c + g * s;
        g = g * c - x * s;
        h = y * s;
        y *= c;
        for (jj = 1;jj <= n;jj++) {
          x = v[jj][j];
          z = v[jj][i];
          v[jj][j] = x * c + z * s;
          v[jj][i] = z * c - x * s;
        }
        z = pythag(f, h);
        w[j] = z;
        if (z) {
          z = 1.0 / z;
          c = f * z;
          s = h * z;
        }
        f = c * g + s * y;
        x = c * y - s * g;
        for (jj = 1;jj <= m;jj++) {
          y = a[jj][j];
          z = a[jj][i];
          a[jj][j] = y * c + z * s;
          a[jj][i] = z * c - y * s;
        }
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free_dvector(rv1, 1, n);
}

double pythag(double a, double b) {
  double absa, absb;
  absa = fabs(a);
  absb = fabs(b);
  if (absa > absb) return absa*sqrt(1.0 + SQR(absb / absa));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + SQR(absa / absb)));
}

void hpsort(int n, float ra[]) {
  unsigned int i, ir, j, l;
  float rra;

  if (n < 2) return;
  l = (n >> 1) + 1;
  ir = n;
  for (;;) {
    if (l > 1) {
      rra = ra[--l];
    } else {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1) {
        ra[1] = rra;
        break;
      }
    }
    i = l;
    j = l + l;
    while (j <= ir) {
      if (j < ir && ra[j] < ra[j+1]) j++;
      if (rra < ra[j]) {
        ra[i] = ra[j];
        i = j;
        j <<= 1;
      } else j = ir + 1;
    }
    ra[i] = rra;
  }
}

#define MAXIT 60
#define UNUSED (-1.11e30)
double zriddr(double (*func)(double xv), double xt1, double xt2, double xacc) {
  int j;
  double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

  /* printf("UNUSED,MAXIT %10.4e, %d\n",UNUSED,MAXIT); */

  fl = (*func)(xt1);
  fh = (*func)(xt2);
  /* printf("fh,fl %10.4e, %10.4e\n",fh,fl); */
  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
    xl = xt1;
    xh = xt2;
    ans = UNUSED;
    for (j = 1;j <= MAXIT;j++) {
      xm = 0.5 * (xl + xh);
      fm = func(xm);
      s = sqrt(fm * fm - fl * fh);
      if (s == 0.0) return ans;
      xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s);
      if (fabs(xnew - ans) <= xacc) return ans;
      ans = xnew;
      fnew = func(ans);
      if (fnew == 0.0) return ans;
      if (SIGN(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      } else if (SIGN(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      } else if (SIGN(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      } else nrerror("never get here.");
      if (fabs(xh - xl) <= xacc) return ans;
    }
    nrerror("ERROR: zriddr exceeded maximum iterations");
  } else {
    if (fl == 0.0) return xt1;
    if (fh == 0.0) return xt2;
    fprintf(stderr, "\nERROR: zriddr, root must be bracketed by limits");
  }
  return 0.0;
}
#undef MAXIT
#undef UNUSED

void nrerror(char error_text[]) {
  fprintf(stderr, "\nERROR: run-time error:\n");
  printf("%s\n", error_text);
  exit(-1);
}

void bp_filter(double bw, double wc, int nps, int nfft, double beta, fcomplex *bpf) {
  /*
    calculate bandpass FIR filter coefficients using a Kaiser window
    frequencies are normalized 0 to be in the range -PI to PI.

    v1.2 clw 26-Sep-2003
    
    bw		normalized filter bandwidth (range 0 to 2PI)
    wc		normalized center frequency (range -PI to +PI)
    nps		number of samples in the filter (must be odd)
    nfft 	number of points in the output fcomplex bpf array
    beta 	Kaiser window coefficient beta
    bpf		bandpass filter coefficients
   
    filter has zero phase delay since it is centered around sample 0
  */

  fcomplex sc, cv;
  double w0;
  int j;

  printf("bp_filter: bw:%8.4f  wc:%8.4f  nfft: %d  nps: %d   Kaiser beta:%8.3f\n", bw, wc, nfft, nps, beta);
  zero_1d((void *)bpf, nfft, sizeof(fcomplex));

  for (j = -nps / 2; j < 0; j++) {
    w0 = bessi0(beta * sqrt(1.0 - SQR(2.0 * j / nps))) / bessi0(beta);
    sc.re = cos(j * wc);
    sc.im = sin(j * wc);			/* carrier to shift filter */
    cv.re = w0 * sin(j * bw / 2.) / (j * PI);	/* L.P. filter coeff. */
    cv.im = 0.0;
    bpf[j + nfft].re = sc.re * cv.re - sc.im * cv.im;
    bpf[j + nfft].im = sc.re * cv.im + sc.im * cv.re;
  }

  bpf[0].re = bw / (2.*PI);
  bpf[0].im = 0.0;

  for (j = 1; j <= nps / 2; j++) {
    w0 = bessi0(beta * sqrt(1.0 - SQR(2.0 * j / nps))) / bessi0(beta);
    sc.re = cos(wc * j);
    sc.im = sin(wc * j);			/* carrier to shift filter */
    cv.re = w0 * sin(j * bw / 2.) / (j * PI);	/* L.P. filter coeff. */
    cv.im = 0.0;
    bpf[j].re = sc.re * cv.re - sc.im * cv.im;
    bpf[j].im = sc.re * cv.im + sc.im * cv.re;
  }
}

#define BETA1	 .5
void cpx_interp(fcomplex **d1, fcomplex **d2, int nx1, int ny1, int n_ovr, int iflg) {
  /*
    18-Nov-2004 clw
    interpolate a complex image patch, assume there is no carrier in x or y directions:

    d1	input fcomplex data patch
    d2	output fcomplex data patch n_ovr*nr1 columns n_ovr*naz1 rows
    nx1	input patch columns
    ny1	input patch rows
    n_ovr	interpolation factor
    iflg	initialization flag iflg = 0 initializes FFT plans and filter coefficients
  */

  static fcomplex **bpf;
  fcomplex *azbpf, *rbpf, *ctmp;
  double re1, re2, im1, im2;
  double bw = 1.0;
  int fir_len;
  int i, j, i1, j1, nx2, ny2;
#ifdef FFTW
  static fftw_plan pfr, pfaz;
  static fftwnd_plan p1, p2;
#else
  int ndim, isign;
  float *dt;
  unsigned int nfft[3];
#endif

  nx2 = nx1 * n_ovr;		/* size of oversampled SLC patch */
  ny2 = ny1 * n_ovr;

  if (iflg == 1) {
    rbpf = (fcomplex *)calloc_1d(nx2, sizeof(fcomplex));		/* range BPF */
    azbpf = (fcomplex *)calloc_1d(ny2, sizeof(fcomplex));		/* azimuth BPF */
    ctmp = (fcomplex *)calloc_1d(MAX(nx2, ny2), sizeof(fcomplex));
    bpf  =  (fcomplex **)calloc_2d(nx2, ny2, sizeof(fcomplex));	/* 2D BPF */
    bw = bw * TWO_PI / n_ovr;	/* 2D bandpass filter bandwidth normalized to 2 PI */
    fir_len = nx1 / 4 + 1;
    printf("interp. filter bandwidth: %10.3f  FIR length: %d\n", bw / TWO_PI, fir_len);

#ifdef FFTW
    pfr = fftw_create_plan(nx2, FFTW_FORWARD, FFTW_ESTIMATE);
    pfaz = fftw_create_plan(ny2, FFTW_FORWARD, FFTW_ESTIMATE);

    p1 = fftw2d_create_plan(ny2, nx2, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    p2 = fftw2d_create_plan(ny2, nx2, FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);

    bp_filter(bw, 0.0, fir_len, nx2, BETA1, ctmp);
    fftw_one(pfr, (fftw_complex *)ctmp, (fftw_complex *)rbpf);

    bp_filter(bw, 0.0, fir_len, ny2, BETA1, ctmp);
    fftw_one(pfaz, (fftw_complex *)ctmp, (fftw_complex *)azbpf);
#else
    bp_filter(bw, 0.0, fir_len, nx2, BETA1, rbpf);
    four1((float *)rbpf - 1, nx2, FFT_F);

    bp_filter(bw, 0.0, fir_len, ny2, BETA1, azbpf);
    four1((float *)azbpf - 1, ny2, FFT_F);
#endif

    for (i = 0; i < ny2; i++) {			/* generate 2-D look filter */
      re2 = azbpf[i].re;
      im2 = azbpf[i].im;
      for (j = 0; j < nx2; j++) {
        re1 = rbpf[j].re;
        im1 = rbpf[j].im;
        bpf[i][j].re = re1 * re2 - im1 * im2;
        bpf[i][j].im = im1 * re2 + re1 * im2;
      }
    }

    for (i = 0; i < ny2; i++) {
      for (j = 0; j < nx2; j++) {
        bpf[i][j].re *= n_ovr * (double)n_ovr / (double)(nx2 * ny2);
        bpf[i][j].im *= n_ovr * (double)n_ovr / (double)(nx2 * ny2);
      }
    }

#ifdef DEBUG
    wrp_cpx(bpf, 0, nx2, ny2, t1f);
    exit(-1);
#endif

#ifdef FFTW
    free(pfr);
    free(pfaz);
#endif
    free(rbpf);
    free(azbpf);
    free(ctmp);
    return;
  }

  zero_2d((void **)d2, nx2, ny2, sizeof(fcomplex));

  for (i = 0; i < ny1; i++) {			/* cut the deck and over sample */
    i1 = i - ny1 / 2;
    if (i1 < 0)i1 += ny1;
    for (j = 0; j < nx1; j++) {
      j1 = j - nx1 / 2;
      if (j1 < 0)j1 += nx1;
      d2[n_ovr*i1][n_ovr*j1] = d1[i1][j1];
    }
  }

#ifdef DEBUG
  wrp_cpx(d2, 0, nx2, ny2, t1f);
#endif

#ifdef FFTW
  fftwnd_one(p1, (fftw_complex *)d2[0], NULL);
#else
  dt = (float *) & d2[0][0];
  dt--;						/* decrement addresses so that indices start at 1 */
  ndim = 2;
  isign = 1;
  nfft[1] = ny2; 		/* NR FFT initialization */
  nfft[2] = nx2;
  nfft[0] = 0;
  fourn(dt, nfft, ndim, isign);			/* 2D forward FFT of region */
#endif

#ifdef DEBUG
  wrp_cpx(d2, 0, nx2, ny2, t1f);
#endif

  for (i = 0; i < ny2; i++) {			/* multiply by bpf */
    for (j = 0; j < nx2; j++) {
      re1 = d2[i][j].re;
      im1 = d2[i][j].im;
      d2[i][j].re = re1 * bpf[i][j].re - im1 * bpf[i][j].im;
      d2[i][j].im = re1 * bpf[i][j].im + im1 * bpf[i][j].re;
    }
  }

#ifdef DEBUG
  wrp_cpx(d2, 0, nx2, ny2, t1f);
#endif

#ifdef FFTW					/* inverse FFT, image is still with image center at 0,0 */
  fftwnd_one(p2, (fftw_complex *)d2[0], NULL);
#else
  dt = (float *) & d2[0][0];
  dt--;						/* decrement addresses so that indices start at 1 */
  ndim = 2;
  isign = -1;
  nfft[1] = ny2; 		/* NR FFT initialization */
  nfft[2] = nx2;
  nfft[0] = 0;
  fourn(dt, nfft, ndim, isign);			/* 2D inverse FFT of region */
#endif

#ifdef DEBUG
  wrp_cpx(d2, 0, nx2, ny2, t1f);
  exit(-1);
#endif
}


#define MA 5
/* #define DEBUG1 1 */

double peak4(fcomplex **cc, int rwin2, int azwin2, int pwin, double *rpl, double *azpl, int iflg) {
  /*
    subroutine to determine peak of correlation function by polynomial fit
    15-Nov-2004 clw
      cc		correlation of complex intensity data
      rwin2	width of fsnr array
      azwin2	height of fsnr array
      pwin	peak window size
      rpl		interpolated range peak location
      azpl	interpolated azimuth peak locaiton
      iflg        initialize array allocation iflg = 1
  */

  static double **covar;	/* covariance matrices of L.S. */
  static double *rv, *azv, *sig;	/* range, azimuth positions, snr */
  static double *ccv;		/* correlation values */
  static double chisq;		/* chi square errors */
  static double *ac;		/* LS solution values */
  static int *ia;		/* array of flags for which coefficients in the ls to fit */

  double dr, da;			/* range and azimuth fractional pixel offsets */
  double pkz;			/* peak value */
  double cc_ave;		/* average correlation function value */
  double cc_max;		/* maximum correlation value */
  double snr;			/* correlation SNR */

  int i, j;
  int i1, j1, i2, j2;		/* peak indices in the uninterpolated correlation function */
  int ict;			/* counter of the number of valid values */

  if (iflg == 1) {

#ifdef DEBUG
    printf("\npeak4: initializing data array allocation for polynomial fit\n");
#endif

    covar = (double **)calloc_2d(MA + 1, MA + 1, sizeof(double));	/* covariance matrix */
    ac = (double *)calloc_1d(MA + 1, sizeof(double));		/* solutions of least square coefficients, range,azimuth */
    rv = (double *)calloc_1d(pwin * pwin + 1, sizeof(double));	/* ranges of points in L.S. */
    azv = (double *)calloc_1d(pwin * pwin + 1, sizeof(double));	/* azimuth positions of points in L.S. */
    ccv = (double *)calloc_1d(pwin * pwin + 1, sizeof(double));	/* correlation values */
    sig = (double *)calloc_1d(pwin * pwin + 1, sizeof(double));	/* sigma of data measurements */
    ia = (int *)calloc_1d(MA + 1, sizeof(int));			/* flag array determining which parameters to fit */
    for (j = 0; j <= MA; j++)ia[j] = 1;

    return (0.0);
  }

  cc_ave = 0.0;
  cc_max = 0.0;
  i1 = 0;
  j1 = 0;

  for (i = 0; i < azwin2; i++) {
    for (j = 0; j < rwin2; j++) {
      cc_ave += fabs(cc[i][j].re);
      if (cc[i][j].re > cc_max) {
        i1 = i;
        j1 = j;
        cc_max = cc[i][j].re;
      }
    }
  }
  if (cc_max == 0.0) {
    *azpl = 0.0;
    *rpl = 0.0;
    return (0.0);
  }
#ifdef DEBUG1
  printf("correlation function  peak i1: %d  j1: %d pwin: %d:  ave: %10.3e\n", i1, j1, pwin, cc_ave / (azwin2*rwin2 - ict));
#endif
  ict = 1;
  for (i = -pwin / 2; i <= pwin / 2; i++) {		/* fill interp. array */
    i2 = i1 + i;
    if (i2 < 0)i2 += azwin2;			/* correlation function is circular periodic */
    if (i2 >= azwin2)i2 -= azwin2;
    for (j = -pwin / 2; j <= pwin / 2; j++) {
      j2 = j1 + j;
      if (j2 < 0)j2 += rwin2;			/* correlation function is circular periodic */
      if (j2 >= rwin2)j2 -= rwin2;
      if (fabs(cc[i2][j2].re) > 0.0) {
        azv[ict] = i;
        rv[ict] = j;
        ccv[ict] = fabs(cc[i2][j2].re);
        sig[ict] = 1.0;
        cc_ave -= cc[i2][j2].re;		/* subtract peak values from average */
#ifdef DEBUG1
        printf("%3d", nint(100.*ccv[ict] / cc_max));
        /*printf("i: %2d  j: %2d  i1: %3d  j1: %3d  i2: %4d  j2: %4d ccv: %10.3e\n",i,j,i1,j1,i2,j2,ccv[ict]); */
#endif
        ict++;
      }
    }
#ifdef DEBUG1
    printf("\n");
#endif
  }

  ict--;

  if (ict >= MA) {
    lfit(rv, azv, ccv,  sig, ict, ac, ia, MA, covar, &chisq,  fp_peak);  	/* L.S. fit */
    da = -ac[3] / (2.*ac[5]);
    dr = -ac[2] / (2.*ac[4]);
    pkz = ac[1] + dr * ac[2] + da * ac[3] + dr * dr * ac[4] + da * da * ac[5];

    if (j1 < rwin2 / 2) *rpl = j1 + dr;		/* wrap offsets */
    else *rpl = j1 + dr - rwin2;

    if (i1 < azwin2 / 2) *azpl = i1 + da;
    else *azpl = i1 + da - azwin2;

    cc_ave = cc_ave / (azwin2 * rwin2 - ict);
    snr = pkz / cc_ave;
  } else {
    da = 0.0;
    dr = 0.0;
    snr = 0.0;
  }
#ifdef DEBUG1
  printf("coefficients: %12.4e  %12.4e %12.4e  %12.4e %12.4e\n", ac[1], ac[2], ac[3], ac[4], ac[5]);
  printf("peak index range:   %6d  frac: %8.6f  total: %10.4f \n", j1, dr, *rpl);
  printf("peak index azimuth: %6d  frac: %8.6f  total: %10.4f\n", i1, da, *azpl);
  printf("peak value: %10.3e  ave: %10.3e  SNR: %10.3f\n\n", pkz, cc_ave, snr);
#endif
  return (snr);
}

void phg_est(fcomplex **cpx, int roff, int azoff, int nr, int naz, double *rpg, double *azpg) {
  /*
     calculate range and azimuth phase gradients/sample for an interferogram chip
     clw 18-Jan-2002
     v1.1	add test for null values clw 11-Nov-2005
   
     cpx	complex 2D array
     roff	range offset to start of window to evaluate gradients
     azoff	azimuth offset to start of window to evaluate gradients
     nr		range widow size
     naz	azimuth window size
     rpg	range phase gradient   (radians/sample)
     azpg	azimuth phase gradient (radians/sample)
   
  */

  int k, n, k1, k2, n1, n2;

  double psr_re, psr_im;
  double psaz_re, psaz_im;

  k1 = 1 + azoff;	/* begining and ending indices for detection of phase gradients */
  k2 = k1 + naz - 1;

  n1 = 1 + roff;
  n2 = n1 + nr - 1;

  psr_re = 0.0;
  psr_im = 0.0; 	/* init sums */
  psaz_re = 0.0;
  psaz_im = 0.0;

  for (k = k1; k < k2; k++) {
    for (n = n1; n < n2; n++) {
      if (cpx[k][n].re == 0.0 || cpx[k-1][n].re == 0.0)continue;	/* prevent error in gradient calc at the edges */
      if (cpx[k][n-1].re == 0.0)continue;
      psr_re += (cpx[k][n].re * cpx[k][n - 1].re + cpx[k][n].im * cpx[k][n - 1].im);
      psr_im += (cpx[k][n].im * cpx[k][n - 1].re - cpx[k][n].re * cpx[k][n - 1].im);

      psaz_re += (cpx[k][n].re * cpx[k - 1][n].re + cpx[k][n].im * cpx[k - 1][n].im);
      psaz_im += (cpx[k][n].im * cpx[k - 1][n].re - cpx[k][n].re * cpx[k - 1][n].im);
    }

  }

  *rpg = 0.0;
  *azpg = 0.0;

  if (psr_re != 0.0)
    *rpg = atan2(psr_im, psr_re);	/* evaluate phase gradients */
  if (psaz_re != 0.0)
    *azpg = atan2(psaz_im, psaz_re);
}

void phg_apply(fcomplex **cpx, int nr, int naz, double rpg, double azpg) {
  /*
    apply a phase gradient to a rectangular patch v1.0 24-Dec-2004 clw

    cpx	fcomplex 2D array
    nr	number of pixels in range
    naz   number of lines in azimuth
    rpg	range phase gradient/pixel
    azpg  azimuth phase gradient/pixel

    v1.1  use identities for cos(a+b) and sin(a+b) to improve efficiency  30-Mar-2005 clw

  */
  double re1, im1;
  double re2, im2;
  double cr, sr;
  double c1, s1;

#ifdef DIRECT_METHOD
  double phi, c2, s2;
#endif

  int k, n;

  cr = cos(rpg);			/* cos and sin of the phase step/range pixel */
  sr = sin(rpg);

  for (k = 1; k < naz; k++) {
    c1 = cos(k * azpg); 		/* initial cos and sin at the start of the line */
    s1 = sin(k * azpg);

    for (n = 1; n < nr; n++) {
      re1 = cpx[k][n].re;
      im1 = cpx[k][n].im;

      re2 = c1 * cr - s1 * sr;		/* use cos(a+b) and sin(a+b) identities to calculate sin and cos from previous values */
      im2 = s1 * cr + c1 * sr;
      c1 = re2;				/* update values for next iteration */
      s1 = im2;

#ifdef DIRECT_METHOD
      phi = n * rpg + k * azpg;
      c2 = cos(phi);
      s2 = sin(phi);
      if (k % 100 == 0 && n % 100 == 0) {
        printf("phg_apply k: %6d  n: %6d cos: %12.9f  sin: %12.9f  err_cos: %12.5e  err_sin: %12.5e\n", k, n, re2, im2, re2 - c2, im2 - s2);
      }
#endif

      cpx[k][n].re = re1 * re2 - im1 * im2;
      cpx[k][n].im = im1 * re2 + re1 * im2;
    }
  }
}

int rd_cpx(fcomplex **a, int width, int nls, int roff, int azoff, int nr, int naz, int dat_type, FILE *datf) {
  /*
      read fcomplex data covering a patch from a data file, either short or floating point complex
      patches are read from a region starting at (line:azoff1, sample: roff1) in the large file.
      nominal patch size is naz lines x nr samples
      return value is lines actually read from the file and is set to 0 if no lines read
   
      all entries in the fcomplex data array are set to (0.0,0.0) if no data read for that location
   
      20-Feb-2002 clw
   
      a		output 2D fcomplex point phase derived from complex data
      width	width of input data file, number of samples/line
      nls		number of lines in the data file
      roff         range offset to starting range sample
      azoff	offset to starting azimuth line of data subset that will be processed
      nr		number of samples/line of the patch
      naz          number of azimuth lines in the patch
      dat_type	type=0 floating point complex, type=1 short complex type=2, float data, set imaginary part to 0
      datf		floating point format FILE data structure
   */

  fcomplex *xtmp;
  scomplex *stmp;
  float *ftmp;

  int i, j, jj, naz1, nr1;
  double scale = .0001;

  xtmp = (fcomplex *)malloc(width * sizeof(fcomplex));
  stmp = (scomplex *)malloc(width * sizeof(scomplex));
  ftmp = (float *)malloc(width * sizeof(float));

  if (xtmp == NULL || stmp == NULL || ftmp == NULL) {
    fprintf(stderr, "\nERROR subroutine rd_cpx: memory allocation error for tmp arrays\n\n");
    exit( -1);
  }

  if (roff >= width) {
    fprintf(stderr, "\nERROR subroutine rd_cpx: starting range sample offset exceeds width: %d  %d\n\n", roff, width);
    exit( -1);
  }

  if (azoff >= nls) {
    fprintf(stderr, "\nERROR subroutine rd_cpx: starting azimuth line exceeds number of lines in the file: %d  %d\n\n", azoff, nls);
    exit( -1);
  }

  if (roff < 0 || azoff < 0) {
    fprintf(stderr, "\nERROR subroutine rd_cpx: invalid starting offset roff or azoff < 0: roff: %d   azoff: %d\n", roff, azoff);
    exit( -1);
  }

  if ((naz + azoff) > nls)
    naz1 = nls - azoff;	/* check if blank lines need to be written at the end of the array */
  else
    naz1 = naz;

  if ((nr + roff) > width)
    nr1 = width - roff;	/* check if blank samples need to be written at the end of the lines */
  else
    nr1 = nr;

  switch (dat_type) {
  case FCMPLX_DATA:
    fseek(datf, (off_t)sizeof(fcomplex)*azoff*width, SEEK_SET);
    break;
  case SCMPLX_DATA:
    fseek(datf, (off_t)sizeof(scomplex)*azoff*width, SEEK_SET);
    break;
  case FLOAT_DATA:
    fseek(datf, (off_t)sizeof(float)*azoff*width, SEEK_SET);
    break;
  default:
    fprintf(stderr, "\nERROR rd_cpx: invalid data file type: %d\n\n", dat_type);
    exit( -1);
  }

  for (i = 0; i < naz1; i++) {
    switch (dat_type) {
    case FCMPLX_DATA:
      fread((char *)xtmp, sizeof(float), 2*width, datf);
      if (feof(datf)) {
        fprintf(stderr, "\nERROR subroutine rd_cpx: unexpected end of data file at line: %d\n", azoff + i);
        exit( -1);
      }
      for (j = 0; j < nr1; j++) {
        jj = j + roff;				/* offset in range */
        a[i][j] = xtmp[jj];
      }
      break;

    case SCMPLX_DATA:
      fread((char *)stmp, sizeof(short), 2*width, datf);
      if (feof(datf)) {
        fprintf(stderr, "\nERROR subroutine rd_cpx: unexpected end of data file at line: %d\n", azoff + i);
        exit( -1);
      }
      for (j = 0; j < nr1; j++) {
        jj = j + roff;				/* offset in range */
        a[i][j].re = stmp[jj].re * scale;
        a[i][j].im = stmp[jj].im * scale;
      }
      break;

    case FLOAT_DATA:
      fread((char *)ftmp, sizeof(float), width, datf);
      if (feof(datf)) {
        fprintf(stderr, "\nERROR subroutine rd_cpx: unexpected end of data file at line: %d\n", azoff + i);
        exit( -1);
      }
      for (j = 0; j < nr1; j++) {
        jj = j + roff;				/* offset in range */
        a[i][j].re = ftmp[jj];
        a[i][j].im = 0.0;
      }
      break;

    default:
      fprintf(stderr, "\nERROR rd_cpx: invalid data file type: %d\n\n", dat_type);
      exit( -1);
    }

    if (nr1 < nr) {					/* fill in missing range samples */
      for (j = nr1; j < nr; j++) {
        a[i][j].re = 0.0;
        a[i][j].im = 0.0;
      }
    }
  }

  if (naz1 < naz) {					/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++) {
        a[i][j].re = 0.0;
        a[i][j].im = 0.0;
      }
    }
  }

  free(xtmp);
  free(stmp);
  free(ftmp);
  return (naz1);
}

void wrp_cpx(fcomplex **tcpx, int ip, int nr, int naz, FILE *cpxf) {
  /*
    subroutine to write a patch of fcomplex data
    v1.1 10-Jan-2002 clw

     tcpx		2-D array of fcomplex point data
     ip		patch number
     nr		number of range samples/line
     naz		number of azimuth lines
     cpxf		file data structure for output file
  */

  off_t off1;
  size_t nrb;
  int i;

#ifdef fread
  int j;
  unsigned char *b1, btmp;
#endif

  off1 = (off_t)(2 * nr * naz);
  fseek(cpxf, ip*sizeof(float)*off1, SEEK_SET);/* move forward to this patch */

  for (i = 0; i < naz; i++) {
    nrb = fwrite((char *)tcpx[i], sizeof(float), 2 * nr, cpxf);
    if (nrb != (size_t)2*nr) {
      fprintf(stderr, "\nERROR wrp_cpx: error writing patch of fcomplex data ip: %d\n\n", ip);
      exit(-1);
    }
  }

#ifdef fread					/* unswap the bytes after write if necessary */
  for (i = 0; i < naz; i++) {
    for (j = 0; j < nr; j++) {
      b1 = (unsigned char *) & (tcpx[i][j].re);

      btmp = b1[0];
      b1[0] = b1[3];
      b1[3] = btmp;
      btmp = b1[1];
      b1[1] = b1[2];
      b1[2] = btmp;

      b1 = (unsigned char *) & (tcpx[i][j].im);

      btmp = b1[0];
      b1[0] = b1[3];
      b1[3] = btmp;
      btmp = b1[1];
      b1[1] = b1[2];
      b1[2] = btmp;
    }
  }
#endif

  fflush(cpxf);
}

void wrp_flt(float **tflt, int ip, int nr, int naz, FILE *fltf) {
  /*
    subroutine to write a patch of float data
    v1.1 5-Jul-2006 clw

     tflt		2-D array of float point data
     ip		patch number
     nr		number of range samples/line
     naz		number of azimuth lines
     fltf		file data structure for output file
  */

  off_t off1;
  size_t nrb;
  int i;

#ifdef fread
  int j;
  unsigned char *b1, btmp;
#endif

  off1 = (off_t)(nr * naz);
  fseek(fltf, ip*sizeof(float)*off1, SEEK_SET);/* move forward to this patch */

  for (i = 0; i < naz; i++) {
    nrb = fwrite((char *)tflt[i], sizeof(float), nr, fltf);
    if (nrb != (size_t)nr) {
      fprintf(stderr, "\nERROR wrp_flt: error writing patch of float data ip: %d\n\n", ip);
      exit(-1);
    }
  }

#ifdef fread					/* unswap the bytes after write if necessary */
  for (i = 0; i < naz; i++) {
    for (j = 0; j < nr; j++) {
      b1 = (unsigned char *) & (tflt[i][j]);
      btmp = b1[0];
      b1[0] = b1[3];
      b1[3] = btmp;
      btmp = b1[1];
      b1[1] = b1[2];
      b1[2] = btmp;
    }
  }
#endif

  fflush(fltf);
}

void fp_peak(double x1, double x2, double p[]) { /* last parameter is the number of terms */
  p[1] = 1.0;
  p[2] = x1;
  p[3] = x2;
  p[4] = x1 * x1;
  p[5] = x2 * x2;
}

#ifdef WIN32
#include<sys/time.h>
static clock_t c1;

void start_timing() {
  c1 = clock();
}

void stop_timing() {
  static clock_t c2;
  c2 = clock();
  printf("\nreal elapsed time (s): %10.3f\n\n", (double)(c2 - c1) / CLOCKS_PER_SEC);
}
#else
#include <sys/types.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
#include <limits.h>

static struct tms buffer;			/* global variables */
static int user_time, system_time, start_time;

void start_timing() {
  start_time = (int)times(&buffer);
  user_time = (int)buffer.tms_utime;
  system_time = (int)buffer.tms_stime;
}

void stop_timing() {
  static int clk_tck, end_time;

  clk_tck = (int)sysconf(_SC_CLK_TCK);
  end_time = (int)times(&buffer);
  user_time = (int)buffer.tms_utime - user_time;
  system_time = (int)buffer.tms_stime - system_time;

  printf("\nuser time (s):    %10.3f\n", (double)user_time / clk_tck);
  printf("system time (s):  %10.3f\n", (double)system_time / clk_tck);
  printf("elapsed time (s): %10.3f\n\n", (double)(end_time - start_time) / clk_tck);
}
#endif
