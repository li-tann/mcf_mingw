/*
  Copyright Gamma Remote Sensing AG 2008
  v1.1  23-Nov-2006 clw remove definitions of TIFF, JPEG, and PNG
  v1.2  16-May-2008 clw	added support for LARGEFILES on Win32 (MSDOS)
*/
#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <string.h>
#include <math.h> 
#include <errno.h>
#include <locale.h>

#include <sys/types.h>
#include <sys/stat.h>

extern void start_timing();		/* timing routines */
extern void stop_timing();

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

#ifndef MM_PAGE_SZ
#define MM_PAGE_SZ 4096		/* memory map page size (bytes) */
#endif

#ifndef PI
#define PI	3.14159265358979323846
#endif

#define TWO_PI  6.28318530717958647693
#define SQRT2   1.41421356237 	/* square root of 2 */
#define RTD	57.2957795131	/* radians to degrees */
#define DTR	.0174532925199	/* degrees to radians */
#define C	2.99792458e8
#define FTPM    3.2808399;	/* feet per meter */

/*** Image Format Definitions ***/
#define UNKNOWN 	0	
#define SUN		1	/* SUN Raster file*/
#define BMP		2	/* Windows BMP */

/***  Color Table Definitions ***/
#define USER_DEFINED	0
#define GREY		1		
#define MPH		2
#define PHASE		3
#define PHUNW		4
#define PH_4BIT		5
#define GREY_4BIT	6
#define CH_NEUT_LSNR    7		/* charges+neutrons+LSNR */
#define NEUT_LSNR 	8		/* neutrons+low SNR */
#define CH_NEUT_LSNR_CUTS	9  	/* charges+neutrons+LSNR+cuts*/
#define CUTS_LSNR	10		/* cuts+low SNR */
#define CUTS_LAWN	11		/* cuts+lawn */
#define LAWN_LSNR	12		/* lawn+low SNR */
#define LAWN_LSNR_CUTS  13		/* lawn+low SNR + cuts */

/*** definitions added 5-Mar-2001 ***/
#define MPH_OLD		20
#define MPH_HLS		21
#define MPH_HSV		22
#define MPH_SIN		23
#define MPH_TOPO1	24
#define PHASE_OLD	30	
#define PHASE_HLS	31
#define PHASE_HSV	32
#define PHASE_SIN	33
#define PHASE_TOPO1	34

/* Phase unwrapping flag bits */
#ifndef PHASE_UNWRAP_DEF
#define PLUS 1		/* plus residue */
#define MINU 2		/* negative residue */
#define CHG  3		/* Logical OR of PLUS and MINU (PLUS|MINU), used in logical expressions to test if a point is a charge */
#define GUID 4		/* guide point (neutron) */
#define LSNR 8		/* Low SNR */
#define VIST 16		/* visited charge, already part of any tree */ 
#define BRPT 32		/* visited charge already part of the current tree */
#define CUT  64		/* branch cut pixel */
#define LAWN 128	/* marks regions that have been unwrapped */
#define TREE 128	/* temporary flag marks point as part of current network of trees, later this bit reused for LAWN */
#define LSNR_VIST 4	/* temporary flag used to mark LOW SNR regions */
#define NOT_BRPT  223	/* unmark charge as part of current tree (!BRPT = 255 - BRPT)*/
#define PHASE_UNWRAP_DEF  1
#endif

#define WHITE     192
#define BLACK     193
#define RED       194
#define BLUE      195
#define SALMON    196 
#define CYAN      197
#define GREEN     198
#define YELLOW    199
#define MAGENTA   200

/*** Memory Map Page Size */
#ifndef MM_PAGE_SZ
#define MM_PAGE_SZ 4096		/* memory map page size (bytes) */
#endif

/*** Earth Ellipsoid Parameters */
#ifndef ELLIPSE_DEF
#define WGS84_MAJOR 6378137.0
#define WGS84_MINOR 6356752.3141
#define GEM6_MAJOR 6378144.0
#define GEM6_MINOR 6356759.0
#define INTERNATIONAL_MAJOR 6378140.0 
#define INTERNATIONAL_MINOR 6356755.0
#define GRS80_MAJOR WGS84_MAJOR
#define GRS80_MINOR WGS84_MINOR
#define ELLIPSE_DEF    1
#endif

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
static double nintarg=0.0; 
#define nint(a) ( ((nintarg=(a)) >= 0.0 )?(int)(nintarg+0.5):(int)(nintarg-0.5) )
#define aint(a) ((double)(int)(a))
#define log2(a) (log(a)/.693147181)
#define exp2(a) (pow(2.0,a))
#define SUN_MACRO	1
#endif

#ifndef STRUCTURES
typedef struct{double re,im;} dcomplex;	/* double precision complex data type */ 
typedef struct {float re,im;} fcomplex; /* single precision complex data type */
typedef struct {short re,im;} scomplex;	/* short complex data type */
typedef struct {unsigned char re,im;} bcomplex;	/* byte data complex structure */
typedef struct {char re,im;} bcpx;		/* signed byte data complex structure */

typedef struct{double x,y,z;} VEC;	/* vector structure cartesian (X,Y,Z)*/
typedef struct{double t,c,n;} VEC_TCN;	/* vector structure (Track, Cross-Track, Normal) */
typedef struct{double s,c,h;} VEC_SCH;	/* vector structure (Track, Cross-Track, Height) */
typedef struct{double n,e,d;} VEC_NED;	/* vector structure (North, East, Down) */
typedef struct{double v,c,p;} VEC_VCP;	/* vector structure (Velocity, Cross-Track, Perpendicular) */
typedef struct{float az,cr;} Slope; 	/* azimuth,range slope */     

typedef struct{			/* state vector structure */
   VEC pos;			/* position vector */
   VEC vel;			/* velocity vector */
} STATE;			/* state vector */

typedef struct{
  int mjd;		/* modified Julian date, which is (julian date - 2433282.5) example: 1.1.1990 = 14610 */
  int utc;		/* time in msec since midnight UTC */
  int orb_num;		/* orbit number */
  VEC pos;		/* sensor position (xyz) Conventional Terrestrial System (km) */
  VEC vel;		/* sensor velocity (xyz) Conventional Terrestrial System (km/sec) */
  int TAI_UTC_delta;	/* value in seconds of the difference between the Internation Atomic Time and UTC */
} ORRM;

typedef double **MAT;
typedef struct{double lat,lon,alt;} POSITION;	/* position in lat/lon and altitude (WGS-84) */

typedef struct{
  char reckey[7];	/* record key STATE for PRC */ 
  double start;		/* start date of arc in days sincd 1.1.2000 12h */
  double end;		/* end date of arc in days since 1.1. 2000 12h */
  double tdtutc;	/* time difference TDT-UTC in seconds */
 } PRC_HDR; 

typedef struct{
  char reckey[7];		/* record key STTERR for Terrestrial frame */ 
  double ttagd;			/* Julian days since 1.1.2000 12h in Terrestrial Dynamic Time */
  double ttagms;		/*  microseconds since 0:00 TDT */
  double xsat,ysat,zsat; 	/* coordinates in mm Convential Terrestrial system */
  double xdsat,ydsat,zdsat;	/* velocity in micrometers/sec Convential Terrestrial system */
  double roll,pitch,yaw;	/* roll, pitch, yaw angles in deg. */
  double radcor;		/* radial correction for the orbits as derived from altimeter */
 } PRC;
#define STRUCTURES 	1
#endif

#define LR_NORMAL 1		/* normal right to left display */
#define LR_REV	 -1		/* display reverse of each line, left to right */

#ifndef nint
static double nintarg; 
#define nint(a) ( ((nintarg=(a)) >= 0.0 )?(int)(nintarg+0.5):(int)(nintarg-0.5) )
#endif

#define SQR(a)    ( (a) * (a) )
#define CS(a,b)   ( ((b) >= 0.0) ? (a) : (-a))
#define ABS(a)    ( ((a) > 0) ? (a) : (-a) )
#define SGN(a)	  (((a)<0) ? -1 : 1)

#ifndef MAX
#define MAX(a,b)  ( ( (a) > (b) ) ? (a) : (b) )
#endif

#ifndef MIN
#define MIN(a,b)  ( ( (a) < (b) ) ? (a) : (b) )
#endif

#define Max(a,b)  ( ( (a) > (b) ) ? (a) : (b) )
#define Min(a,b)  ( ( (a) < (b) ) ? (a) : (b) )

#define ASSERT(x) if (!(x)) fprintf(stderr," Assert failed: x\n");
#define ALLOC(ptr, type, n)  ASSERT(ptr = (type *)malloc((n)*sizeof(type)))

typedef struct{short x,y;}Xpoint;
typedef struct{
    unsigned char	red;
    unsigned char	green;
    unsigned char	blue;
}COLOR; 

typedef struct Point2Struct {	/* 2d point */
	double x, y;
	} Point2;
typedef Point2 Vector2;

typedef struct {		/* window: a discrete 2-D rectangle */
    int x0, y0;			/* xmin and ymin */
    int x1, y1;			/* xmax and ymax (inclusive) */
} Win;

typedef struct {		/* a polygon edge */
    double x;	/* x coordinate of edge's intersection with current scanline */
    double dx;	/* change in x with respect to y */
    int i;	/* edge number: edge i goes from pt[i] to pt[i+1] */
} Edge;

#ifdef GTK
#ifdef MSDOS
#ifdef G_OS_WIN32
#undef G_OS_WIN32
#endif
#include <glib.h>
#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h>
#define ZOOM_MAG_FACTOR 3
#define MAG_WIN_SIZE 120
#else
#include <gtk/gtk.h>
#include <gdk/gdkx.h>
#include <gdk/gdkkeysyms.h>
#define ZOOM_MAG_FACTOR 3
#define MAG_WIN_SIZE 120
#endif

#define WIN_WIDTH_MIN 	256
#define WIN_HEIGHT_MIN	256
#define WIN_BORDER 4

#define D_NONE		-1	/* no data */
#define D_BYTE		0	/* single unsigned byte */
#define D_SHORT		1	/* 2 byte short integer */
#define D_INT           2	/* 4 byte long integer */
#define D_FLOAT		3	/* single precision 4 bytes */
#define D_DOUBLE	4	/* double precision 8 bytes */
#define D_BCOMPLEX	5	/* byte complex, byte real, byte imaginary */
#define D_SCOMPLEX	6	/* short complex, short real, short imaginary */
#define D_FCOMPLEX	7	/* single precision complex, float real, float imaginary */
#define D_DCOMPLEX	8	/* double precision complex, double real, double imaginary */
#define D_RGB		9	/* RGB triplets 3 bytes/pixel */
#define D_SSLC		10	/* complex SLC (short integer) displayed in as detected greyscale */
#define D_FSLC		11	/* complex SLC (float) displayed in as detected greyscale */	
#define D_DSLC		12	/* complex SLC (double) displayed in as detected greyscale */
#define D_DB		13	/* float data displayed as dB values */
#define D_CC		14	/* float data correlation coefficient */
#define D_PHASE		15      /* phase */
#define D_HEIGHT	16      /* height */
#define D_FLAGS		17	/* phase unwrapping flags */
#define D_FFT		18	/* 2-D FFT of floating point complex data */
#define D_DEM_PAR       19	/* display DEM heights and location with DEM parameters */
#define D_BCPX		20	/* display signed byte complex data */
#define D_FFT_SC	21	/* display 2-D FFT of short complex data */

#define MAX_GDK_COLOR_TABLES   14

typedef struct{
  char *label;		/* image label (usually file name) */
  FILE *fstruct;	/* pointer to file structure for the data */
  int data_type;	/* data type D_BYTE, D_SHORT,... */
  int width;		/* samples/line in the input data file and pixmap */
  int height;		/* number of lines in the pixmap */
  int start_line;	/* starting line in the file (starts with 0) */
  int x_off;		/* x offset (across) of the byte data */
  guchar *data;		/* byte image data passed into structure */
}ImageInfo;

typedef struct{
  int zmf;		/* zoom magnification factor */
  int zws_x;		/* zoom window size in x */
  int zws_y;		/* zoom window size in y */
  int zbpp;		/* bytes per pixel (1, 3(RGB)) */			
  guchar *zdata;	/* magnified data */
}ZoomInfo;
#endif
