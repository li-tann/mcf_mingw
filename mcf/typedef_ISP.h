/*
    Copyright 2008 Gamma Remote Sensing clw/uw
    v2.4 13-Aug-2003 clw  added Doppler rate and acceleration polynomials
    v2.5 20-Sep-2003 clw  modified SLC_PAR sensor and date string lengths
    v2.6 23-Jun-2004 clw  added KOE data structure
    v2.7 27-Feb-2006 clw  added ASAR WSS image parameters and WSS Doppler ADS parameters
    v2.8 11-Feb-2008 clw  move to 6 term offset model
    v2.9 16-May-2008 clw  added support for LARGEFILES on Win32 (MSDOS)
    v3.0 25-Aug-2008 clw  added GPRI data structure 
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

#ifndef G_OS_WIN32
#include <sys/mman.h>
#else
#endif

#ifdef FFTW	/* FFTW v2.1.5 subroutines */
#ifdef SFFTW
#include <sfftw.h>
#include <srfftw.h>
#endif
#ifndef SFFTW
#include <fftw.h>
#include <rfftw.h>
#endif
#endif

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
#define DATA_TYPES    1
#endif

#ifndef MM_PAGE_SZ
#define MM_PAGE_SZ 4096	/* memory map page size (bytes) */
#endif

#define MAX_STATE 64	/* maximum number of state vectors in the SLC parameter file */
#define MAR       6    	/* range offset fit linear coefficients in   a0 + a1*r + a2*az +a3*r*az + a4*r*r + a5*az*az */
#define MAAZ      6    	/* azimuth offset fit linear coefficients in a0 + a1*r + a2*az+ a3*r*az + a4*r*r + a5*az*az */

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

#ifndef ELLIPSE_DEF
#define WGS84_MAJOR 6378137.0
#define WGS84_MINOR 6356752.3141
#define WGS84_E2 6.6943800367e-3

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
static double nintarg = 0.0;
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
} fcomplex;	/* single precision complex data type */
typedef struct {
  short re, im;
} scomplex;	/* short complex data type */
typedef struct {
  unsigned char re, im;
} bcomplex;	/* byte data complex structure */

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

typedef struct {	/* state vector structure */
  VEC pos;		/* position vector */
  VEC vel;		/* velocity vector */
} STATE;		/* state vector */

typedef struct {
  int mjd;		/* modified Julian date, which is (julian date - 2433282.5) example: 1.1.1990 = 14610 */
  int utc;		/* time in msec since midnight UTC */
  int orb_num;		/* orbit number */
  VEC pos;		/* sensor position (xyz) Conventional Terrestrial System (km) */
  VEC vel;		/* sensor velocity (xyz) Conventional Terrestrial System (km/sec) */
  int TAI_UTC_delta;	/* value in seconds of the difference between the Internation Atomic Time and UTC */
} ORRM;

typedef double **MAT;
typedef struct {
  double lat, lon, alt;
} POSITION;	/* position in lat/lon and altitude (WGS-84) */

typedef struct {
  char reckey[7];	/* record key STATE for PRC */
  double start;		/* start date of arc in days sincd 1.1.2000 12h */
  double end;		/* end date of arc in days since 1.1. 2000 12h */
  double tdtutc;	/* time difference TDT-UTC in seconds */
} PRC_HDR;

typedef struct {
  char reckey[7];		/* record key STTERR for Terrestrial frame */
  double ttagd;			/* Julian days since 1.1.2000 12h in Terrestrial Dynamic Time */
  double ttagms;		/*  microseconds since 0:00 TDT */
  double xsat, ysat, zsat; 	/* coordinates in mm Convential Terrestrial system */
  double xdsat, ydsat, zdsat;	/* velocity in micrometers/sec Convential Terrestrial system */
  double roll, pitch, yaw;	/* roll, pitch, yaw angles in deg. */
  double radcor;		/* radial correction for the orbits as derived from altimeter */
} PRC;

typedef struct {	/* Keplerian Orbital Elemens */
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

#define STRUCTURES 	1
#endif

typedef struct {		/* single look complex image parameters */
  char title[128];		/* ascii string with title of the scene */
  char sensor[64];		/* sensor name (RADARSAT, SRL-1, SRL-2, ERS-1, ERS-2, JERS-1,...) */
  char date[32];		/* date in form: YYYY MM DD hh mm ss.ttt UTC */
  double t0, t1, t2;		/* time of image start, center, end  UTC seconds since start of day */
  double tazi;			/* time per azimuth line (s) */
  int hdrsz;			/* header size in bytes for each image line */
  int nr;			/* number of range pixels/line in the image */
  int naz;			/* number of range lines in the scene */
  int nlr;			/* number of range looks, for SLC = 1 */
  int nlaz;			/* number of azimuth looks, for SLC = 1 */
  int nstate;			/* number of state vectors */
  char image_format[16];	/* image format FCOMPLEX, SCOMPLEX, FLOAT, SHORT, BYTE */
  char image_geometry[16];	/* image geometry type SLANT_RANGE, GROUND_RANGE, GEOCODED */
  double rpssf;			/* range pixel spacing scale factor, without resampling rpssf=1.0 */
  double azpssf;		/* azimuth pixel spacing scale factor, without resampling azpssf=1.0 */
  double c_lat, c_lon;		/* latitude, longitude of scene center in decimal degrees */
  double c_head;		/* sub-satellite track heading at scene center (decimal-degrees) */
  double rps;			/* slant range pixel spacing  (meters) */
  double azps;			/* azimuth along track pixel spacing (meters) */
  double r0, r1, r2;		/* near, center, far slant or ground range of image (meters) */
  double sr0[6];		/* first slant range (ground range) polynomial coefficients
				   sr0[0] contains the reference orbit time for the polynomial,
				   slant range = sr0[1] + sr0[2]*(GR-r0) + sr0[3]*(GR-r0)**2 +
						 sr0[4]*(GR-r0)**3 + sr0[5]*(GR-r0)**4
				  (r0 is the ground range of the first pixel) */
  double sr1[6];		/* center slant range (ground range) polynomial coefficients
				   sr1[0] contains the reference orbit time for the polynomial,
				   slant range = sr1[1] + sr1[2]*(GR-r0) + sr1[3]*(GR-r0)**2 +
						sr1[4]*(GR-r0)**3 + sr1[5]*(GR-r0)**4
				   (r0 is the ground range of the first pixel) */
  double sr2[6];		/* last slant range (ground range) polynomial coefficients
				   sr2[0] contains the reference orbit time for the polynomial,
				   slant range = sr2[1] + sr2[2]*(GR-r0) + sr2[3]*(GR-r0)**2 +
						sr2[4]*(GR-r0)**3 + sr2[5]*(GR-r0)**4
				   (r0 is the ground range of the first pixel) */
  double c_inc;			/* incidence angle at the center of the scene (deg.) */
  char az_deskew[8];		/* azimuth deskew (ON, OFF) */
  double az_ang;	        /* nominal azimuth antenna angle (decimal degrees CW about N,
                                   right looking SAR 90.0, left looking: -90.0) */
  double fcen;			/* radar carrier center frequency (Hz) */
  double fadc;			/* sample rate  of radar analog to digital converter (Hz) */
  double chbw;			/* radar range chirp bandwidth (Hz) */
  double prf;			/* radar pulse repetition frequency (Hz) */
  double azpbw;			/* 3 dB azimuth processing bandwidth (Hz) */
  double fdp[4];		/* doppler centroid polynomial coefficients
				   fdp[0] + fdp[1]*(r-r1) + fdp[2]*(r-r1)**2 + fdp[3]*(r-r1)**3
				   (r is the slant range and r1 is the slant range at the center of the image */
  double fdp_dot[4];		/* derivative w.r.t. along-track time of each of the terms in fdp[] */
  double fdp_ddot[4];		/* second derivative w.r.t. along-track time of each of the terms in fdp[] */
  double rx_gain;		/* receiver gain (dB) */
  double cal_gain;		/* calibration gain (dB) */

  double rdist, re_cen;		/* distance of SAR sensor from earth center at center scene,
				   center scene geocentric radius (m)*/
  double el_major, el_minor;	/* earth ellipsoid semi-major, semi-minor axises (m) */
  double t_state;		/* UTC time (sec) since start of day for first state vector */
  double tis;			/* time interval between state vectors (s) */
  STATE state[MAX_STATE];	/* maximum of MAX_STATE state vectors (X,Y,Z) CTS */
} SLC_PAR;

typedef struct {			/* ESA Wide Swath SLC burst parameters */
  /*
    ASAR WS beam illumination order:       1-3-5-2-4
  */
  int bn;		/* WSS beam number (1-->5) */
  int rank;		/* number of pulses in the air */
  int lpb;		/* lines per burst (echoes received for each burst) total time in mode = (rank + lines/burst)/PRF */
  double tdiff;		/* time difference between sensing time of the first input line and the 
			zero Doppler time of the first output image line */
  double tfss;		/* on-orbit sensing time (s) of the first echo source packet for this beam */
  double cycle;		/* ASAR SCANSAR cycle time (s)*/
} WSS_PAR;

/* ASAR WSS product Doppler Grid ADS */
typedef struct {
  double az_time;	/* zero doppler time */
  double dop[100];	/* Doppler centroid estimates Hz. */
  double rng[100];	/* slant range (m) */
} WSS_DOP_ADS;

/* GPRI data structure 20080825 clw */
typedef struct {  
  double az_start_angle;/* starting rotation angle (degrees, + is clockwise looking down the rotation axis) */ 
  double az_angle_step;	/* angular step degrees (+angle is clockwise rotation looking down the rotation axis) */
  double ant_elev_angle;/* antenna elevation angle, + angles are up,*/
  double ref_north;	/* reference point northing or latitude in the projection and datum of the DEM */
  double ref_east;	/* reference point easting or longitude in the projection and datum of the DEM */
  double ref_alt;	/* reference point altitude in the projection and datum of the DEM */
  double scan_heading;	/* heading of central sweep scan line relative to north, + is clockwise from north, 
			   looking down the tower rotation axis. The tower rotation axis is aligned with vertical down towards the earth center. Down is one component of the North, East, Down (NED) coordinate system */
  VEC tx_coord;		/* transmit antenna phase center coordinates (XYZ meters) in the local tower coordinate system */
  VEC rx1_coord;	/* receive antenna 1 phase center coordinates (XYZ meters) in the local tower coordinate system */
  VEC rx2_coord;	/* receive antenna 2 phase center coordinates (XYZ meters) in the local tower coordinate system */
  double tower_roll;	/* tower roll angle rotation about the local Y axis */
  double tower_pitch;	/* tower pitch angle rotation about the local X axis */
  double phase_offset;	/* interferogram phase offset s.t. ph = -4pi/lam* (r2 - r1) + phase_offset */
}GPRI_PAR;

typedef struct {		/* interferogram parameters and SLC offset data */
  char title[128];		/* interferogram title */
  int initr, initaz;		/* integer offsets in pixels, range and azimuth */
  int nofstr, npr;		/* number of points offset from start, number of points range to process */
  int rstr, rend, nr, rsp;	/* starting range sample, ending range sample, number of points in range for offset estimates, range offset sample spacing */
  int azstr, azend, naz, azsp;	/* starting azimuth, ending azimuth, number of azimuth points for offset estimates, azimuth offset sample spacing */
  int rwin, azwin;		/* maximum sizes of window in range, azimuth to estimate offset at each point */
  double thres;			/* fringe SNR threshold to save offset estimate or discard */
  double rpoly[MAR];		/* range offset polynomial as a function of range and azimuth; range coordinate is relative to
				   interferogram starting range pixel (nofstr), azimuth coordinate is relative to SLC-1 */
  double azpoly[MAAZ];		/* azimuth offset polynomial as a function of range in azimuth; range coordinate is relative to
				   interferogram starting range pixel (nofstr), azimuth coordinate is relative to SLC-1 */
  int lbegin;			/* offset to starting line in interferogram (relative to start of SLC-1) */
  int nls;			/* number of interfogram lines */
  int nr_int;			/* width of interferogram (samples) */
  int nrb;			/* offset from start of each line for the first valid interferogram pixel */
  int nrps;			/* number of valid slant range pixels in the interferogram */
  int rlks, azlks;		/* number of interferometric looks in range and azimuth */
  double rps_int, azps_int; 	/* interferogram sample spacing (meters) in slant range and azimuth */
  double rps_res, azps_res;  	/* true ground range and azimuth spacing (meters) in resampled height map */
  double grg_start;		/* starting ground range relative to center swath at altitude=0.0 */
  int ngrg;			/* samples in the cross-track direction */
  int ngaz;			/* samples in the azimuth direction */
} OFF_PAR;

typedef struct {
  /*****************************************************************************************
      	State vector estimate of baseline at center of SLC-1,
  	baseline velocity, (expressed in TCN) (m)
  ******************************************************************************************/
  VEC_TCN base_orb;
  VEC_TCN bdot_orb;
  /*****************************************************************************************
        Baseline at center of SLC-1, baseline velocity, (expressed in TCN) (m)
  	using L.S. fitting of the baseline parameters from GCPs
  *******************************************************************************************/
  VEC_TCN base;
  VEC_TCN bdot;
  /******************************************************************************************/
  double phc;	/* interferometric phase constant estimated using Least Squares model fitting
		of the baseline from GCPs (radians) */
} BASELINE;

typedef struct {
  /*****************************************************************************************
      	State vector or offset derived estimate of baseline at center of SLC-1 in SCH,
  	baseline rates/meter expresed in SCH (m)
  ******************************************************************************************/
  VEC_SCH base_est;
  VEC_SCH dbds_est;	/* derivative with respect to distance dbase/ds */
  double phc_est;	/* estimated phase constant */
  /*****************************************************************************************
        Baseline at center of SLC-1, baseline derivates with respect to position (SCH) (m)
  	using L.S. fitting of the baseline parameters from GCPs
  *******************************************************************************************/
  VEC_SCH base;
  VEC_SCH dbds;
  /******************************************************************************************/
  double phc;	/* interferometric phase constant estimated using Least Squares model fitting
		of the baseline from GCPs (radians) */
} BASE_SCH;

typedef struct {
  int index;	/* GCP index label */
  int irg;	/* slant range pixel position */
  int iaz;	/* azimuth pixel position */
  double alt;	/* altitude for point relative to reference ellipsoid (m)*/
  double ph;	/* unwrapped measured phase (radians) */
  double phr;	/* phase with range phase added in (radians) */
  double lat;	/* latitude of point (deg.) */
  double lng;	/* longitude of point (deg.) */
  double ti;	/* time relative to center of scene (s) */
  double rg;	/* slant range (m) */
  double th;	/* look angle (relative to TCN) (radians) */
  VEC_TCN base;	/* baseline expressed in TCN */
  VEC_TCN lv;	/* look-vector for GCP expressed in TCN */
} GCP;
