#include "display.h"
#include "bmp_image.h"
#include "rasterfile.h"
#define REAL float

#define MAX_STR_LEN 132
#define PROB "min"		/* only work with minimum problems */

#define G_HUGE_UINT  4294967295
#define G_HUGE_INT   2147483647
#define G_SHORT_MAX  32767
#define G_EPSILON	1.e-7
#define G_NULL  -1

#define UNLABELED	0
#define LABELED		1
#define PERM_LABELED	2
#define MAX_PATH_NODES  500000	/* maximum number of nodes along the shortest path */
#define MAX_NODE_ARCS   3	/* maximum number of arcs that can leave or enter a node */

#define EDGE_PIX	4	/* used to determine cost of connecting to edge cost = COST_SCALE_FACTOR*EDGE_PIX*Average correlation */
#define NULL_COST	50	/* arc cost when the weight is 0.0 */
#define EDGE_FRAC	0.1	/* residues within this image fraction to the edge can connect through the super node */
#define COST_SCALE_FACTOR 100
#define CAP_MAX      128	/* maximum arc capacity */

#define MAX_TRIANGLE_LIST_SIZE  30000	/* length of ping-pong triangle lists */

#define RANDOM 1 		/* random selection of excess nodes */
#define G_DEPTH 8		/* depth to grow when searching edge */
#define MAX_EDGE_STEP 32	/* maximum spacing around the hull for arcs */

#ifdef GTS
#include <gts.h>		/* use open source GTS library for triangulation, requires glib from www.gtk.org */
#else
#include "triangle.h"		/* use triangle routine */
#endif

typedef struct {
  unsigned short x, y;
}
G_POS;

typedef struct {
  short b;		/* supply b(i) */
  short e;		/* excess e(i) */
  int pi;		/* node potential PI */
  int pred;		/* predecessor node */
  int d;		/* total cost to this node */
  int hdx;		/* heap index for this node */
  char state;		/* Dijkstra node state, UNLABELED, LABELED, PERM_LABELED */
  unsigned short narc_out;	/* number of outward directed arcs at this node */
  int *arc_out;			/* array containing indices of arcs leaving this node */
  unsigned short narc_in; 	/* number of inward directed arcs at this node */
  int *arc_in;			/* array containing indices of arcs entering this node */
}
G_NODE;

typedef struct {
  int head;		/* head of the arc, node where it is headed */
  int tail;		/* arc tail, node where comes from */
  short flow;		/* flow in this arc */
  /*  short low; 	*/	/* lower bound for arc flow */
  short cap;		/* upper bound on arc capacity */
  short rcap;		/* residual capacity */
  short cost;		/* arc cost per unit flow */
  short rcost;		/* reduced cost of this arc per unit flow */
}
G_ARC;

typedef struct {
  REAL x, y;
}
POINT;

typedef struct {		/* corners of the triangle are indices into the array of points */
  int c[3];
}
TRIANGLE;

typedef struct {		/* triangle neighbors are indices into the array of triangles */
  int t[3];
}
NEIGHBORS;

typedef struct {
  int x, y;
}
COORD;

typedef struct {
  int index;		/* index into the list of triangles */
  COORD crn[3];		/* pairs of x,y coordinates */
  int nbt[3];		/* neighboring triangles crn[0]<->cnr[1] edge has triangle neighbor nb[0],
    			   crn[1]<->crn[2] has triangle neighbor nv[1]... */
  unsigned char flg[3];	/* flag if the point is unwrapped or not, 0=not unwrapped, 1=unwrapped */
}
TRI;

extern void start_timing();
extern void stop_timing();

extern unsigned char *uvector(int nl, int nh);
extern unsigned char **umatrix(int nrl, int nrh, int ncl, int nch);

extern void free_uvector(int *v, int nl, int nh);
extern void free_umatrix(unsigned char **m, int nrl, int nrh, int ncl, int nch);

/* routines to read Sun raster and BMP format images */
extern int rd_image(char *fn, unsigned char *image, int im_format, COLOR *ct, int *nc, int *nbits, int *width, int *height);
extern int init_rdras(char *f1, FILE **f1f, RASTER_HDR **r1, COLOR **ct, int *ncolors, int *nbits, int *width, int *height);
extern int init_rdbmp(char *f_ras, FILE **rf, BITMAPFILEHEADER **bfh, BITMAPINFOHEADER **bih, COLOR **ct, int *nc, int *nbits, int *width, int *height);
extern int rdras(FILE *f1f, RASTER_HDR *r1, int iline, unsigned char *data);
extern int rdbmp(FILE *rf, BITMAPFILEHEADER *bfh, BITMAPINFOHEADER *bih, int iline, unsigned char *data);
