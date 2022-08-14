/*
   unwrap phase using triangulation and MCF (successive shortest paths)
   clw Copyright Gamma Remote Sensing AG 6-Aug-2008
*/
#include "mcf2.h"
#define PATCH_OVERLAP	512		/* nominal overlap of patches */
#define PATCH_R_SIZE	1024		/* nominal range patch size */
#define PATCH_AZ_SIZE	1024		/* nominal azimuth patch size */
#define NL_PAMB		7		/* number of lines used to determine phase ambiguity number (odd) */

extern void m_ssp2(G_NODE *nodes, G_ARC *arcs, int net_node, int net_arc);
extern void m_tri_init2(POINT *pt, TRIANGLE *tr, NEIGHBORS *nb, float *wn, short *chrg, int numberoftriangles,
                          int net_node, int net_arc, G_NODE *nodes, G_ARC *arcs, int dwgt_mode);

extern void mssptr(float **a, float **wgt, unsigned char **msk, float **unw, int naz, int nr, int mflg, int tmode);
extern void **calloc_2d(size_t ncols, size_t nlines, size_t size);	/* 2D memory array allocation and initializaton */
extern void *calloc_1d(size_t nv, size_t size);			/* 1D memory array allocation and initializaton */
extern void zero_2d(void **a, int width, int nlines, size_t size);	/* 2D memory initialization */
extern void zero_1d(void *a, int width, size_t size);			/* 1D memory initialization */
extern void wrp_flt(float **tunw, int ip, int nr, int naz, FILE *fltf);

int rdp_cpx_phase(float **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *cpxf);
int rdp_cpx(fcomplex **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *fltf);
int rdp_float(float **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *fltf);
int rdp_msk_ras(unsigned char **msk, int width, int nls, int roff, int azoff, int nr, int naz,
                RASTER_HDR *r1, COLOR *ct1, int nc1, FILE *mskf);
int rdp_msk_bmp(unsigned char **msk, int width, int nls, int roff, int azoff, int nr, int naz,
                BITMAPFILEHEADER *bfh1, BITMAPINFOHEADER *bih1, COLOR *ct1, int nc1, FILE *mskf);

int main(int argc, char **argv) {

  double ph0, ph1;	/* wrapped and unwrapped phases at phase reference point */

  int i, j, l, k;		/* loop counters */
  int i1, i2;		/* azimuth offset for overlapping patches */
  int j1, j2;		/* range offset for overlapping patches */
  int ip;		/* patch counter */
  int width;		/* width of interferogram data file */
  int wflg = 0;		/* weight flag */
  int mflg = 0;		/* mask flag */
  int nr1, naz1;		/* number of range and azimuth pixels to process */
  int roff = 0, azoff = 0;	/* range and azimuth offsets for unwrapping */
  int roff1, azoff1;	/* range and azimuth offsets to a patch */
  int roff2, azoff2;
  int nls;		/* number of lines in the interferogram file */
  int mw, mhgt;		/* mask width and height */
  int nbits1, nc1;	/* mask raster file parameters */
  int im_flag1 = UNKNOWN;
  int rinit, azinit;	/* phase reference point coordinates */
  int pflg2 = 0;		/* phase flag determines if the starting phase is set to 0, 0(default): use initial value, 1: set to 0.0 */
  int tmode = 0;		/* triangulation mode 0: filled triangular mesh   1: Delaunay triangular mesh */
  int npat_r = 0;		/* number of range patches (across) */
  int npat_az = 0;	/* number of azimuth patches (down) */
  int nr;		/* number of range samples in the patches */
  int naz;		/* number of azimuth samples in the patches */
  int ovrlp = PATCH_OVERLAP;	/* patch overlap */
  int ovrlp2, ovrlp4;
  int naz_rd;		/* number of lines read to a patch */
  int nls_wgt;		/* number of lines in the weight file */
  int nl1, nl2;		/* starting and ending lines in the patch to write out */
  int ns1;		/* number of output segments */
  int minm, maxm;	/* minimum and maximum multiples */
  int mc;		/* number of elements in the list */
  int nh;		/* number of histogram values */
  int mx_amb;		/* phase unwrapping ambiguity determined from overlap data */
  int nl_pamb = NL_PAMB;	/* number of lines to check for phase ambiguity */

  fcomplex **a3, *ax3;

  float *ftmp;		/* single line buffer for floating point data */
  float **a, *a1;	/* phase modulo 2PI before unwrapping 2D array */
  float **u1;		/* unwrapped phase 2D array */
  float **u2;
  float **u3, *unw3;	/* output unwrapped phase */
  float **wgt, *wgt1;	/* weights */
  float **w3 = NULL, *wgt3 = NULL; /* weight data */
  float *mwgt;		/* array containing weights of TWO_PI in the overlap */
  float *hist;		/* histogram of multiples of 2PI */

  int *m2pi;		/* array containing multiples of TWO_PI in the overlap */

  unsigned char **m1, *msk1;	/* mask 2D array */
  unsigned char **m2, *msk2;
  unsigned char **m3 = NULL, *msk3 = NULL;
  unsigned char *uctmp;		/* array for reading 8-bit image data for mask */

#ifdef fread
  unsigned char btmp, *b1;	/* used for byte swapping */
#endif

  char bmp_ext[] = "bmp\0";
  char ras_ext[] = "ras\0";
  char mt1[255];

  RASTER_HDR *r1;
  BITMAPFILEHEADER *bfh1;
  BITMAPINFOHEADER *bih1;
  COLOR *ct1;

  FILE *intf = NULL, *wgtf = NULL, *mskf = NULL, *unwf = NULL;

#ifdef DEBUG_PATCH
  FILE *utmpf;
#endif

  printf("*** Phase unwrapping using Minimum Cost Flow (MCF) and triangulation ***\n");
  printf("*** Copyright 2008, Gamma Remote Sensing, v1.7 clw/uw 6-Aug-2008 ***\n");

  if (argc < 6) {
    printf("\nusage: %s <interf> <wgt> <mask> <unw> <width> [tri_mode] [roff] [loff] [nr] [nlines] [npat_r] [npat_az] [ovrlap] [r_init] [az_init] [init_flag]\n\n", argv[0]) ;
    printf("input parameters:\n");
    printf("  interf     (input) interferogram (*.int,*.flt)(fcomplex)\n");
    printf("  wgt        (input) weight factors (0 -> 1.0) file (float)(enter - for uniform weight)\n");
    printf("  mask       (input) validity mask (SUN raster or BMP format, value 0 -> pixel not used)(enter - if no mask)\n");
    printf("  unw        (output) unwrapped phase image (*.unw)(float)\n");
    printf("  width      number of samples/row\n");
    printf("  tri_mode   triangulation mode\n");
    printf("               0: filled triangular mesh (default)\n");
    printf("               1: Delaunay triangulation\n");
    printf("  roff       offset to starting range of section to unwrap (default: 0)\n");
    printf("  loff       offset to starting line of section to unwrap (default: 0)\n");
    printf("  nr         number of range samples of section to unwrap (default(-): width - roff)\n");
    printf("  nlines     number of lines of section to unwrap (default(-): total number of lines - loff)\n");
    printf("  npat_r     number of patches in range\n");
    printf("  npat_az    number of patches in azimuth\n");
    printf("  ovrlap     overlap between patches in pixels (overlap >= %d, default(-): %d)\n", nl_pamb, ovrlp);
    printf("  r_init     phase reference point range offset (default(-): roff)\n");
    printf("  az_init    phase reference point azimuth offset (default(-): loff)\n");
    printf("  init_flag  flag to set phase at reference point\n");
    printf("               0: use initial point phase value (default)\n");
    printf("               1: set phase to 0.0 at initial point\n\n");
    exit(-1);
  }

  start_timing();
  intf = fopen(argv[1], "r+b");
  if (intf == NULL) {
    fprintf(stderr, "\nERROR: cannot open complex interferogram file: %s\n\n", argv[1]);
    exit(-1);
  }
  printf("input interferogram file: %s\n", argv[1]);

  if (strcmp(argv[2], "-") != 0) {
    wgtf = fopen(argv[2], "r+b");
    if (wgtf == NULL) {
      fprintf(stderr, "\nERROR: cannot open weights file: %s\n\n", argv[2]);
      exit(-1);
    }
    printf("weight file: %s\n", argv[2]);
    wflg = 1;
  }

  if (strcmp(argv[3], "-") != 0) {
    mskf = fopen(argv[3], "r+b");
    if (mskf == NULL) {
      fprintf(stderr, "\nERROR: cannot open mask file (Sun raster/BMP format) file: %s\n\n", argv[3]);
      exit(-1);
    }
    mflg = 1;
    printf("mask file: %s\n", argv[3]);
  }

  unwf = fopen(argv[4], "w+b");
  if (unwf == NULL) {
    fprintf(stderr, "\nERROR: cannot open output unwrapped phase file: %s\n\n", argv[4]);
    exit(-1);
  }
  printf("unwrapped phase file: %s\n", argv[4]);

#ifdef DEBUG_PATCH
  utmpf = fopen("temp.unw", "w+b");
  if (utmpf == NULL) {
    fprintf(stderr, "\nERROR: cannot open temporary unwrapped phase file\n\n");
    exit(-1);
  }
  printf("temporary unwrapped phase file: %s\n", "temp.unw");
#endif

  sscanf(argv[5], "%d", &width);				/* width of interferogram file */
  printf("\nline width (samples): %d\n", width);
  fseek(intf, (off_t)0L, SEEK_END);
  nls = (int)(ftell(intf) / (off_t)(width * 2 * sizeof(float)));
  rewind(intf);
  printf("number of lines in the interferogram: %d\n", nls);

  if (wflg == 1) {
    fseek(wgtf, (off_t)0L, SEEK_END);
    nls_wgt = (int)(ftell(wgtf) / (off_t)(width * sizeof(float)));
    rewind(wgtf);
    printf("number of lines in the weight file %d\n", nls_wgt);
  }

  if (argc > 6) {
    sscanf(argv[6], "%d", &tmode);
    if (tmode < 0 || tmode > 1) {
      fprintf(stderr, "\nERROR: invalid triangulation mode: %d\n\n", tmode);
      exit(-1);
    }
    switch (tmode) {
    case 0:
      printf("triangulation mode: filled triangular mesh\n");
      break;
    case 1:
      printf("triangulation mode: Delaunay triangulation\n");
      break;
    default:
      fprintf(stderr, "\nERROR: invalid triangulation mode: %d\n\n", tmode);
      exit(-1);
    }
  }

  if (argc > 7) {
    sscanf(argv[7], "%d", &roff);
    if (roff < 0) {
      fprintf(stderr, "\nERROR: range offset < 0 : %d\n\n", roff);
      exit(-1);
    }
    if (roff >= width) {
      fprintf(stderr, "\nERROR: range offset >= width of the interferogram: %d\n\n", roff);
      exit(-1);
    }
  }

  if (argc > 8) {
    sscanf(argv[8], "%d", &azoff);
    if (azoff < 0) {
      fprintf(stderr, "\nERROR: azimuth offset < 0: %d\n\n", azoff);
      exit(-1);
    }
    if (azoff >= nls) {
      fprintf(stderr, "\nERROR: azimuth offset >= number of lines in the interferogram: %d\n\n", azoff);
      exit(-1);
    }
  }

  nr1 = width - roff;
  if (argc > 9) {
    if ((strcmp(argv[9], "-") != 0) && (strcmp(argv[9], "0") != 0)) {
      sscanf(argv[9], "%d", &nr1);
      if (nr1 <= 0) {
        fprintf(stderr, "\nERROR: number of range samples <= 0: %d\n\n", nr1);
        exit(-1);
      }
      if (nr1 > width - roff) {
        fprintf(stderr, "\nERROR: insufficent samples for desired interferogram segment width: %d\n\n", nr1);
        exit(-1);
      }
    }
  }

  naz1 = nls - azoff;
  if (argc > 10) {
    if ((strcmp(argv[10], "-") != 0) && (strcmp(argv[10], "0") != 0)) {
      sscanf(argv[10], "%d", &naz1);
      if (naz1 <= 0) {
        fprintf(stderr, "\nERROR: number of lines <= 0: %d\n\n", naz1);
        exit(-1);
      }
      if (naz1 > nls - azoff) {
        fprintf(stderr, "\nERROR: insufficent lines for desired interferogram segment height: %d\n\n", naz1);
        exit(-1);
      }
    }
  }

  if (argc > 11) {
    if (strcmp(argv[11], "-") != 0) {
      sscanf(argv[11], "%d", &npat_r);
      if (npat_r <= 0) {
        fprintf(stderr, "\nERROR:number of range patches <= 0: %d\n\n", npat_r);
        exit(-1);
      }
    }
  }

  if (argc > 12) {
    if (strcmp(argv[12], "-") != 0) {
      sscanf(argv[12], "%d", &npat_az);
      if (npat_az <= 0) {
        fprintf(stderr, "\nERROR:number of azimuth patches <= 0: %d\n\n", npat_az);
        exit(-1);
      }
    }
  }

  if (argc > 13) {
    if (strcmp(argv[13], "-") != 0) {
      sscanf(argv[13], "%d", &ovrlp);
      if (ovrlp < 7) {
        fprintf(stderr, "\nERROR: overlap between patches less than minimum required value: %d\n\n", nl_pamb);
        exit(-1);
      }
    }
  }

  if (npat_r == 0) {
    if (ovrlp > PATCH_R_SIZE) {
      fprintf(stderr, "\nERROR: enter patch overlap that is smaller than than the range patch width: %d\n\n", PATCH_R_SIZE);
      exit(-1);
    }
    npat_r = MAX(ceil((double)(nr1 - ovrlp) / (PATCH_R_SIZE - ovrlp)), 1);
    nr = ceil((double)(nr1 + (npat_r - 1) * ovrlp) / npat_r);
  } else {
    nr = ceil((double)(nr1 + (npat_r - 1) * ovrlp) / npat_r);
  }

  if (npat_az == 0) {
    if (ovrlp > PATCH_AZ_SIZE) {
      fprintf(stderr, "\nERROR: enter patch overlap that is smaller than than the azimuth patch length: %d\n\n", PATCH_AZ_SIZE);
      exit(-1);
    }
    npat_az = MAX(ceil((double)(naz1 - ovrlp) / (PATCH_AZ_SIZE - ovrlp)), 1);
    naz = ceil((double)(naz1 + (npat_az - 1) * ovrlp) / npat_az);
  } else {
    naz = ceil((double)(naz1 + (npat_az - 1) * ovrlp) / npat_az);
  }

  printf("\nnumber of range samples to unwrap: %d  number of azimuth lines to unwrap: %d\n", nr1, naz1);
  printf("number of range patches:   %5d   range patch width:   %5d\n", npat_r, nr);
  printf("number of azimuth patches: %5d   azimuth patch width: %5d\n", npat_az, naz);
  printf("patch overlap: %d\n", ovrlp);
  ovrlp2 = ovrlp / 2;
  ovrlp4 = ovrlp / 4;

  rinit = roff;
  if (argc > 14) {
    if (strcmp(argv[14], "-") != 0) {
      sscanf(argv[14], "%d", &rinit);
      if ((rinit  < roff) || (rinit >= (roff + nr1))) {
        fprintf(stderr, "\nERROR: range position of phase reference point outside of image segment\n");
        exit(-1);
      }
    }
  }

  azinit = azoff;
  if (argc > 15) {
    if (strcmp(argv[15], "-") != 0) {
      sscanf(argv[15], "%d", &azinit);
      if ((azinit < azoff) || (azinit >= (azoff + naz1))) {
        fprintf(stderr, "\nERROR: azimuth position of phase reference point outside of image segment\n");
        exit(-1);
      }
    }
  }
  printf("phase reference point offset range: %d  azimuth: %d\n", rinit, azinit);

  if (argc > 16) {
    sscanf(argv[16], "%d", &pflg2);
    switch (pflg2) {
    case 0:
      printf("using initial point value as starting phase\n");
      break;
    case 1:
      printf("setting intial phase value to 0.0\n");
      break;
    default:
      fprintf(stderr, "\nERROR: invalid value for init phase flag: %d\n\n", pflg2);
      exit(-1);
    }
  }

  if (mflg == 1) {			/* open Sun raster or BMP format mask file */
    strcpy(mt1, argv[3]);
    if (strstr(mt1, ras_ext) != NULL)im_flag1 = SUN;
    if (strstr(mt1, bmp_ext) != NULL)im_flag1 = BMP;
    if (im_flag1 == UNKNOWN) {
      fprintf(stderr, "\nERROR: only SUN raster and BMP formats are supported: %s\n\n", mt1);
      exit(-1);
    }
    if (im_flag1 == SUN) init_rdras(mt1, &mskf, &r1, &ct1, &nc1, &nbits1, &mw, &mhgt);
    else init_rdbmp(mt1, &mskf, &bfh1, &bih1, &ct1, &nc1, &nbits1, &mw, &mhgt);
    if (nbits1 != 8)fprintf(stderr, "\nERROR: mask file can only be 8 bits/pixel\n");
    printf("\ninput phase unwrapping mask file: %s\n", mt1);
    printf("input phase unwrapping width: %5d   height: %5d\n", mw, mhgt);
    if (mw != width) {
      fprintf(stderr, "\nERROR: width of mask and input interferogram do not match: %d %d\n\n", mw, width);
      exit(-1);
    }
    if (mhgt != nls) {
      fprintf(stderr, "\nERROR: height of mask and input interferogram do not match: %d %d\n\n", mhgt, nls);
      exit(-1);
    }
  }

  ftmp = (float *)malloc(width * sizeof(float));
  uctmp = (unsigned char *)malloc(width + 3);

  a1 = (float *)malloc(nr * naz * sizeof(float));			/* allocate memory and initialize arrays */
  a =  (float **)malloc(sizeof(fcomplex *) * naz);

  u1 = (float **)calloc_2d(nr, naz, sizeof(float));		/* 2D memory array allocation and initializaton */
  u2 = (float **)calloc_2d(nr, naz, sizeof(float));		/* 2D memory array allocation and initializaton */

  wgt1 = (float *)malloc(nr * naz * sizeof(float));
  wgt =  (float **)malloc(sizeof(float *) * naz);

  msk1 = (unsigned char *)malloc(nr * naz);
  m1 = (unsigned char **)malloc(sizeof(unsigned char *) * naz);

  msk2 = (unsigned char *)malloc(nr * naz);
  m2 = (unsigned char **)malloc(sizeof(unsigned char *) * naz);

  m2pi = (int *)calloc(sizeof(int), ovrlp * MAX(naz, nr));		/* arrays for resolving phase ambiguity between patches */
  mwgt = (float *)calloc(sizeof(float), ovrlp * MAX(nr, naz));

  if (a1 == NULL || a == NULL || wgt1 == NULL || wgt == NULL || uctmp == NULL || m2pi == NULL || mwgt == NULL ||
      msk1 == NULL || m1 == NULL || ftmp == NULL || msk2 == NULL || m2 == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for 2D arrays\n\n");
    exit(-1);
  }

  unw3 = (float *)calloc(sizeof(float), naz * width);		/* output phase */
  u3 = (float **)malloc(sizeof(float *) * naz);
  if (unw3 == NULL || u3 == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for output unwrapped phase\n\n");
    exit(-1);
  }
  printf("memory allocation for output unwrapped phase buffer: %10.3f Mb\n", (sizeof(float)*naz*width) / 1048576.);

  for (i = 0; i < naz; i++) {		/* assign row pointers */
    a[i]  = a1 + nr * i;
    u3[i] = unw3 + i * width;
    m1[i] = msk1 + nr * i;
    m2[i] = msk2 + nr * i;
    wgt[i] = wgt1 + nr * i;
  }

  printf("\nimage segment range offset:  %8d  azimuth offset: %8d\n", roff, azoff);
  printf("image segment range samples: %8d  azimuth lines:  %8d\n", nr, naz);

  if (wflg == 0) {			/* no weight file, set all weights to 1.0 */
    for (i = 0; i < naz; i++) {
      for (j = 0; j < nr; j++) wgt[i][j] = 1.0;
    }
  }

  roff1 = roff;
  azoff1 = azoff;

  if (azoff != 0) {				/* write blank lines if necessary */
    printf("writing %d blank lines at the beginning of the output file\n", azoff);
    for (j = 0; j < azoff; j++) {
      fwrite((char *)u3[0], sizeof(float), width, unwf);
    }
  }

  for (l = 0; l < npat_az; l++) {
    for (k = 0; k < npat_r; k++) {
      ip = l * npat_r + k;
      printf("\n***** PATCH: %d *****\n", ip + 1);
      roff1 = roff + k * (nr - ovrlp);
      azoff1 = azoff + l * (naz - ovrlp);

      naz_rd = rdp_cpx_phase(a, width, nls, roff1, azoff1, nr, naz, intf);

      printf("phase patch  (r,az): %5d  %5d   r_offset:%5d  az_offset:%5d  lines read:%5d\n", k + 1, l + 1, roff1, azoff1, naz_rd);
      if (naz_rd <= 0) {
        fprintf(stderr, "\nERROR: no data read from complex interferogram file\n");
        exit(-1);
      }
      if (wflg == 1) {
        naz_rd = rdp_float(wgt, width, nls, roff1, azoff1, nr, naz, wgtf);
        printf("weight patch (r,az): %5d  %5d   r_offset:%5d  az_offset:%5d  lines read:%5d\n", k + 1, l + 1, roff1, azoff1, naz_rd);
      }

      if (mflg == 1) {						/* read raster image mask file */
        if (im_flag1 == SUN)naz_rd = rdp_msk_ras(m1, mw, mhgt, roff1, azoff1, nr, naz, r1, ct1, nc1, mskf);
        else naz_rd = rdp_msk_bmp(m1, mw, mhgt, roff1, azoff1, nr, naz, bfh1, bih1, ct1, nc1, mskf);
        printf("mask patch   (r,az): %5d  %5d   r_offset:%5d  az_offset:%5d  lines read:%5d\n", k + 1, l + 1, roff1, azoff1, naz_rd);
      }

      zero_2d((void **)u1, nr, naz, sizeof(float));		/* initialize output unwrapped phase to 0.0 */
      mssptr(a, wgt, m1, u1, naz, nr, mflg, tmode); 		/* phase unwrapping with MCF SSP */
      ph0 = 0.0;						/* default value for phase integration constant */

#ifdef DEBUG_PATCH
      wrp_flt(u1, ip, nr, naz, utmpf);
      * / 			/* debugging of determining of the phase unwrapping constant */
#endif

      if (ip == 0)goto copy_patch;				/* for the first patch, do not adjust the phase */

      if ((l == 0) || (k != 0)) {				/* for the first row read adjacent patch */
        roff2 = roff + (k - 1) * (nr - ovrlp);
        azoff2 = azoff1;
        printf("\nreading patch adjacent to current patch: %d  roff2: %d  azoff2: %d\n", ip, roff2, azoff2);

        for (i = 0; i < naz; i++) {
          for (j = 0; j < nr; j++)u2[i][j] = u3[i][j + roff2];
        }

        if (mflg == 1) {
          if (im_flag1 == SUN)naz_rd = rdp_msk_ras(m2, mw, mhgt, roff2, azoff2, nr, naz, r1, ct1, nc1, mskf);
          else naz_rd = rdp_msk_bmp(m2, mw, mhgt, roff2, azoff2, nr, naz, bfh1, bih1, ct1, nc1, mskf);
          printf("adjacent patch mask (r,az): %5d  %5d   r_offset2:%5d  az_offset2:%5d  lines read:%5d\n", k - 1, l, roff2, azoff2, naz_rd);
        }

        mc = 0; 	      /* counter of valid points in overlap */
        for (i = ovrlp2; i < naz - ovrlp2; i++) {
          for (j = -nl_pamb / 2; j < nl_pamb / 2; j++) {
            j1 = (ovrlp - ovrlp2) + j;
            j2 = nr - ovrlp / 2 + j;

#ifdef DEBUG_AMB
            if (i % 20 == 0)printf("i,j,j1,j2,u1[i][j1],u2[i][j2]: %5d %5d %5d %5d %10.3f %10.3f\n", i, j, j1, j2, u1[i][j1], u2[i][j2]);
#endif

            if ((mflg == 1) && ((m1[i][j1] == 0) || (m2[i][j2] == 0))) continue;	/* if point masked out, do not add to list */
            if (a[i][j1] != 0.0) {
              m2pi[mc] = nint((u1[i][j1] - u2[i][j2]) / TWO_PI);
              mwgt[mc] = wgt[i][j1];
              mc++;
            }
          }
        }
      } else {   							/* derive phase constant from unwrapped phase in the row above  */
        roff2 = roff + k * (nr - ovrlp);
        azoff2 = azoff + (l - 1) * (naz - ovrlp);

        printf("\nusing patch above current patch: %d  roff2: %d  azoff2: %d\n", ip - npat_r + 1, roff2, azoff2);

        for (i = 0; i < naz; i++) {
          for (j = 0; j < nr; j++)u2[i][j] = u3[i][j + roff2];	/* u3 array contains results from previous row */
        }

        if (mflg == 1) {
          if (im_flag1 == SUN)naz_rd = rdp_msk_ras(m2, mw, mhgt, roff2, azoff2, nr, naz, r1, ct1, nc1, mskf);
          else naz_rd = rdp_msk_bmp(m2, mw, mhgt, roff2, azoff2, nr, naz, bfh1, bih1, ct1, nc1, mskf);
          printf("above patch mask (r,az): %5d  %5d   r_offset:%5d  az_offset:%5d  lines read:%5d\n", k, l - 1, roff2, azoff2, naz_rd);
        }

        mc = 0; 	      /* counter of valid points in overlap */
        for (i = -nl_pamb / 2; i < nl_pamb / 2; i++) {
          i1 = (ovrlp - ovrlp2) + i;
          i2 = naz - ovrlp / 2 + i;
          for (j = ovrlp2; j < nr - ovrlp2; j++) {

#ifdef DEBUG_AMB
            if (j % 20 == 0)printf("i,j,i1,i2,u1[i1][j],u2[i2][j]: %5d %5d %5d %5d %10.3f %10.3f\n", i, j, i1, i2, u1[i1][j], u2[i2][j]);
#endif

            if ((mflg == 1) && ((m1[i1][j] == 0) || (m2[i2][j] == 0))) continue;	/* if point masked out, do not add to list */
            if (a[i1][j] != 0.0) {
              m2pi[mc] = nint((u1[i1][j] - u2[i2][j]) / TWO_PI);
              mwgt[mc] = wgt[i1][j];
              mc++;
            }
          }
        }
      }

      printf("number of elements in overlap region: %d\n", mc);
      if (mc == 0) {
        printf("WARNING: no valid samples in the overlap region for estimation of the phase ambiguity\n\n");
        goto copy_patch;
      }

      minm = m2pi[0];
      maxm = m2pi[0];

      for (i = 1; i < mc; i++) {
        if (minm > m2pi[i])minm = m2pi[i];
        if (maxm < m2pi[i])maxm = m2pi[i];
      }
      printf("minimum ambiguity multiples of 2PI: %d\n", minm);
      printf("maximum ambiguity multiples of 2PI: %d\n\n", maxm);

      nh = maxm - minm + 1;
      hist = calloc(sizeof(float), nh);

      for (i = 0; i < mc; i++) {
        hist[m2pi[i] - minm] += mwgt[i];
      }

      mx_amb = 0;
      for (i = 0; i < nh; i++) {
        if (hist[i] > hist[mx_amb])mx_amb = i;
        printf("ambiguity: %4d   sum weights: %12.3f\n", i + minm, hist[i]);
      }

      mx_amb += minm;
      ph0 = mx_amb * TWO_PI;
      printf("selected 2PI ambiguity: %d   phase offset: %10.3f radians\n\n", mx_amb, ph0);
      free(hist);

copy_patch:
      if (k == 0)j1 = 0;					/* first patch in the row uses all the range samples at the start */
      else j1 = ovrlp - ovrlp2;

      if (k == npat_r - 1) j2 = roff + nr1 - roff1;	/* test edge on last patch in the row */
      else j2 = nr;

      printf("\npatch starting and ending range sample to copy: %d %d\n", j1, j2 );
      if (mflg == 1) {					/* check if mask set to 0, then output = 0.0 */
        for (i = 0; i < naz; i++) {
          for (j = j1; j < j2; j++) {
            if (m1[i][j] != 0) u3[i][j + roff1] = u1[i][j] - ph0;
            else u3[i][j + roff1] = 0.0;
          }
        }
      } else {
        for (i = 0; i < naz; i++) {
          for (j = j1; j < j2; j++) {
            u3[i][j + roff1] = u1[i][j] - ph0;
          }
        }
      }
    }

    if (l == 0)nl1 = 0;					/* test for first row */
    else nl1 = ovrlp - ovrlp2;

    if (l == npat_az - 1)nl2 = azoff + naz1 - azoff1;	/* test for last row */
    else nl2 = naz - ovrlp2;

    printf("unwrapped phase buffer row starting and ending line to copy: %d %d\n", nl1, nl2);
    fwrite((unsigned char *)&u3[nl1][0], sizeof(float), (nl2 - nl1)*width, unwf);

#ifdef fread						/* unswap the bytes after write if necessary */
    for (i = nl1 * width; i < nl2*width; i++) {
      b1 = (unsigned char *)(unw3 + i);
      btmp = b1[0];
      b1[0] = b1[3];
      b1[3] = btmp;
      btmp = b1[1];
      b1[1] = b1[2];
      b1[2] = btmp;
    }
#endif

    if (ferror(unwf)) {
      fprintf(stderr, "\nERROR: write error to output unwrapped phase file\n\n");
      exit(-1);
    }
    fflush(unwf);
  }

  if (azoff + naz1 < nls) {
    nl1 = 0;
    nl2 = nls - (azoff + naz1);
    printf("\nwriting %d blank lines at the end of the file\n", nl2 - nl1);
    for (i = 0; i < width*naz; i++)unw3[i] = 0.0;		/* zero out output buffer array */
    for (j = 0; j < (nl2 - nl1); j++) {
      fwrite((char *)&u3[0][0], sizeof(float), width, unwf);
    }
    fflush(unwf);
  }

  fclose(unwf);

  unwf = fopen(argv[4], "r+b"); 				/* open for update */
  if (unwf == NULL) {
    fprintf(stderr, "\nERROR: cannot reopen output unwrapped phase file: %s\n\n", argv[4]);
    exit(-1);
  }
  printf("output unwrapped phase file: %s\n\n", argv[4]);
  free (u1[0]);
  free(u2[0]);
  free(msk1);
  free(msk2);
  free(m1);
  free(m2);
  free(a1);
  free(a);

  ax3 = (fcomplex *)malloc(sizeof(fcomplex) * width * naz);
  a3 = (fcomplex **)malloc(sizeof(fcomplex *) * naz);

  if (ax3 == NULL || a3 == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for interferogram data buffer\n\n");
    exit(-1);
  }
  printf("memory allocation for interferogram data buffer: %10.3f Mb\n", (sizeof(fcomplex)*naz*width) / 1048576.);
  for (i = 0; i < naz; i++)a3[i] = ax3 + i * width;

  if (wflg == 1) {
    wgt3 = (float *)malloc(sizeof(float) * width * naz);
    w3 = (float **)malloc(sizeof(float *) * naz);
    if (wgt3 == NULL || w3 == NULL) {
      fprintf(stderr, "\nERROR: memory allocation error for weight data buffer\n\n");
      exit(-1);
    }
    for (i = 0; i < naz; i++)w3[i] = wgt3 + i * width;
  }

  if (mflg == 1) {
    msk3 = (unsigned char *)calloc(sizeof(unsigned char), width * naz);
    m3 = (unsigned char **)malloc(sizeof(unsigned char *) * naz);

    if (msk3 == NULL || m3 == NULL) {
      fprintf(stderr, "\nERROR: memory allocation error for mask file buffer\n\n");
      exit(-1);
    }
    printf("memory allocation for mask buffer: %10.3f Mb\n", (naz*width) / 1048576.);
    for (i = 0; i < naz; i++)m3[i] = msk3 + i * width;
  }

  printf("\n***** correcting phase values from initial conditions *****\n");
  ns1 = ceil((double)nls / naz);
  printf("number of input/output segments: %d\n", ns1);

  rdp_float(u3, width, nls, 0, azinit, width, 1, unwf); 	/* read line to get the initial value at ref. point */
  ph0 = u3[0][rinit];

  if (mflg == 1) {
    if (im_flag1 == SUN)rdp_msk_ras(m3, mw, mhgt, 0, azinit, width, 1, r1, ct1, nc1, mskf);
    else rdp_msk_bmp(m3, mw, mhgt, 0, azinit, width, 1, bfh1, bih1, ct1, nc1, mskf);
    if (m3[0][rinit] == 0) {
      printf("\nWARNING: phase reference point masked out and set to 0.0\n\n");
      ph0 = 0.0;
    }
  }

  if (wflg == 1) {
    rdp_float(w3, width, nls, 0, azinit, width, 1, wgtf);
    if (w3[0][rinit] == 0.0) {
      printf("\nWARNING: weight equals 0.0, phase reference set to 0.0\n\n");
      ph0 = 0.0;
    }
  }

  rdp_cpx(a3, width, nls, 0, azinit, width, 1, intf); 		/* see if valid data present at ref. point in the interferogram */
  if (a3[0][rinit].re == 0.0 && a3[0][rinit].im == 0.0) {
    printf("\nWARNING: phase reference point interferogram value (0.0,0.0) (no data)\n\n");
    ph0 = 0.0;
  }

  printf("phase at reference point:  range: %d  azimuth: %d   %10.3f radians\n", rinit + 1, azinit + 1, ph0);
  switch (pflg2) {
  case 0:
    ph1 = fmod(ph0, TWO_PI);
    if (ph1 >  PI)ph1 -= TWO_PI;
    if (ph1 < -PI)ph1 += TWO_PI;
    ph0 = ph0 - ph1;			/* multiples of TWO_PI between wrapped and unwrapped phase */
    printf("phase initialization flag: 0, global phase offset: %10.3f radians\n", ph0);
    break;
  case 1:
    printf("phase initialization flag: 1, global phase offset: %10.3f radians\n", ph0);
    break;
  default:
    fprintf(stderr, "\nERROR: invalid value for init phase flag: %d\n\n", pflg2);
    exit(-1);
  }

  for (k = 0; k < ns1; k++) {
    azoff2 = k * naz;
    naz_rd = rdp_cpx(a3, width, nls, 0, azoff2, width, naz, intf);
    naz_rd = rdp_float(u3, width, nls, 0, azoff2, width, naz, unwf);
    printf("segment: %4d  number of lines: %d  phase offset (rad): %10.3f\n", k + 1, naz_rd, ph0);

    for (i = 0; i < naz_rd; i++) {
      for (j = 0; j < width; j++) {
        if ((a3[i][j].re == 0.0) && (a3[i][j].im == 0.0))u3[i][j] = 0.0;	/* check for 0.0 in the input data */
        if (u3[i][j] != 0.0)u3[i][j] -= ph0;
      }
    }

    if (mflg == 1) {
      if (im_flag1 == SUN)rdp_msk_ras(m3, mw, mhgt, 0, azoff2, width, naz, r1, ct1, nc1, mskf);
      else rdp_msk_bmp(m3, mw, mhgt, 0, azoff2, width, naz, bfh1, bih1, ct1, nc1, mskf);

      for (i = 0; i < naz_rd; i++) {
        for (j = 0; j < width; j++) {
          if (m3[i][j] == 0)u3[i][j] = 0.0;	/* check mask for 0 mask value */
        }
      }
    }

    if (wflg == 1) {				/* check if weight == 0.0 */
      naz_rd = rdp_float(w3, width, nls, 0, azoff2, width, naz, wgtf);
      for (i = 0; i < naz_rd; i++) {
        for (j = 0; j < width; j++) {
          if (w3[i][j] == 0.0)u3[i][j] = 0.0;
        }
      }
    }

    fseek(unwf, (off_t)sizeof(float)*k*naz*width, SEEK_SET);	/* seek to starting line */
    fwrite((unsigned char *)&u3[0][0], sizeof(float), naz_rd*width, unwf);
    if (ferror(unwf)) {
      fprintf(stderr, "\nERROR: Error writing to output file: %s\n\n", argv[4]);
      exit(-1);
    }
    fflush(unwf);
  }
  stop_timing();
  return(0);
}

int rdp_cpx_phase(float **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *cpxf) {
  /*
     read a complex data covering a patch from a data file
     patches are read from a region starting at (line:azoff, sample: roff) in the large file.
     nominal patch size is naz lines x nr samples
     return value is lines actually read from the file and is set to 0 if no lines read

     all entries in the float data array are set to 0.0 if no data read for that location

     28-Nov-2001 clw

     a		output 2D floating point phase derived from complex data
     width	width of input data file, number of samples/line
     roff         range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz          number of azimuth lines in the patch
     cpxf		complex floating point format FILE data structure

  */
  int i, j, jj, naz1, nr1;
  fcomplex *ctmp;

  ctmp = (fcomplex *)malloc(width * 2 * sizeof(float));
  if (ctmp == NULL) {
    fprintf(stderr, "\nERROR subroutine rdp_cpx_phase: memory allocation error for temp complex array\n\n");
    exit(-1);
  }

  if (roff >= width) {
    fprintf(stderr, "\nERROR subroutine rdf_cpx_phase: starting range sample offset exceeds width: %d  %d\n\n", roff, width);
    exit(-1);
  }

  if (azoff >= nls) {
    fprintf(stderr, "\nERROR subroutine rdf_cpx_phase: starting azimuth line exceeds number of lines in the file: %d  %d\n\n", azoff, nls);
    exit(-1);
  }

  if (roff < 0 || azoff < 0) {
    fprintf(stderr, "\nERROR subroutine rdf_cpx_phase:: invalid starting offset roff or azoff < 0: roff: %d   azoff: %d\n", roff, azoff);
    exit(-1);
  }

  if ((naz + azoff) > nls)naz1 = nls - azoff;	/* check if blank lines need to be written at the end of the array */
  else naz1 = naz;

  if ((nr + roff) > width)nr1 = width - roff;	/* check if blank samples need to be written at the end of the lines */
  else nr1 = nr;

  fseek(cpxf, (off_t)sizeof(fcomplex)*azoff*width, SEEK_SET);	/* seek to starting line */

  for (i = 0; i < naz1; i++) {
    fread((char *)ctmp, sizeof(float), 2*width, cpxf);
    if (feof(cpxf)) {
      fprintf(stderr, "\nERROR subroutine rd_phase_cpx: unexpected end of fcomplex file at line: %d\n", azoff + i);
      exit(-1);
    }

    for (j = 0; j < nr1; j++) {				/* convert to phase */
      jj = j + roff;					/* offset in range */
      if ((ctmp[jj].re == 0.0) && (ctmp[jj].im == 0.0))a[i][j] = 0.0;
      else a[i][j] = atan2((double)ctmp[jj].im, (double)ctmp[jj].re);
    }
    if (nr1 < nr) {					/* fill in missing range samples */
      for (j = nr1; j < nr; j++)a[i][j] = 0.0;
    }
  }

  if (naz1 < naz) {					/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++)a[i][j] = 0.0;
    }
  }

  free(ctmp);
  return (naz1);
}

int rdp_float(float **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *fltf) {
  /*
     read floating point data covering a patch from a data file
     patches are read from a region starting at (line:azoff1, sample: roff1) in the large file.
     nominal patch size is naz lines x nr samples
     return value is lines actually read from the file and is set to 0 if no lines read

     all entries in the float data array are set to 0.0 if no data read for that location

     28-Nov-2001 clw

     a		output 2D floating point phase derived from complex data
     width	width of input data file, number of samples/line
     roff         range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz          number of azimuth lines in the patch
     cpxf		complex floating point format FILE data structure

  */
  int i, j, jj, naz1, nr1;
  float *ftmp;

  ftmp = (float *)malloc(width * sizeof(float));
  if (ftmp == NULL) {
    fprintf(stderr, "\nERROR subroutine rdp_float: memory allocation error for temp float array\n\n");
    exit(-1);
  }

  if (roff >= width) {
    fprintf(stderr, "\nERROR subroutine rdp_float: starting range sample offset exceeds width: %d  %d\n\n", roff, width);
    exit(-1);
  }

  if (azoff >= nls) {
    fprintf(stderr, "\nERROR subroutine rdp_float: starting azimuth line exceeds number of lines in the file: %d  %d\n\n", azoff, nls);
    exit(-1);
  }

  if (roff < 0 || azoff < 0) {
    fprintf(stderr, "\nERROR subroutine rdp_float: invalid starting offset roff or azoff < 0: roff: %d   azoff: %d\n", roff, azoff);
    exit(-1);
  }

  if ((naz + azoff) > nls)naz1 = nls - azoff;	/* check if blank lines need to be written at the end of the array */
  else naz1 = naz;

  if ((nr + roff) > width)nr1 = width - roff;	/* check if blank samples need to be written at the end of the lines */
  else nr1 = nr;

  fseek(fltf, (off_t)sizeof(float)*azoff*width, SEEK_SET);	/* seek to start of line */

  for (i = 0; i < naz1; i++) {
    fread((char *)ftmp, sizeof(float), width, fltf);
    if (feof(fltf)) {
      fprintf(stderr, "\nERROR subroutine rd_float: unexpected end of floating point file at line: %d\n", azoff + i);
      exit(-1);
    }

    for (j = 0; j < nr1; j++) {				/* convert to phase */
      jj = j + roff;					/* offset in range */
      a[i][j] = ftmp[jj];
    }
    if (nr1 < nr) {					/* fill in missing range samples */
      for (j = nr1; j < nr; j++)a[i][j] = 0.0;
    }
  }

  if (naz1 < naz) {					/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++)a[i][j] = 0.0;
    }
  }

  free(ftmp);
  return (naz1);
}

int rdp_cpx(fcomplex **a, int width, int nls, int roff, int azoff, int nr, int naz, FILE *fltf) {
  /*
     read fcomplex data covering a patch from a data file
     patches are read from a region starting at (line:azoff1, sample: roff1) in the large file.
     nominal patch size is naz lines x nr samples
     return value is lines actually read from the file and is set to 0 if no lines read

     all entries in the fcomplex data array are set to (0.0,0.0) if no data read for that location

     13-Dec-2001 clw

     a		output 2D fcomplex point phase derived from complex data
     width	width of input data file, number of samples/line
     roff         range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz          number of azimuth lines in the patch
     cpxf		complex floating point format FILE data structure

  */
  int i, j, jj, naz1, nr1;
  fcomplex *xtmp;

  xtmp = (fcomplex *)malloc(width * sizeof(fcomplex));
  if (xtmp == NULL) {
    fprintf(stderr, "\nERROR subroutine rdp_cpx: memory allocation error for temp fcomplex array\n\n");
    exit(-1);
  }

  if (roff >= width) {
    fprintf(stderr, "\nERROR subroutine rdp_cpx: starting range sample offset exceeds width: %d  %d\n\n", roff, width);
    exit(-1);
  }

  if (azoff >= nls) {
    fprintf(stderr, "\nERROR subroutine rdp_cpx: starting azimuth line exceeds number of lines in the file: %d  %d\n\n", azoff, nls);
    exit(-1);
  }

  if (roff < 0 || azoff < 0) {
    fprintf(stderr, "\nERROR subroutine rdp_cpx: invalid starting offset roff or azoff < 0: roff: %d   azoff: %d\n", roff, azoff);
    exit(-1);
  }

  if ((naz + azoff) > nls)naz1 = nls - azoff;	/* check if blank lines need to be written at the end of the array */
  else naz1 = naz;

  if ((nr + roff) > width)nr1 = width - roff;	/* check if blank samples need to be written at the end of the lines */
  else nr1 = nr;

  fseek(fltf, (off_t)sizeof(fcomplex)*azoff*width, SEEK_SET);	/* seek to start of line */

  for (i = 0; i < naz1; i++) {
    fread((char *)xtmp, sizeof(float), 2*width, fltf);
    if (feof(fltf)) {
      fprintf(stderr, "\nERROR subroutine rdp_cpx: unexpected end of fcomplex file at line: %d\n", azoff + i);
      exit(-1);
    }

    for (j = 0; j < nr1; j++) {
      jj = j + roff;					/* offset in range */
      a[i][j] = xtmp[jj];
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
  return (naz1);
}

int rdp_msk_ras(unsigned char **msk, int width, int nls, int roff, int azoff, int nr, int naz,
                RASTER_HDR *r1, COLOR *ct1, int nc1, FILE *mskf) {
  /*
     read mask data covering a patch from a Sun raster file
     patches are read from a region starting at (line:azoff1, sample: roff1) in the large file.
     nominal patch size is naz lines x nr samples
     return value is lines actually read from the file and is set to 0 if no lines read

     all entries in the float data array are set to 0 if no data read for that location
     4-Dec-2001 clw

     msk		output mask derived from raster file
     width	width of input data file, number of samples/line
     roff         range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz          number of azimuth lines in the patch
     r1           Sun rasterfile data structure
     ct1		color table
     nc1          number of entries in the color table

  */
  unsigned char *uctmp;
  int dn, naz1, nr1;
  int i, j, jj;

  uctmp = (unsigned char *)malloc(width);
  if (uctmp == NULL) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_ras: memory allocation error for temp char array\n\n");
    exit(-1);
  }

  if (roff >= width) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_ras: starting range sample offset exceeds width: %d  %d\n\n", roff, width);
    exit(-1);
  }

  if (azoff >= nls) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_ras: starting azimuth line exceeds number of lines in the file: %d  %d\n\n", azoff, nls);
    exit(-1);
  }

  if (roff < 0 || azoff < 0) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_ras: invalid starting offset roff or azoff < 0: roff: %d   azoff: %d\n", roff, azoff);
    exit(-1);
  }

  if ((naz + azoff) > nls)naz1 = nls - azoff;	/* check if blank lines need to be written at the end of the array */
  else naz1 = naz;

  if ((nr + roff) > width)nr1 = width - roff;	/* check if blank samples need to be written at the end of the lines */
  else nr1 = nr;

  for (i = 0; i < naz1; i++) {
    rdras(mskf, r1, i + azoff, (unsigned char *)uctmp);

    for (j = 0; j < nr1; j++) {				/* convert to phase */
      jj = j + roff;
      dn = uctmp[jj];					/* include offset in range */
      if (dn >= nc1) {
        fprintf(stderr, "\nERROR subroutine rdp_msk_ras: invalid entry in mask image file: %d  number of color table entries: %d\n\n", dn, nc1);
        exit(-1);
      }

      if ((ct1[dn].red != 0) || (ct1[dn].green != 0) || (ct1[dn].blue != 0)) msk[i][j] = 255;
      else msk[i][j] = 0;
    }
    if (nr1 < nr) {					/* fill in missing range samples */
      for (j = nr1; j < nr; j++)msk[i][j] = 0;
    }
  }

  if (naz1 < naz) {					/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++)msk[i][j] = 0;
    }
  }

  free(uctmp);
  return (naz1);
}

int rdp_msk_bmp(unsigned char **msk, int width, int nls, int roff, int azoff, int nr, int naz,
                BITMAPFILEHEADER *bfh1, BITMAPINFOHEADER *bih1, COLOR *ct1, int nc1, FILE *mskf) {
  /*
     read mask data covering a patch from a BMP format file
     patches are read from a region starting at (line:azoff1, sample: roff1) in the large file.
     nominal patch size is naz lines x nr samples
     return value is lines actually read from the file and is set to 0 if no lines read

     all entries in the float data array are set to 0 if no data read for that location

     10-Dec-2001 clw

     msk		output mask derived from raster file
     width	width of input data file, number of samples/line
     roff         range offset to starting range sample
     azoff	offset to starting azimuth line of data subset that will be processed
     nr		number of samples/line of the patch
     naz          number of azimuth lines in the patch
     bfh1		BMP BITMAPFILEHEADER
     bih1		BMP BITMAPINFOHEADER
     ct1		color table
     nc1          number of entries in the color table

  */
  unsigned char *uctmp;
  int dn, naz1, nr1;
  int i, j, jj;

  uctmp = (unsigned char *)malloc(width);
  if (uctmp == NULL) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_bmp: memory allocation error for temp char array\n\n");
    exit(-1);
  }

  if (roff >= width) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_bmp: starting range sample offset exceeds width: %d  %d\n\n", roff, width);
    exit(-1);
  }

  if (azoff >= nls) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_bmp: starting azimuth line exceeds number of lines in the file: %d  %d\n\n", azoff, nls);
    exit(-1);
  }

  if (roff < 0 || azoff < 0) {
    fprintf(stderr, "\nERROR subroutine rdp_msk_bmp: invalid starting offset roff or azoff < 0: roff: %d   azoff: %d\n", roff, azoff);
    exit(-1);
  }

  if ((naz + azoff) > nls)naz1 = nls - azoff;	/* check if blank lines need to be written at the end of the array */
  else naz1 = naz;

  if ((nr + roff) > width)nr1 = width - roff;	/* check if blank samples need to be written at the end of the lines */
  else nr1 = nr;

  for (i = 0; i < naz1; i++) {
    rdbmp(mskf, bfh1, bih1, i + azoff, (unsigned char *)uctmp);

    for (j = 0; j < nr1; j++) {				/* convert to phase */
      jj = j + roff;
      dn = uctmp[jj];					/* include offset in range */
      if (dn >= nc1) {
        fprintf(stderr, "\nERROR subroutine rdp_msk_bmp: invalid entry in mask image file: %d  number of color table entries: %d\n\n", dn, nc1);
        exit(-1);
      }

      if ((ct1[dn].red != 0) || (ct1[dn].green != 0) || (ct1[dn].blue != 0)) msk[i][j] = 255;
      else msk[i][j] = 0;
    }
    if (nr1 < nr) {					/* fill in missing range samples */
      for (j = nr1; j < nr; j++)msk[i][j] = 0;
    }
  }

  if (naz1 < naz) {					/* fill blank lines at the bottom of the buffer */
    for (i = naz1; i < naz; i++) {
      for (j = 0; j < nr; j++)msk[i][j] = 0;
    }
  }

  free(uctmp);
  return (naz1);
}











