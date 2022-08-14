/*
 Copyright 2006 clw Gamma Remote Sensing AG
 20-May-2006 Display and Graphics subroutines v1.2
*/
#include "display.h"
#include "bmp_image.h"
#include "rasterfile.h"
#define PPI 75		/* screen pixels per inch */
#define IPM 39  	/* approx. inches/meter */

#define COLOR_MODEL 0	/* color model definition: 0: "old-model" used in software */
/* color model definition: 1: HLS Color Model (double hexcone) */
/* color model definition: 2: HSV Color Model (single hexcone) */

extern void swap(char*, int, int);
int init_rdras(char *f_ras, FILE **r1f, RASTER_HDR **r1, COLOR **ct, int *ncolors, int *nbits, int *width, int *height);
int init_wrras(char *f_ras, FILE **r1f, RASTER_HDR **r1, COLOR *ct, int ncolors, int nbits, int width, int height);
int rdras(FILE *rf, RASTER_HDR *r1, int iline, unsigned char *data);
int wrras(FILE *rf, RASTER_HDR *r1, int iline, unsigned char *data);

int init_rdbmp(char *f_ras, FILE **rf, BITMAPFILEHEADER **bfh, BITMAPINFOHEADER **bih, COLOR **ct, int *nc, int *nbits, int *width, int *height);
int init_wrbmp(char *f_ras, FILE **rf, BITMAPFILEHEADER **bfh, BITMAPINFOHEADER **bih, COLOR *ct, int nc, int nbits, int width, int height);
int rdbmp(FILE *rf, BITMAPFILEHEADER *bfh, BITMAPINFOHEADER *bih, int iline, unsigned char *data);
int wrbmp(FILE *rf, BITMAPFILEHEADER *bfh, BITMAPINFOHEADER *bih, int iline, unsigned char *data);

int ct_8bit(COLOR **ct1, int ct_num, int *nc);		/* generate 8-bit color tables */
int ct_8bit_par(COLOR **ct1, int ct_num, int *nc, double h0, double hrange, double sat, double imin, double imax);		/* more flexible program to generate 8-bit color tables */
int ct_16bit(COLOR **ct1, int ct_num, int *nc);		/* generate 16-bit color tables */
int rd_image(char *fn, unsigned char *image, int im_format, COLOR *ct1, int *nc, int *nbits, int *width, int *height);
int wr_image(char *fn, unsigned char *image, int im_format, COLOR *ct1, int nc, int nbits, int width, int height);

void concave(int nvert, Point2 *point, Win *win, Point2 *ptx, int *np);

double value(double n1, double n2, double hue);
int hls_to_rgb(double hue, double lightness, double saturation, double *r, double *g, double *b);
int hsv_to_rgb(double hue, double lightness, double saturation, double *r, double *g, double *b);
int rgb_to_hsv(double red, double green, double blue, double *hue, double *value, double *saturation);

int init_rdras(char *fname, FILE **rf1, RASTER_HDR **rs1, COLOR **ct1, int *ncolors, int *nbits, int *width, int *height) {
  /*
      initialize reading a SUN raster file image, read the header, and extract color table and parameters
      18-Nov-1999 clw copyright Gamma Remote Sensing AG
  */
  int i;
  RASTER_HDR *r1;
  COLOR *ct;
  FILE *rf;

  rf = fopen(fname, FOPEN_RDONLY_BINARY);
  if (rf == 0) {
    fprintf(stderr, "\nERROR: cannot open SUN raster file: %s\n\n", fname);
    exit( -1);
  }
  *rf1 = rf;

  r1 = (RASTER_HDR *)malloc(sizeof(RASTER_HDR));		/* allocate memory for header */
  if (r1 == NULL)
    fprintf(stderr, "\nERROR: memory allocation error for SUN raster file header\n\n");
  *rs1 = r1;

  fread((char *)r1, sizeof(char), sizeof(RASTER_HDR), rf);	/* read the header */

#ifdef CPU_LITTLE_END
  swap((char *)r1, sizeof(int), sizeof(RASTER_HDR) / sizeof(int));  /* SUN raster file header always BIG ENDIAN */
#endif

  if (r1->ras_magic != RAS_MAGIC) {
    fprintf(stderr, "\nERROR: image file not SUN raster format: %s\n\n", fname);
    return ( -1);
  }
  *width = r1->ras_width;
  *height = r1->ras_height;
  *nbits = r1->ras_depth;

  if ((*nbits != 8) && (*nbits != 24)) {
    fprintf(stderr, "\nERROR: Support only for 8 and 24 bit SUN raster files: %d bits/sample in file: %s\n", r1->ras_depth, fname);
    exit( -1);
  }

  if (r1 -> ras_depth == 8) {
    if (r1->ras_maptype == RMT_EQUAL_RGB) {
      *ncolors = r1->ras_maplength / 3;
      ct = (COLOR *)malloc(*ncolors * sizeof(COLOR));
      if (ct == NULL) {
        fprintf(stderr, "\nERROR: memory allocation error for raster color table for file: %s\n\n", fname);
        exit( -1);
      }
      *ct1 = ct;		/* copy pointer for color table */
      for (i = 0; i < *ncolors; i++)
        ct[i].red = fgetc(rf);
      for (i = 0; i < *ncolors; i++)
        ct[i].green = fgetc(rf);
      for (i = 0; i < *ncolors; i++)
        ct[i].blue = fgetc(rf);
    } else {
      *ncolors = 256;
      ct = (COLOR *)malloc(*ncolors * sizeof(COLOR));
      if (ct == NULL) {
        fprintf(stderr, "\nERROR: memory allocation error for raster color table file: %s\n\n", fname);
        exit( -1);
      }
      *ct1 = ct;			/* copy pointer for color table */
      for (i = 0; i < *ncolors; i++)
        ct[i].red = i;
      for (i = 0; i < *ncolors; i++)
        ct[i].green = i;
      for (i = 0; i < *ncolors; i++)
        ct[i].blue = i;
    }
  } else
    *ncolors = 0;
  return (0);
}

int rdras(FILE *rf, RASTER_HDR *r1, int iline, unsigned char *data) {
  /*
      read a line of a SUN raster format  image
      18-Nov-1999 clw copyright Gamma Remote Sensing AG
  */
  off_t offset;
  size_t nb, nbl;

  nbl = (size_t)(r1->ras_width * r1->ras_depth / 8);		/* number of bytes/line */
  if (nbl % 2 != 0)
    nbl = 2 * (nbl / 2 + 1); 	        /* pad out to even number of bytes width as required by SUN definition */
  offset = (off_t)(sizeof(RASTER_HDR) + r1->ras_maplength) + (off_t)iline * (off_t)nbl;
  fseek(rf, (off_t)(offset), SEEK_SET);

  nb = fread((unsigned char *)data, sizeof(char), nbl, rf);
  /*   for(i=0; i < nbl; i++)printf("iline: %5d  i:%5d  %5d\n",iline,i,(int)data[i]); */

  if (nb != nbl) {
    fprintf(stderr, "\nERROR: reading SUN raster file at line %d\n\n", iline);
    exit( -1);
  }
  return (0);
}

int init_wrras(char *fname, FILE **rf1, RASTER_HDR **rs1, COLOR *ct, int ncolors, int nbits, int width, int height) {
  /*
      initialize SUN raster image, opening file, generating header, and writing to disk
      18-Nov-1999 clw  copyright Gamma Remote Sensing AG
  */

  int i;
  int nbl;
  RASTER_HDR *r1;
  FILE *rf;

  rf = fopen(fname, FOPEN_WRONLY_BINARY);
  if (rf == 0) {
    fprintf(stderr, "\nERROR: cannot open output SUN raster file: %s\n\n", fname);
    exit( -1);
  }
  *rf1 = rf;

  r1 = (RASTER_HDR *)malloc(sizeof(RASTER_HDR));			/* allocate memory for header */
  if (r1 == NULL) {
    fprintf(stderr, "\nERROR: memory allocation error for SUN raster file header\n\n");
    exit( -1);
  }
  *rs1 = r1;

  r1->ras_magic = RAS_MAGIC;
  r1->ras_type = RT_STANDARD;
  r1->ras_width = width;			/* width in samples */
  r1->ras_height = height;			/* number of lines */
  r1->ras_depth = nbits;			/* bits per sample */

  nbl = width * nbits / 8;
  if (nbl % 2 != 0)
    nbl = 2 * (nbl / 2 + 1); 	        /* pad out to even number of bytes width as required by SUN definition */
  r1->ras_length = nbl * height;			/* size of the file in bytes */

  if (nbits == 8) {
    r1->ras_maptype = RMT_EQUAL_RGB;
    r1->ras_maplength = ncolors * 3;
  } else {
    r1->ras_maptype = RMT_NONE;
    r1->ras_maplength = 0;
  }

#ifdef CPU_LITTLE_END
  swap((char *)r1, sizeof(int), sizeof(RASTER_HDR) / sizeof(int));  /* SUN raster file header always BIG ENDIAN */
#endif

  fwrite((char *)r1, 1, sizeof(RASTER_HDR), rf);

#ifdef CPU_LITTLE_END
  swap((char *)r1, sizeof(int), sizeof(RASTER_HDR) / sizeof(int));  /* SUN raster file header always BIG ENDIAN */
#endif

  if (nbits == 8) {
    for (i = 0; i < ncolors; i++)
      fputc(ct[i].red, rf);
    for (i = 0; i < ncolors; i++)
      fputc(ct[i].green, rf);
    for (i = 0; i < ncolors; i++)
      fputc(ct[i].blue, rf);
  }

  fflush(rf);
  return (0);
}

int wrras(FILE *rf, RASTER_HDR *r1, int iline, unsigned char *data) {
  /*
      write a line of data to a SUN raster file
      18-Nov-1999 clw copyright Gamma Remote Sensing AG
  */

  off_t offset;
  size_t nb, nbl, dtv;
  unsigned char tail[1] = {0};

  nbl = r1->ras_width * r1->ras_depth / 8;
  dtv = nbl;
  if (nbl % 2 != 0)
    nbl = 2 * (nbl / 2 + 1); 		/* pad out to even number of bytes width as required by SUN definition */
  offset = (off_t)(sizeof(RASTER_HDR) + r1->ras_maplength) + (off_t)iline * (off_t)nbl;
  fseek(rf, (off_t)(offset), SEEK_SET);

  nb = fwrite((unsigned char *)data, sizeof(char), dtv, rf);
  if (nb != dtv) {
    fprintf(stderr, "\nERROR subroutine wrras: writing SUN raster file at line %d\n\n", iline);
    exit( -1);
  }
  if (nbl > dtv)
    fwrite((unsigned char *)tail, sizeof(char), (nbl - dtv), rf);	/* pad out to the width of the output image file/line */
  return (0);
}

int init_rdbmp(char *fname, FILE **rf1, BITMAPFILEHEADER **bfh1, BITMAPINFOHEADER **bih1, COLOR **ct1, int *ncolors, int *nbits, int *width, int *height) {
  /*
      initilize reading a BMP format image, extract color tables, and decode parameters
      18-Nov-1999 clw copyright Gamma Remote Sensing AG
  */
  int j;
  BMP_RGBQUAD *colors;
  BITMAPFILEHEADER *bfh;
  BITMAPINFOHEADER *bih;
  COLOR *ct;
  FILE *rf;

  rf = fopen(fname, FOPEN_RDONLY_BINARY);
  if (rf == 0) {
    fprintf(stderr, "\nERROR: cannot open BMP raster file: %s\n\n", fname);
    exit( -1);
  }
  *rf1 = rf;

  bfh = (BITMAPFILEHEADER*)malloc(sizeof(BITMAPFILEHEADER));
  if (bfh == NULL)
    fprintf(stderr, "\nERROR: memory allocation error for BMP file information header\n\n");
  *bfh1 = bfh;

  bih = (BITMAPINFOHEADER*)malloc(sizeof(BITMAPINFOHEADER));
  if (bfh == NULL)
    fprintf(stderr, "\nERROR: memory allocation error for BMP file information header\n\n");
  *bih1 = bih;

  fread((char *)&bfh->bfType, 1, 2, rf);
  if (bfh->bfType[0] != 'B' || bfh->bfType[1] != 'M') {
    fprintf(stderr, "\nERROR: input image not BMP format: %s\n\n", fname);
    exit( -1);
  }
  fread((char *)&bfh->bfSize, 1, 4, rf);
  fread((char *)&bfh->bfReserved1, 1, 2, rf);
  fread((char *)&bfh->bfReserved2, 1, 2, rf);
  fread((char *)&bfh->bfOffBits, 1, 4, rf);

#ifndef CPU_LITTLE_END
  swap((char *)&bfh->bfSize, 4, 1);
  swap((char *)&bfh->bfReserved1, 2, 1);
  swap((char *)&bfh->bfReserved2, 2, 1);
  swap((char *)&bfh->bfOffBits, 4, 1);
#endif

  fread((char *)bih, 1, sizeof(BITMAPINFOHEADER), rf);

#ifndef CPU_LITTLE_END 			/* always use little endian byte order */
  swap((char *)&bih->biSize, 4, 1);
  swap((char *)&bih->biWidth, 4, 1);
  swap((char *)&bih->biHeight, 4, 1);
  swap((char *)&bih->biPlanes, 2, 1);
  swap((char *)&bih->biBitCount, 2, 1);
  swap((char *)&bih->biCompression, 4, 1);
  swap((char *)&bih->biXPelsPerMeter, 4, 1);
  swap((char *)&bih->biYPelsPerMeter, 4, 1);
  swap((char *)&bih->biClrUsed, 4, 1);
  swap((char *)&bih->biClrImportant, 4, 1);
#endif

  *width = bih->biWidth;
  *height = bih->biHeight;

  if (bih->biCompression != BI_RGB) {
    fprintf(stderr, "\nERROR: compressed BMP images not supported: %s\n\n", fname);
    exit( -1);
  }
  *nbits = bih->biBitCount;

  if (bih->biBitCount != 24) {
    if (bih->biClrUsed == 0)
      bih->biClrUsed = 1 << bih->biBitCount;
    *ncolors = bih->biClrUsed;

    colors = (BMP_RGBQUAD *)malloc(sizeof(BMP_RGBQUAD) * bih->biClrUsed);
    ct = (COLOR *)malloc(sizeof(COLOR) * bih->biClrUsed);
    *ct1 = ct;
    if (colors == NULL || ct == NULL) {
      fprintf(stderr, "\nERROR: memory allocation failure for BMP color table size: %d  file: %s\n\n", bih->biClrUsed, fname);
      exit( -1);
    }

    fread((char *)colors, 1, bih->biClrUsed*sizeof(BMP_RGBQUAD), rf);
    for (j = 0; j < *ncolors; j++) {
      ct[j].blue = colors[j].blue;
      ct[j].green = colors[j].green;
      ct[j].red = colors[j].red;
    }
    free(colors);
  } else
    *ncolors = 0;

#ifdef DEBUG
  printf("\ninit BMP read: width: %d   height: %d   number of bytes/line: %d  number of colors: %d\n", (int)iline, *width, *height, (int)nbl, *ncolors);
#endif

  return (0);
}

int rdbmp(FILE *rf, BITMAPFILEHEADER *bfh, BITMAPINFOHEADER *bih, int iline, unsigned char *data) {
  /*
       read a line from a BMP format image into array data
       17-Jan-2000 clw copyright Gamma Remote Sensing AG
  */
  off_t offset;
  size_t nb, nbl;
  double bpl;

  bpl = bih->biWidth * bih->biBitCount / 8.0;
  nbl = (int)ceil(bpl);
  if (nbl % 4 != 0)
    nbl = 4 * (nbl / 4 + 1);	/* number of bytes/line must end on a 4 byte boundary */
  offset = (off_t)bfh->bfOffBits + (off_t)((bih->biHeight - iline - 1) * nbl);
  fseek(rf, (off_t)(offset), SEEK_SET);

  nb = fread((unsigned char *)data, sizeof(char), nbl, rf);
  if (nb != nbl) {
    fprintf(stderr, "\nERROR: reading BMP raster file at line %d\n\n", iline);
    exit( -1);
  }
  return (0);
}

int init_wrbmp(char *fname, FILE **rf1, BITMAPFILEHEADER **bfh1, BITMAPINFOHEADER **bih1, COLOR *ct, int ncolors, int nbits, int width, int height) {
  /*
      initialize writing a BMP format image, generating headers, opening output file, and writting to disk
      19-Aug-2003 clw copyright Gamma Remote Sensing AG
  */
  int j;
  BMP_RGBQUAD *colors;
  int nbl;
  double bpl;
  BITMAPFILEHEADER *bfh;
  BITMAPINFOHEADER *bih;
  FILE *rf;

  rf = fopen(fname, FOPEN_WRONLY_BINARY);
  if (rf == 0) {
    fprintf(stderr, "\nERROR: cannot open BMP raster file: %s\n\n", fname);
    exit( -1);
  }
  *rf1 = rf;

  bfh = (BITMAPFILEHEADER*)malloc(sizeof(BITMAPFILEHEADER));
  if (bfh == NULL)
    fprintf(stderr, "\nERROR: memory allocation error for BMP file information header\n\n");
  *bfh1 = bfh;

  bih = (BITMAPINFOHEADER*)malloc(sizeof(BITMAPINFOHEADER));
  if (bfh == NULL)
    fprintf(stderr, "\nERROR: memory allocation error for BMP file information header\n\n");
  *bih1 = bih;

  bpl = (width * nbits) / 8.0;		/* number of bytes/line must end on a 4 byte boundary */
  nbl = (int)ceil(bpl);
  if (nbl % 4 != 0)
    nbl = 4 * (nbl / 4 + 1);

  bfh->bfType[0] = 'B';
  bfh->bfType[1] = 'M';
  bfh->bfSize = 14 + 40 + ncolors * sizeof(BMP_RGBQUAD) + nbl * height;
#ifdef DEBUG

  printf("initialize write BMP file: ncolors: %d  width: %d   height: %d   bytes/line: %d\n", ncolors, width, height, nbl);
#endif

  bfh->bfReserved1 = 0;
  bfh->bfReserved2 = 0;
  bfh->bfOffBits = (14 + 40 + (ncolors * sizeof(BMP_RGBQUAD)));

#ifndef CPU_LITTLE_END
  swap((char *)&(bfh->bfSize), 4, 1);
  swap((char *)&(bfh->bfReserved1), 2, 1);
  swap((char *)&(bfh->bfReserved2), 2, 1);
  swap((char *)&(bfh->bfOffBits), 4, 1);
#endif

  fwrite((char *)&(bfh->bfType), 1, 2, rf);		/* write out the BMP file header */
  fwrite((char *)&(bfh->bfSize), 1, 4, rf);
  fwrite((char *)&(bfh->bfReserved1), 1, 2, rf);
  fwrite((char *)&(bfh->bfReserved2), 1, 2, rf);
  fwrite((char *)&(bfh->bfOffBits), 1, 4, rf);

#ifndef CPU_LITTLE_END
  swap((char *)&(bfh->bfSize), 4, 1);
  swap((char *)&(bfh->bfReserved1), 2, 1);
  swap((char *)&(bfh->bfReserved2), 2, 1);
  swap((char *)&(bfh->bfOffBits), 4, 1);
#endif

  bih->biSize = 40;
  bih->biWidth = width;
  bih->biHeight = height;
  bih->biPlanes = 1;
  bih->biBitCount = nbits;
  bih->biCompression = BI_RGB;
  bih->biSizeImage = nbl * height;
  bih->biXPelsPerMeter = PPI * IPM;
  bih->biYPelsPerMeter = PPI * IPM;
  bih->biClrUsed = ncolors;
  bih->biClrImportant = 0;

#ifndef CPU_LITTLE_END 			/* always use little endian byte order */
  swap((char *)&bih->biSize, 4, 1);
  swap((char *)&bih->biWidth, 4, 1);
  swap((char *)&bih->biHeight, 4, 1);
  swap((char *)&bih->biPlanes, 2, 1);
  swap((char *)&bih->biBitCount, 2, 1);
  swap((char *)&bih->biCompression, 4, 1);
  swap((char *)&bih->biXPelsPerMeter, 4, 1);
  swap((char *)&bih->biYPelsPerMeter, 4, 1);
  swap((char *)&bih->biClrUsed, 4, 1);
  swap((char *)&bih->biClrImportant, 4, 1);
#endif

  fwrite((char *)bih, 1, sizeof(BITMAPINFOHEADER), rf);

#ifndef CPU_LITTLE_END 			/* always use little endian byte order */
  swap((char *)&bih->biSize, 4, 1);
  swap((char *)&bih->biWidth, 4, 1);
  swap((char *)&bih->biHeight, 4, 1);
  swap((char *)&bih->biPlanes, 2, 1);
  swap((char *)&bih->biBitCount, 2, 1);
  swap((char *)&bih->biCompression, 4, 1);
  swap((char *)&bih->biXPelsPerMeter, 4, 1);
  swap((char *)&bih->biYPelsPerMeter, 4, 1);
  swap((char *)&bih->biClrUsed, 4, 1);
  swap((char *)&bih->biClrImportant, 4, 1);
#endif

  if (bih->biBitCount != 24) {
    colors = (BMP_RGBQUAD *)malloc(sizeof(BMP_RGBQUAD) * ncolors);
    if (colors == NULL) {
      fprintf(stderr, "\nERROR: memory allocation failure for BMP color table size: %d\n\n", bih->biClrUsed);
      exit( -1);
    }

    for (j = 0; j < ncolors; j++) {
      colors[j].red = ct[j].red;
      colors[j].green = ct[j].green;
      colors[j].blue = ct[j].blue;
      colors[j].reserved = 0;
    }
    fwrite((char *)colors, 1, ncolors*sizeof(BMP_RGBQUAD), rf);
  } else
    ncolors = 0;
  fflush(rf);
  return (0);
}

int wrbmp(FILE *rf, BITMAPFILEHEADER *bfh, BITMAPINFOHEADER *bih, int iline, unsigned char *data) {
  /*
       write a line to a BMP format image file. Pad out to the required length to satisfy BMP
       format constraints.
       17-Jan-2000 clw copyright Gamma Remote Sensing AG

  */
  off_t offset;
  size_t nb, nbl, dtv;

  double bpl;
  unsigned char tail[3] = {0, 0, 0};

  bpl = (bih->biWidth * bih->biBitCount) / 8.0;		/* number of bytes/line must end on a 4 byte boundary */
  nbl = (size_t)((int)ceil(bpl));
  dtv = nbl;						/* number of bytes in the data/line */
  if (nbl % 4 != 0)
    nbl = 4 * (nbl / 4 + 1); 		/* number of bytes per line in the output BMP file */

  offset = (off_t)bfh->bfOffBits + (off_t)((bih->biHeight - iline - 1) * nbl);
  fseek(rf, (off_t)(offset), SEEK_SET);

  nb = fwrite((unsigned char *)data, sizeof(char), dtv, rf);
  if (nb != dtv) {
    fprintf(stderr, "\nERROR subroutine wrbmp: writing BMP raster file at line %d\n\n", iline);
    exit( -1);
  }
  if (nbl > dtv)
    fwrite((unsigned char *)tail, sizeof(char), (nbl - dtv), rf);	/* pad out to the number of bytes/line in the BMP image */
  return (0);
}

#define NC8 256

int ct_8bit(COLOR **ct1, int ct_num, int *nc) {
  /*
       generate 8-bit color tables for Gamma Software
       8-Dec-1999  clw  copyright Gamma Remote Sensing AG

       ct_8bit      generate 8 bit color tables
       ct1		pointer to color table
       ct_num	color table number:
    		0  greyscale 0 --> 255
    		1  RMG mag/phase color table 4 bits intensity, 4 bits phase
                    2  RMG phase 8 bits
    		3  phase unwrapping color table (RMG colors + flags)
    		4  4-bit phase only
    		5  4-bit greyscale only
       nc		number of colors in the table

  */
  int i, j, k;
  double intensity, saturation, hue, imin = 0., imax = 1.;
  double rval, gval, bval;
  COLOR *ct;
  COLOR	red = {255, 50, 50}, blue = {50, 50, 255}, white = {255, 255, 255};
  COLOR salmon = {235, 140, 100}, green = {0, 255, 0}, yellow = {255, 255, 0};
  COLOR black = {0, 0, 0}, cyan = {0, 255, 255}, magenta = {255, 0, 255};


  ct = (COLOR *)malloc(NC8 * sizeof(COLOR));
  if (ct == 0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit: error in allocation of color table\n\n");
    exit( -1);
  }
  *ct1 = ct;

  for (i = 0; i < NC8; i++) {
    ct[i].red = (unsigned char)0;
    ct[i].green = (unsigned char)0;
    ct[i].blue = (unsigned char)0;
  }

  switch (ct_num) {
  case GREY: 			/* greyscale */
    for (i = 0; i < NC8; i++) {
      ct[i].red = (unsigned char)i;
      ct[i].green = (unsigned char)i;
      ct[i].blue = (unsigned char)i;
    }
    *nc = NC8;
    break;

  case MPH: 			/* 4 bits intensity 4 bits phase using "old RMG", HLS, or HSV color model */

    if (COLOR_MODEL == 0) {			/* old color model */
      for (j = 0; j < 16; j++) {
        k = 16 * j;
        ct[j + 240].red = 3 * MIN(MIN(k, 85), 255 - k);
        ct[j + 240].green = 3 * MIN(MAX(k - 85, 85 - k), 85);
        ct[j + 240].blue = 3 * MIN(MAX(k - 170, 170 - k), 85);
        ct[j + 240].red = (int)((int)ct[j + 240].red / 255. * 200. + 55.);
        ct[j + 240].green = (int)((int)ct[j + 240].green / 255. * 200. + 55.);
        ct[j + 240].blue = (int)((int)ct[j + 240].blue / 255. * 200. + 55.);
      }

      for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
          ct[i*16 + j].red = (int)((int)ct[j + 240].red * i / 15.);
          ct[i*16 + j].green = (int)((int)ct[j + 240].green * i / 15.);
          ct[i*16 + j].blue = (int)((int)ct[j + 240].blue * i / 15.);
        }
      }

    } else {				/* HSV or HLS color model */
      if (COLOR_MODEL == 1) {		/* HLS color model (double hexcone) */
        saturation = 0.75;		/* define fixed saturation value level */
        /* saturation = 1.0 is not ideal for the  */
        /* representation of the intensity  */

        /* with HLS model it is not recommended to use the   */
        /* full intensity range (0.0 to 1.0) because color  */
        /* information gets lost near the two extremes 0.0 and 1.0 */
        imin = 0.15;			/* minimum intensity value */
        imax = 0.85;			/* maximum intensity value */
      } else {				/* COLOR_MODEL == 2, HSV color model (single hexcone) */
        saturation = 0.75;		/* define fixed saturation value level */
        /* saturation = 1.0 is not ideal for the  */
        /* representation of the intensity  */

        /* with HSV model high intensity values are no problem */
        imin = 0.15;			/* minimum intensity value */
        imax = 1.00;			/* maximum intensity value (typically 1.0) */
      }

      for (i = 0; i < 16; i++) {		/* intensity levels */
        if (i == 0) {
          intensity = 0.;			/* include black */
        } else {
          intensity = (double)(i * (imax - imin) / 14.);	/* calculate intensity value level */
        }
        for (j = 0; j < 16; j++) {		/* hue levels */
          hue = (double)j * 360. / 16. + 180.;

          if (COLOR_MODEL == 1) {	/* HLS color model (double hexcone) */
            k = hls_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
          } else {				/* COLOR_MODEL == 2, HSV color model (single hexcone) */
            k = hsv_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
          }
          ct[i*16 + j].red = rval;
          ct[i*16 + j].green = gval;
          ct[i*16 + j].blue = bval;
        }
      }
    }

    *nc = NC8;
    break;

  case PHASE: 		/* 8 bits phase using "old RMG", HLS, or HSV color model */
    if (COLOR_MODEL == 0) {			/* old color model */
      for (j = 0; j < 256; j++) {
        ct[j].red = 3 * MIN(MIN(j, 85), 255 - j);
        ct[j].green = 3 * MIN(MAX(j - 85, 85 - j), 85);
        ct[j].blue = 3 * MIN(MAX(j - 170, 170 - j), 85);
        ct[j].red = (int)((int)ct[j].red / 255. * 200. + 55.);
        ct[j].green = (int)((int)ct[j].green / 255. * 200. + 55.);
        ct[j].blue = (int)((int)ct[j].blue / 255. * 200. + 55.);
      }
    } else {				/* HSV or HLS color model */
      if (COLOR_MODEL == 1) {		/* HLS color model (double hexcone) */
        saturation = 0.75;		/* define fixed saturation value level */
        /* saturation = 1.0 is not ideal for the  */
        /* representation of the intensity  */

        /* with HLS model it is not recommended to use the   */
        /* full intensity range (0.0 to 1.0) because color  */
        /* information gets lost near the two extremes 0.0 and 1.0 */
        intensity = 0.5;			/* fixed intensity value */
      } else {				/* COLOR_MODEL == 2, HSV color model (single hexcone) */
        saturation = 0.75;		/* define fixed saturation value level */
        /* saturation = 1.0 is not ideal for the  */
        /* representation of the intensity  */

        /* with HSV model high intensity values are no problem */
        intensity = 0.5;			/* fixed intensity value */
      }

      for (j = 0; j < 255; j++) {		/* hue levels */
        hue = (double)j * 360. / 256. + 180.;
        if (COLOR_MODEL == 1) {	/* HLS color model (double hexcone) */
          k = hls_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
        } else {				/* COLOR_MODEL == 2, HSV color model (single hexcone) */
          k = hsv_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
        }
        ct[j].red = rval;
        ct[j].green = gval;
        ct[j].blue = bval;
      }
    }
    *nc = NC8;
    break;

  case PHUNW:
    for (j = 0;j < 96; j++) {
      ct[j].red = Min(Min(7 * j, 224), 7 * (96 - j)) + 31;
      ct[j].green = Min(7 * Max(j - 32, 32 - j), 224) + 31;
      ct[j].blue = Min(7 * Max(64 - j, j - 64), 224) + 31;

      ct[j + 96].red = ct[j].red * .6;	/* dimmed out for non-unwrapped areas */
      ct[j + 96].green = ct[j].green * .6;
      ct[j + 96].blue = ct[j].blue * .6;
    }
    ct[WHITE] = white;
    ct[BLACK] = black;
    ct[RED] = red;
    ct[BLUE] = blue;
    ct[SALMON] = salmon;
    ct[CYAN] = cyan;
    ct[GREEN] = green;
    ct[YELLOW] = yellow;
    ct[MAGENTA] = magenta;
    *nc = NC8;
    break;

  case PH_4BIT: 			/* 4-bit phase only */
    if (COLOR_MODEL == 0) {			/* old color model */
      for (j = 0; j < 16; j++) {
        k = 16 * j;
        ct[j + 240].red = 3 * MIN(MIN(k, 85), 255 - k);
        ct[j + 240].green = 3 * MIN(MAX(k - 85, 85 - k), 85);
        ct[j + 240].blue = 3 * MIN(MAX(k - 170, 170 - k), 85);
        ct[j + 240].red = (int)((int)ct[j + 240].red / 255. * 200. + 55.);
        ct[j + 240].green = (int)((int)ct[j + 240].green / 255. * 200. + 55.);
        ct[j + 240].blue = (int)((int)ct[j + 240].blue / 255. * 200. + 55.);
      }
      for (i = 0; i < 16; i++) {
        for (j = 0; j < 16; j++) {
          ct[i*16 + j].red = (int)((int)ct[j + 240].red);
          ct[i*16 + j].green = (int)((int)ct[j + 240].green);
          ct[i*16 + j].blue = (int)((int)ct[j + 240].blue);
        }
      }
    } else {				/* HSV or HLS color model */
      if (COLOR_MODEL == 1) {		/* HLS color model (double hexcone) */
        saturation = 0.75;		/* define fixed saturation value level */
        /* saturation = 1.0 is not ideal for the  */
        /* representation of the intensity  */

        /* with HLS model it is not recommended to use the   */
        /* full intensity range (0.0 to 1.0) because color  */
        /* information gets lost near the two extremes 0.0 and 1.0 */
        intensity = 0.5;			/* fixed intensity value */
      } else {				/* COLOR_MODEL == 2, HSV color model (single hexcone) */
        saturation = 0.75;		/* define fixed saturation value level */
        /* saturation = 1.0 is not ideal for the  */
        /* representation of the intensity  */

        /* with HSV model high intensity values are no problem */
        intensity = 0.5;			/* fixed intensity value */
      }

      for (i = 0; i < 16; i++) {		/* intensity levels */
        for (j = 0; j < 16; j++) {		/* hue levels */
          hue = (double)j * 360. / 16. + 180.;

          if (COLOR_MODEL == 1) {	/* HLS color model (double hexcone) */
            k = hls_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
          } else {				/* COLOR_MODEL == 2, HSV color model (single hexcone) */
            k = hsv_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
          }
          ct[i*16 + j].red = rval;
          ct[i*16 + j].green = gval;
          ct[i*16 + j].blue = bval;
        }
      }
    }

    *nc = NC8;
    break;

  case GREY_4BIT: 				/* 4-bit greyscale only */
    for (i = 0; i < 16; i++) {
      for (j = 0; j < 16; j++) {
        ct[i*16 + j].red = (int)(255. * i / 15.);
        ct[i*16 + j].green = (int)(255. * i / 15.);
        ct[i*16 + j].blue = (int)(255. * i / 15.);
      }
    }
    *nc = NC8;
    break;

  case CH_NEUT_LSNR:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LSNR) > 0)
        ct[i] = salmon;
      if ((int)( (unsigned char)i & GUID) > 0)
        ct[i] = white;
      if ((int)( (unsigned char)i & PLUS) > 0)
        ct[i] = red;
      if ((int)( (unsigned char)i & MINU) > 0)
        ct[i] = blue;
    }
    *nc = NC8;
    break;

  case NEUT_LSNR:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LSNR) > 0)
        ct[i] = salmon;
      if ((int)( (unsigned char)i & GUID) > 0)
        ct[i] = white;
    }
    *nc = NC8;
    break;

  case CH_NEUT_LSNR_CUTS:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LSNR) > 0)
        ct[i] = salmon;
      if ((int)( (unsigned char)i & CUT) > 0)
        ct[i] = yellow;
      if ((int)( (unsigned char)i & GUID) > 0)
        ct[i] = white;
      if ((int)( (unsigned char)i & PLUS) > 0)
        ct[i] = red;
      if ((int)( (unsigned char)i & MINU) > 0)
        ct[i] = blue;
    }
    *nc = NC8;
    break;

  case CUTS_LSNR:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LSNR) > 0)
        ct[i] = salmon;
      if ((int)( (unsigned char)i & CUT) > 0)
        ct[i] = yellow;
    }
    *nc = NC8;
    break;

  case CUTS_LAWN:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LAWN) > 0)
        ct[i] = cyan;
      if ((int)( (unsigned char)i & CUT) > 0)
        ct[i] = yellow;
    }
    *nc = NC8;
    break;

  case LAWN_LSNR:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LAWN) > 0)
        ct[i] = cyan;
      if ((int)( (unsigned char)i & LSNR) > 0)
        ct[i] = salmon;
    }
    *nc = NC8;
    break;

  case LAWN_LSNR_CUTS:
    for (i = 0; i < 256; i++) {
      if ((int)( (unsigned char)i & LAWN) > 0)
        ct[i] = cyan;
      if ((int)( (unsigned char)i & LSNR) > 0)
        ct[i] = salmon;
      if ((int)( (unsigned char)i & CUT) > 0)
        ct[i] = yellow;
    }
    *nc = NC8;
    break;

  default:
    printf("\nWARNING: ct_8bit: invalid color table selected: %d\n\n", ct_num);
    return ( -1);
    break;
  }
  return (0);
}

int ct_8bit_par(COLOR **ct1, int ct_num, int *nc, double h0, double hrange, double sat, double imin, double imax) {		/* more flexible program to generate 8-bit color tables */
  /*
       generate 8-bit color tables for Gamma Software
       5-Mar-2001  uw copyright Gamma Remote Sensing AG

       ct_8bit_par  generate 8 bit color tables (with increased flexibility)
       ct1	pointer to color table
       ct_num	color table type number (as defined in display.h):
    		GREY		1
    		MPH_OLD		20
    		MPH_HLS		21
    		MPH_HSV		22
    		MPH_SIN		23
    		MPH_TOPO1	24
    		PHASE_OLD	30
    		PHASE_HLS	31
    		PHASE_HSV	32
    		PHASE_SIN	33
    		PHASE_TOPO1	34
       nc		number of colors in the table
       h0		starting hue (0.0 ... 360.0)
       hrange	ending hue = h0 + hrange, negative hrange values
       		correspond to an color cycle sense inversion
    		hrange (-360.0 ... 360.0)
       sat		saturation value (0.0 ... 1.0)
       imin		minimum intensity value (0.0 ... 1.0)
       imax		maximum intensity value (0.0 ... 1.0)
       		for imax < imin : inversion of brightness
  */

  int i, j, k, ncols;
  double intensity, saturation, hue;
  double rval, gval, bval;
  COLOR *ct;

  /*** check validity of input parameters ***/
  if (h0 < 0.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: h0 < 0.0 not valid: %12.6ef\n\n", h0);
    exit( -1);
  }
  if (h0 > 360.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: h0 > 360.0 not valid  %12.6ef\n\n", h0);
    exit( -1);
  }
  if (hrange < -360.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: hrange < -360.0 not valid: %12.6ef\n\n", hrange);
    exit( -1);
  }
  if (hrange > 360.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: hrange >  360.0 not valid  %12.6ef\n\n", hrange);
    exit( -1);
  }
  if (sat < 0.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: sat < 0.0 not valid: %12.6ef\n\n", sat);
    exit( -1);
  }
  if (sat > 1.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: sat > 1.0 not valid  %12.6ef\n\n", sat);
    exit( -1);
  }
  if (imin < 0.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: imin < 0.0 not valid: %12.6ef\n\n", imin);
    exit( -1);
  }
  if (imin > 1.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: imin > 1.0 not valid  %12.6ef\n\n", imin);
    exit( -1);
  }
  if (imax < 0.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: imax < 0.0 not valid: %12.6ef\n\n", imax);
    exit( -1);
  }
  if (imax > 1.0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: imax > 1.0 not valid  %12.6ef\n\n", imax);
    exit( -1);
  }

  ncols = 256;
  *nc = ncols;
  ct = (COLOR *)malloc(ncols * sizeof(COLOR));
  if (ct == 0) {
    fprintf(stderr, "\nERROR subroutine ct_8bit_par: error in allocation of color table\n\n");
    exit( -1);
  }
  *ct1 = ct;

  for (i = 0; i < ncols; i++) {
    ct[i].red = (unsigned char)0;
    ct[i].green = (unsigned char)0;
    ct[i].blue = (unsigned char)0;
  }

  switch (ct_num) {

  case GREY: 			/* greyscale */
    for (i = 0; i < ncols; i++) {
      ct[i].red = (unsigned char)i;
      ct[i].green = (unsigned char)i;
      ct[i].blue = (unsigned char)i;
    }
    break;

  case MPH_OLD: 			/* 4 bits intensity 4 bits phase using "old RMG", HLS, or HSV color model */
    /*** sat,h0,hrange values are not used ***/
    for (j = 0; j < 16; j++) {
      hue = (double)h0 + j * (hrange / 16.);
      while (hue > 360.0)
        hue -= 360.;	/* limit to (0.0 ... 360.) range */
      while (hue < 0.0)
        hue += 360.;
      hue = hue * 256. / 360.;

      ct[j + 240].red = 3 * MIN(MIN(hue, 85), 255 - hue);
      ct[j + 240].green = 3 * MIN(MAX(hue - 85, 85 - hue), 85);
      ct[j + 240].blue = 3 * MIN(MAX(hue - 170, 170 - hue), 85);
      ct[j + 240].red = (int)((int)ct[j + 240].red / 255. * 200. + 55.);
      ct[j + 240].green = (int)((int)ct[j + 240].green / 255. * 200. + 55.);
      ct[j + 240].blue = (int)((int)ct[j + 240].blue / 255. * 200. + 55.);
    }

    i = 0;				/* case hue == 0.0 */
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */

    for (i = 1; i < 16; i++) {		/*** 15 intensity levels ***/
      intensity = (double)imin + (i - 1) * (imax - imin) / 14.;	/* calculate double intensity value (0.0...1.0) */
      ct[i].red = (int)(intensity * 255.);
      ct[i].green = (int)(intensity * 255.);
      ct[i].blue = (int)(intensity * 255.);

      for (j = 0; j < 16; j++) {	/* hue */
        ct[i*16 + j].red = (int)((int)ct[j + 240].red * intensity);
        ct[i*16 + j].green = (int)((int)ct[j + 240].green * intensity);
        ct[i*16 + j].blue = (int)((int)ct[j + 240].blue * intensity);
      }
    }
    break;

  case MPH_HLS: 			/* 4 bits intensity 4 bits phase using HLS color model */

    /* HLS color model (double hexcone) */
    saturation = sat;
    /* saturation = 1.0 is not ideal for the  */
    /* representation of the intensity  */
    /* with HLS model it is not recommended to use the   */
    /* full intensity range (0.0 to 1.0) because color  */
    /* information gets lost near the two extremes 0.0 and 1.0 */

    i = 0;				/* case hue == 0.0 */
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */

    for (i = 1; i < 16; i++) {		/* 15 intensity levels */
      intensity = (double)imin + (i - 1) * (imax - imin) / 14.;	/* calculate double intensity value (0.0...1.0) */
      ct[i].red = (int)((double)intensity * 255.);
      ct[i].green = (int)((double)intensity * 255.);
      ct[i].blue = (int)((double)intensity * 255.);

      for (j = 0; j < 16; j++) {	/* hue levels */
        hue = (double)h0 + j * (hrange / 16.);
        while (hue > 360.0)
          hue -= 360.;			/* limit to (0.0 ... 360.) range */
        while (hue < 0.0)
          hue += 360.;

        /* HLS color model (double hexcone) */
        k = hls_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
        ct[i*16 + j].red = rval;
        ct[i*16 + j].green = gval;
        ct[i*16 + j].blue = bval;
      }
    }

    break;

  case MPH_HSV: 				/* 4 bits intensity 4 bits phase using HSV color model */
    /* HSV color model (single hexcone) */
    saturation = sat;
    /* saturation = 1.0 is not ideal for the  */
    /* representation of the intensity  */
    /* With HSV model it is not recommended to use the   */
    /* full intensity range (0.0 to 1.0) because color  */
    /* information gets lost near the lower extreme 0.0 */

    i = 0;				/*** case hue == 0.0 ***/
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;	/* 0: black (= no data) */

    for (i = 1; i < 16; i++) {		/*** 15 intensity levels ***/
      intensity = (double)imin + (i - 1) * (imax - imin) / 14.;	/* calculate double intensity value (0.0...1.0) */
      ct[i].red = (int)((double)intensity * 255.);
      ct[i].green = (int)((double)intensity * 255.);
      ct[i].blue = (int)((double)intensity * 255.);

      for (j = 0; j < 16; j++) {		/* hue levels */
        hue = (double)h0 + j * (hrange / 16.);
        while (hue > 360.0)
          hue -= 360.;	/* limit to (0.0 ... 360.) range */
        while (hue < 0.0)
          hue += 360.;

        /* HSV color model (single hexcone) */
        k = hsv_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
        ct[i*16 + j].red = rval;
        ct[i*16 + j].green = gval;
        ct[i*16 + j].blue = bval;
      }
    }

    break;

  case MPH_SIN: 			/* 4 bits intensity 4 bits phase using SIN color model */
    /* saturation value is not used  */
    i = 0;				/*** case hue == 0.0 ***/
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */

    for (i = 1; i < 16; i++) {		/*** 15 intensity levels ***/
      intensity = (double)imin + (i - 1) * (imax - imin) / 14.;	/* calculate double intensity value (0.0...1.0) */
      ct[i].red = (int)((double)intensity * 255.);
      ct[i].green = (int)((double)intensity * 255.);
      ct[i].blue = (int)((double)intensity * 255.);

      for (j = 0; j < 16; j++) {	/* hue levels */
        hue = (double)h0 + j * (hrange / 16.);
        while (hue > 360.0)
          hue -= 360.;			/* limit to (0.0 ... 360.) range */
        while (hue < 0.0)
          hue += 360.;

        ct[i*16 + j].red = (unsigned char)((double)intensity * 255. * (0.5 + 0.5 * sin(hue * DTR - TWO_PI / 3)));
        ct[i*16 + j].green = (unsigned char)((double)intensity * 255. * (0.5 + 0.5 * sin(hue * DTR )));
        ct[i*16 + j].blue = (unsigned char)((double)intensity * 255. * (0.5 + 0.5 * sin(hue * DTR + TWO_PI / 3)));
      }
    }

    break;

  case PHASE_OLD: 			/* 8 bits phase using "old RMG" color model */
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */

    for (j = 1; j < 256; j++) {
      k = j - 1;
      hue = (double)h0 + k * (hrange / 255.);
      while (hue > 360.0)
        hue -= 360.;			/* limit to (0.0 ... 360.) range */
      while (hue < 0.0)
        hue += 360.;
      hue = hue * 256. / 360.;
      ct[j].red = 3 * MIN(MIN(hue, 85), 255 - hue);
      ct[j].green = 3 * MIN(MAX(hue - 85, 85 - hue), 85);
      ct[j].blue = 3 * MIN(MAX(hue - 170, 170 - hue), 85);
      ct[j].red = (int)((int)ct[j].red / 255. * 200. + 55.);
      ct[j].green = (int)((int)ct[j].green / 255. * 200. + 55.);
      ct[j].blue = (int)((int)ct[j].blue / 255. * 200. + 55.);
    }

    break;

  case PHASE_HLS: 			/* 8 bits phase using HLS color model */

    /* HLS color model (double hexcone) */
    saturation = sat;
    /* saturation = 1.0 is not ideal for the  */
    /* representation of the intensity  */
    /* with HLS model it is not recommended to use the   */
    /* full intensity range (0.0 to 1.0) because color  */
    /* information gets lost near the two extremes 0.0 and 1.0 */
    intensity = (imin + imax) / 2.;	/* fixed intensity value */

    i = 0;				/*** case hue == 0.0 ***/
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */


    for (j = 1; j < 256; j++) {		/* 255 hue levels */
      k = j - 1;
      hue = (double)h0 + k * (hrange / 255.);
      while (hue > 360.0)
        hue -= 360.;			/* limit to (0.0 ... 360.) range */
      while (hue < 0.0)
        hue += 360.;

      /* HLS color model (double hexcone) */
      k = hls_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
      ct[j].red = rval;
      ct[j].green = gval;
      ct[j].blue = bval;
    }

    break;

  case PHASE_HSV: 			/* 8 bits phase using HSV color model */
    /* HSV color model (single hexcone) */
    saturation = sat;
    /* saturation = 1.0 is not ideal for the  */
    /* representation of the intensity  */
    /* with HLS model it is not recommended to use the   */
    /* full intensity range (0.0 to 1.0) because color  */
    /* information gets lost near the two extremes 0.0 and 1.0 */
    intensity = (imin + imax) / 2.;	/* fixed intensity value */

    i = 0;				/*** case hue == 0.0 ***/
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */


    for (j = 1; j < 256; j++) {		/* 255 hue levels */
      k = j - 1;
      hue = (double)h0 + k * (hrange / 255.);
      while (hue > 360.0)
        hue -= 360.;			/* limit to (0.0 ... 360.) range */
      while (hue < 0.0)
        hue += 360.;

      /* HLS color model (double hexcone) */
      k = hsv_to_rgb(hue, intensity, saturation, &rval, &gval, &bval);
      ct[j].red = rval;
      ct[j].green = gval;
      ct[j].blue = bval;
    }

    break;

  case PHASE_SIN: 			/* 8 bits phase using SIN color model */
    intensity = (imin + imax) / 2.;	/* fixed intensity value */

    i = 0;				/*** case hue == 0.0 ***/
    ct[0].red = 0;
    ct[0].green = 0;
    ct[0].blue = 0;			/* 0: black (= no data) */

    for (j = 1; j < 256; j++) {		/* 255 hue levels */
      k = j - 1;
      hue = (double)h0 + k * (hrange / 255.);
      while (hue > 360.0)
        hue -= 360.;			/* limit to (0.0 ... 360.) range */
      while (hue < 0.0)
        hue += 360.;

      ct[j].red = (unsigned char)((double)intensity * 255. * (0.5 + 0.5 * sin(hue * DTR - TWO_PI / 3)));
      ct[j].green = (unsigned char)((double)intensity * 255. * (0.5 + 0.5 * sin(hue * DTR )));
      ct[j].blue = (unsigned char)((double)intensity * 255. * (0.5 + 0.5 * sin(hue * DTR + TWO_PI / 3)));
    }

    break;

  default:
    printf("\nWARNING: ct_8bit: invalid color table selected: %d\n\n", ct_num);
    return ( -1);
    break;
  }
  return (0);
}

#define NC16 65536
int ct_16bit(COLOR **ct1, int ct_num , int *nc) {
  /*
      generate 16 bit color tables
      17-Jan-2000 clw  copyright Gamma Remote Sensing AG

      ct_16bit  generate 16 bit color tables
      ct1		pointer to color table
      ct_num	color table number:
    		0  greyscale 0 --> 255
    		1  RMG mag/phase color table 8 bits intensity, 8 bits phase
                    2  RMG phase 16 bit
      nc		number of colors in the table
  */
  int i, j;
  COLOR *ct;

  COLOR black = {0, 0, 0};

  ct = (COLOR *)malloc(NC16 * sizeof(COLOR));
  if (ct == 0) {
    fprintf(stderr, "\nERROR subroutine ct_16bit: error in allocation of color table\n\n");
    exit( -1);
  }
  *ct1 = ct;

  for (i = 0; i < NC16; i++)
    ct[i] = black;

  switch (ct_num) {
  case GREY: 			/* greyscale */
    for (i = 0; i < 65536; i++) {
      ct[i].red = (unsigned char)(i / 256);
      ct[i].green = (unsigned char)(i / 256);
      ct[i].blue = (unsigned char)(i / 256);
    }
    *nc = 65536;
    break;

  case MPH: 			/* RMG 8 bits intensity 8 bits phase */
    for (j = 0; j < 256; j++) {
      ct[j + 65280].red = 3 * MIN(MIN(j, 85), 255 - j);
      ct[j + 65280].green = 3 * MIN(MAX(j - 85, 85 - j), 85);
      ct[j + 65280].blue = 3 * MIN(MAX(j - 170, 170 - j), 85);
      ct[j + 65280].red = (int)((int)ct[j + 65280].red / 255. * 223. + 32.);
      ct[j + 65280].green = (int)((int)ct[j + 65280].green / 255. * 223. + 32.);
      ct[j + 65280].blue = (int)((int)ct[j + 65280].blue / 255. * 223. + 32.);
    }

    for (i = 0; i < 256; i++) {
      for (j = 0; j < 256; j++) {
        ct[i*256 + j].red = (int)((int)ct[j + 65280].red * i / 65535.);
        ct[i*256 + j].green = (int)((int)ct[j + 65280].green * i / 65535.);
        ct[i*256 + j].blue = (int)((int)ct[j + 65280].blue * i / 65535.);
      }
    }
    *nc = 65535;
    break;

  case PHASE:
    for (j = 0; j < 65536; j++) {			/* RMG colors */
      ct[j].red = 3 * MIN(MIN(j, 21845), 65535 - j) / 256;
      ct[j].green = 3 * MIN(MAX(j - 21845, 21845 - j), 21845) / 256;
      ct[j].blue = 3 * MIN(MAX(j - 43690, 43690 - j), 21845) / 256;
    }
    *nc = 65535;
    break;

  default:
    printf("\nWARNING: ct_16bit: invalid color table selected: %d\n\n", ct_num);
    return ( -1);
    break;
  }
  return (0);
}

int rd_image(char *fn, unsigned char *image, int im_format, COLOR *ct, int *nc, int *nbits, int *width, int *height) {
  /*
        read a SUN raster or BMP format image. This routine allocates memory for the image buffer, The color
        table is returned in array ct1 with nc elements.
        17-Jan-2000 clw copyright Gamma Remote Sensing AG

    	fn		file name for output image, type defined by extension *.ras, *.bmp etc...
    	image		pointer to image array
  	im_format	image format (see display.h for supported formats)
  	ct		pointer to color table array
  	nc		pointer to location containing number of color table entries
    	ct_type		integer specifying color table for 8 bit or 16 bit images, 0=greyscale, 1=RMG, ...;
    	nbits		number of bits per sample, 8, 16, 24
    	width		width of image
  	height		number of lines
  */
  int i;
  int nc1;					/* number of colors in the color tables */
  int ofb;
  int width1, height1;
  int nbits1;

  RASTER_HDR *r1;		/* SUN raster file headers */
  BITMAPFILEHEADER *bfh1;	/* BMP format image headers */
  BITMAPINFOHEADER *bih1;
  FILE *fnf;

  switch (im_format) {
  case SUN:
    init_rdras(fn, &fnf, &r1, &ct, &nc1, &nbits1, &width1, &height1);
    break;
  case BMP:
    init_rdbmp(fn, &fnf, &bfh1, &bih1, &ct, &nc1, &nbits1, &width1, &height1);
    break;
  default:
    fprintf(stderr, "\nERROR rd_image: unsupported image format: %d\n\n", im_format);
    exit( -1);
  }

  *width = width1;
  *height = height1;
  *nbits = nbits1;
  *nc = nc1;

  switch (nbits1) {		/* allocate image buffer */
  case 8:
    image = (unsigned char *)malloc(sizeof(char) * width1 * height1);
    break;
  case 16:
    image = (unsigned char *)malloc(sizeof(unsigned short) * width1 * height1);
    break;
  case 24:
    image = (unsigned char *)malloc(sizeof(COLOR) * width1 * height1);
    break;
  default:
    fprintf(stderr, "\nERROR rd_image: invalid number of bits/pixel: %d\n", nbits1);
    exit( -1);
  }

  if (image == NULL) {
    fprintf(stderr, "\nERROR rd_image: memory allocation error for image data, width: %d   height:%d\n", width1, height1);
  }

  for (i = 0; i < height1; i++) {
    ofb = width1 * i * nbits1 / 8;
    switch (im_format) {
    case SUN:
      rdras(fnf, r1, i, image + ofb);
      break;
    case BMP:
      rdbmp(fnf, bfh1, bih1, i, image + ofb);
      break;
    default:
      fprintf(stderr, "\nERROR rd_image: unsupported image format: %d\n\n", im_format);
      exit( -1);
    }
  }
  fclose(fnf);		/* close the file */
  return (0);
}

int wr_image(char *fn, unsigned char *image, int im_format, COLOR *ct, int nc, int nbits, int width, int height) {
  /*
       This subroutine writes an image file stored in array image. 8, 16, or 24 bits/image are
       supported.
       v1.1 updated byte offset for > 2GB 30-May-2006
       clw copyright Gamma Remote Sensing AG

    	fn		file name for output image;
    	image		pointer to image array
  	im_format	image format (see display.h for supported formats)
  	ct		pointer to color table array
  	nc		pointer to location containing number of color table entries
    	nbits		number of bits per sample, 8, 16, 24
    	width		width of image
  	height		number of lines
  */
  int i;
  off_t ofb;			/* offset number of bytes */

  RASTER_HDR *r1;		/* SUN raster file headers */
  BITMAPFILEHEADER *bfh1;	/* BMP format image headers */
  BITMAPINFOHEADER *bih1;
  FILE *fnf;

  switch (nbits) {
  case 8:
    switch (im_format) {
    case SUN:
      init_wrras(fn, &fnf, &r1, ct, nc, nbits, width, height);
      break;
    case BMP:
      init_wrbmp(fn, &fnf, &bfh1, &bih1, ct, nc, nbits, width, height);
      break;
    default:
      fprintf(stderr, "\nERROR wr_image: unsupported image format in 8 bit mode: %d\n\n", im_format);
      exit( -1);
    }
    break;

  case 16:
    switch (im_format) {
    case BMP:
      init_wrbmp(fn, &fnf, &bfh1, &bih1, ct, nc, nbits, width, height);
      break;
    default:
      fprintf(stderr, "\nERROR wr_image:: unsupported image format in 16 bit mode: %d\n\n", im_format);
      exit( -1);
    }
    break;

  case 24:
    ct = NULL;
    nc = 0;
    switch (im_format) {
    case SUN:
      init_wrras(fn, &fnf, &r1, ct, nc, nbits, width, height);
      break;
    case BMP:
      init_wrbmp(fn, &fnf, &bfh1, &bih1, ct, nc, nbits, width, height);
      break;
    default:
      fprintf(stderr, "\nERROR wr_image: unsupported image format in 24 bit mode: %d\n\n", im_format);
      exit( -1);
    }
    break;

  default:
    fprintf(stderr, "\nERROR wr_image: unsupported number of bits per sample: %d\n\n", nbits);
    exit( -1);
  }

  for (i = 0; i < height; i++) {
    ofb = (off_t)width * i * (nbits / 8);
    switch (im_format) {
    case SUN:
      wrras(fnf, r1, i, image + ofb);
      break;
    case BMP:
      wrbmp(fnf, bfh1, bih1, i, image + ofb);
      break;
    default:
      fprintf(stderr, "\nERROR: unsupported image format: %d\n\n", im_format);
      exit( -1);
    }
  }
  fclose(fnf);		/* close the file */
  return (0);
}


/*
 * Concave Polygon Scan Conversion
 * by Paul Heckbert
 * from "Graphics Gems", Academic Press, 1990
 *
 * concave: scan convert nvert-sided concave non-simple polygon with vertices at
 * (point[i].x, point[i].y) for i in [0..nvert-1] within the window win by
 * calling spanproc for each visible span of pixels.
 * Polygon can be clockwise or counterclockwise.
 * Algorithm does uniform point sampling at pixel centers.
 * Inside-outside test done by Jordan's rule: a point is considered inside if
 * an emanating ray intersects the polygon an odd number of times.
 * spanproc should fill in pixels from xl to xr inclusive on scanline y,
 * e.g:
 *	spanproc(y, xl, xr)
 *	int y, xl, xr;
 *	{
 *	    int x;
 *	    for (x=xl; x<=xr; x++)
 *		pixel_write(x, y, pixelvalue);
 *	}
 *
 *  Paul Heckbert	30 June 81, 18 Dec 89
*/

static int n;			/* number of vertices */
static Point2 *pt;		/* vertices */

static int nact;		/* number of active edges */
static Edge *active;		/* active edge list:edges crossing scanline y */

static void cdelete(int i);		/* remove edge i from active list */
static void cinsert(int i, int y);	/* append edge i to end of active list */

int compare_ind(const void *a, const void *b);
int compare_active(const void *a, const void *b);

void concave(int nvert, Point2 *point, Win *win, Point2 *ptx, int *np) {
  int k, y0, y1, y, i, j, xl, xr, l;
  int *ind;			/* list of vertex indices, sorted by pt[ind[j]].y */

  n = nvert;
  pt = point;
  if (n <= 0)
    return ;
  ALLOC(ind, int, n);
  ALLOC(active, Edge, n);

  /* create y-sorted array of indices ind[k] into vertex list */
  for (k = 0; k < n; k++) ind[k] = k;
  qsort(ind, n, sizeof ind[0], compare_ind);	/* sort ind by pt[ind[k]].y */

  nact = 0;				/* start with empty active list */
  k = 0;				/* ind[k] is next vertex to process */
  y0 = (int) MAX(win->y0, ceil(pt[ind[0]].y - .5));		/* ymin of polygon */
  y1 = (int) MIN(win->y1, floor(pt[ind[n - 1]].y - .5));	/* ymax of polygon */

  for (y = y0; y <= y1; y++) {		/* step through scanlines */
    /* scanline y is at y+.5 in continuous coordinates */

    /* check vertices between previous scanline and current one, if any */
    for (; k < n && pt[ind[k]].y <= y + .5; k++) {
      /* to simplify, if pt.y=y+.5, pretend it's above */
      /* invariant: y-.5 < pt[i].y <= y+.5 */
      i = ind[k];
      /*
       * insert or delete edges before and after vertex i (i-1 to i,
       * and i to i+1) from active list if they cross scanline y
       */
      j = i > 0 ? i - 1 : n - 1;	/* vertex previous to i */
      if (pt[j].y <= y - .5)		/* old edge, remove from active list */
        cdelete(j);
      else if (pt[j].y > y + .5)	/* new edge, add to active list */
        cinsert(j, y);
      j = i < n - 1 ? i + 1 : 0;	/* vertex next after i */
      if (pt[j].y <= y - .5)		/* old edge, remove from active list */
        cdelete(i);
      else if (pt[j].y > y + .5)	/* new edge, add to active list */
        cinsert(i, y);
    }

    /* sort active edge list by active[j].x */
    qsort(active, nact, sizeof active[0], compare_active);

    /* draw horizontal segments for scanline y */
    for (j = 0; j < nact; j += 2) {			/* draw horizontal segments */
      /* span between j & j+1 is inside, span between j+1 & j+2 is outside */
      xl = (int) ceil(active[j].x - .5);		/* left end of span */
      if (xl < win->x0) xl = win->x0;
      xr = (int) floor(active[j + 1].x - .5);		/* right end of span */
      if (xr > win->x1) xr = win->x1;

#ifdef DRAW_SPAN
      if (xl <= xr) (*spanproc)(y, xl, xr);		/* draw pixels in span */
#endif
      if (xl <= xr) {
        for (l = xl; l <= xr; l++) {
          ptx[*np].x = (double)l;
          ptx[*np].y = (double)y;
          (*np)++;
        }
      }
      active[j].x += active[j].dx;	/* increment edge coords */
      active[j + 1].x += active[j + 1].dx;
    }
  }
}

static void cdelete(int i) {		/* remove edge i from active list */
  int j;

  for (j = 0; j < nact && active[j].i != i; j++)
    ;
  if (j >= nact)
    return ;	/* edge not in active list; happens at win->y0*/
  nact--;
  memcpy((char *)&active[j], (char *)&active[j + 1], (nact - j)*sizeof active[0]);
}

static void cinsert(int i, int y) {		/* append edge i to end of active list */
  int j;
  double dx;
  Point2 *p, *q;

  j = i < n - 1 ? i + 1 : 0;
  if (pt[i].y < pt[j].y) {
    p = &pt[i];
    q = &pt[j];
  } else	{
    p = &pt[j];
    q = &pt[i];
  }
  /* initialize x position at intersection of edge with scanline y */
  active[nact].dx = dx = (q->x - p->x) / (q->y - p->y);
  active[nact].x = dx * (y + .5 - p->y) + p->x;
  active[nact].i = i;
  nact++;
}

/* comparison routines for qsort */
int compare_ind(const void *a, const void *b) {
  return pt[*(int *)a].y <= pt[*(int *)b].y ? -1 : 1;
}

int compare_active(const void *a, const void *b) {
  return ((Edge *)a)->x <= ((Edge *)b)->x ? -1 : 1;
}

double value(double n1, double n2, double hue) {
  double va;
  if (hue > 360.)
    hue -= 360.;
  if (hue < 0.)
    hue += 360.;
  if (hue < 60.)
    va = n1 + (n2 - n1) * hue / 60.0;
  else if (hue < 180.)
    va = n2;
  else if (hue < 240.)
    va = n1 + (n2 - n1) * (240. - hue) / 60.;
  else
    va = n1;
  return (va);
}

int hls_to_rgb(double hue, double lightness, double saturation, double *red, double *green, double *blue) {

  double m1, m2;

  if (lightness <= 0.5)
    m2 = lightness * (1. + saturation);
  else
    m2 = lightness + saturation - lightness * saturation;

  m1 = 2. * lightness - m2;

  if (saturation <= 0.0) {
    *red = (double)255. * lightness;
    *green = (double)255. * lightness;
    *blue = (double)255. * lightness;
  } else {
    *red = (double)255. * value(m1, m2, hue + 120.);
    *green = (double)255. * value(m1, m2, hue);
    *blue = (double)255. * value(m1, m2, hue - 120.);
  }
  return (1);

}

int hsv_to_rgb(double hue, double value, double saturation, double *red, double *green, double *blue) {
  int i = 0;
  double h, f, p = 0., q = 0., t = 0.;

  if (saturation <= 0.0) {
    *red = (double)255. * value;
    *green = (double)255. * value;
    *blue = (double)255. * value;
  } else {
    if (hue == 360.0) hue = 0.0;
    h = hue / 60.;
    i = (int)h;
    f = h - (double)i;
    p = value * (1. - saturation);
    q = value * (1. - saturation * f);
    t = value * (1. - (saturation * (1. - f)));
  }


  switch (i) {
  case 0:
    *red = (double)255. * value;
    *green = (double)255. * t;
    *blue = (double)255. * p;
    break;
  case 1:
    *red = (double)255. * q;
    *green = (double)255. * value;
    *blue = (double)255. * t;
    break;
  case 2:
    *red = (double)255. * p;
    *green = (double)255. * value;
    *blue = (double)255. * t;
    break;
  case 3:
    *red = (double)255. * p;
    *green = (double)255. * q;
    *blue = (double)255. * value;
    break;
  case 4:
    *red = (double)255. * t;
    *green = (double)255. * p;
    *blue = (double)255. * value;
    break;
  case 5:
    *red = (double)255. * value;
    *green = (double)255. * p;
    *blue = (double)255. * q;
    break;
  default:
    printf("\nWARNING: ct_8bit: invalid HSV case: %d \n\n", i);
    return ( -1);
    break;
  }

  return (1);

}

int rgb_to_hsv(double red, double green, double blue, double *hue, double *value, double *saturation) {
  /*
    convert from RGB color system to Hue Saturation Value system  ref. Fundementals of
    Computer Graphics. Foley + Van Dam pp. 615 4-Dec-2001
  */
  double maxc, minc;
  double h, s, v;
  double rc, gc, bc;
  double delta;

  maxc = MAX(red, green);
  maxc = MAX(maxc, blue);

  minc = MIN(red, green);
  minc = MIN(minc, blue);

  v = maxc;		/* value */

  if (maxc > 0.0)
    s = (maxc - minc) / maxc;	/* saturation */
  else
    s = 0.0;

  h = 0.0;		/* hue undefined if s = 0.0 */
  if (s != 0.0) {
    delta = maxc - minc;
    rc = (maxc - red) / delta;
    gc = (maxc - green) / delta;
    bc = (maxc - blue) / delta;

    if (red == maxc)
      h = bc - gc;
    if (green == maxc)
      h = 2.0 + rc - bc;
    if (blue == maxc)
      h = 4.0 + gc - rc;

    h = h * 60.0;
    if (h < 0.0)
      h += 360.0;
  }
  *saturation = s;
  *hue = h;
  *value = v;
  return 1;
}

