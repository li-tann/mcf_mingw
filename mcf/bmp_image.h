#ifndef _bmp_h
#define _bmp_h

#define BI_RGB      0L
#define BI_RLE8     1L			/* run length encoding */
#define BI_RLE4     2L

/* all data are always in LITTLE ENDIAN byte order (INTEL compliant) */

typedef struct{
  char bfType[2];	
  unsigned int bfSize; 
  short bfReserved1; 
  short bfReserved2; 
  int bfOffBits; 
}BITMAPFILEHEADER;

/* The BITMAPFILEHEADER data structure contains the following fields:

Field           Description
bfType          Specifies the type of file. It must be BM.
bfSize          Specifies the size in DWORDs of the file. 
bfReserved1     Is reserved and must be set to zero. 
bfReserved2     Is reserved and must be set to zero. 
bfOffBits       Specifies in bytes the offset from the BITMAPFILEHEADER 
                of the actual bitmap in the file. 

*/

typedef struct{
  int biSize; 
  int biWidth;
  int biHeight;
  short biPlanes;
  short biBitCount;
  int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  int biClrUsed;
  int biClrImportant;
}BITMAPINFOHEADER;

/* The BITMAPINFOHEADER structure has the following fields:

Field           Description
biSize          Specifies the number of bytes required by the 
                BITMAPINFOHEADER structure. 
biWidth         Specifies the width of the bitmap in pixels. 
biHeight        Specifies the height of the bitmap in pixels. 
biPlanes        Specifies the number of planes for the target device and
                must be set to 1. 
biBitCount      Specifies the number of bits per pixel. This value must 
                be 1, 4, 8, or 24. 
biCompression   Specifies the type of compression for a compressed      
                bitmap. It can be one of the following values:.
                Value           Meaning
                BI_RGB          Specifies that the bitmap is not 
                                compressed.
                BI_RLE8         Specifies a run-length encoded format 
                                for bitmaps with 8 bits per pixel. The 
                                compression format is a two-byte 
                                format consisting of a count byte 
                                followed by a byte containing a color 
                                index. See the following 'Comments' 
                                section for more information.
                BI_RLE4         Specifies a run-length encoded format 
                                for bitmaps with 4 bits per pixel. The 
                                compression format is a two-byte 
                                format consisting of a count byte 
                                followed by two word-length color 
                                indexes. See the following 'Comments' 
                                section for more information.
biSizeImage     Specifies the size in bytes of the image. 
biXPelsPerMeter Specifies the horizontal resolution in pixels per meter                         of the target device for the bitmap. An application can                         use this value to select a bitmap from a resource group                         that best matches the characteristics of the current                            device. biYPelsPerMeter Specifies the vertical                                  resolution in pixels per meter of the target device for                         the bitmap. 
biClrUsed       Specifies the number of color indexes in the color table 
                actually used by the bitmap. If this value is 0, the    
                bitmap uses the maximum number of colors corresponding  
                to the value of the biBitCount field. See the           
                description of the BITMAPINFO data structure earlier in         
                this chapter for more information on the maximum sizes 
                of the color table. If biClrUsed is nonzero, then the 
                biClrUsed field specifies the actual number of colors 
                which the graphics engine or device driver will access 
                if the biBitCount field is less than 24. If the         
                biBitCount field is set to 24, the biClrUsed field      
                specifies the size of the reference color table used to 
                optimize performance of Windows color palettes.
                If the bitmap is a 'packed' bitmap (that is, a bitmap in 
                which the bitmap array immediately follows the  
                BITMAPINFO header and which is referenced by a single   
                pointer), the biClrUsed field must be set to 0 or to the 
                actual size of the color table. 
biClrImportant  Specifies the number of color indexes that are considered 
                important for displaying the bitmap. If this value is 0,        
                then all colors are important. 
*/

typedef struct{
   unsigned char    blue;
   unsigned char    green;
   unsigned char    red;
   unsigned char    reserved;
} BMP_RGBQUAD;

/*
The RGBQUAD structure contains the following fields:

Field           Description
rgbBlue         Specifies the intensity of blue in the color. 
rgbGreen        Specifies the intensity of green in the color. 
rgbRed          Specifies the intensity of red in the color. 
rgbReserved     Is not used and must be set to zero. 

*/
#endif
