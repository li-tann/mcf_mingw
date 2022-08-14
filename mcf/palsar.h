/*  Copyright Gamma Remote Sensing AG 2006 
    v1.0 clw 11-Jun-2006
    v1.1 clw 06-Oct-2006	Added parameters for ERSDAC data
    v1.2 clw 02-May-2007	Added mode and polarization paramters extracted from Aux. data
*/

#define HH	0		/* polarization definitions */
#define HV	1
#define VH	2
#define VV	4
/*
  The SAR echo consists of a 412 or 292 byte header + 2*echo samples + 2*right fill samples 
*/
#define SDSH_SIZE_EOC		720	/* size of SAR data file descriptor record in bytes for EOC data */
#define SDSH_SIZE_ERSDAC	16252	/* size of SAR data file descriptor record in bytes for ERSDAC data*/
#define NSTATE_CEOS_EOC		28	/* number of state vectors in the CEOS leader file for EOC data */
#define NSTATE_CEOS_ERSDAC	15	/* number of state vectors in the CEOS leader file for ERSDAC data */
#define SLH_SIZE_EOC		412	/* signal data line header size in bytes EOC data */
#define SLH_SIZE_ERSDAC		292	/* signal data line header size in bytes ERSDAC data */

typedef struct{
  int mode;		/* radar mode 0:init, 1:standby2, 2:standby3, 3:standby, 4:calib, 5:obs standby (calib), 6: obs from aux. data */
  int obs_mode;		/* OBS mode: 0:Stripmap/Direct downlink, 1:ScanSAR, 2:Polarimetry, from aux. data */
  int Mpps;		/* 1Mpps timer from observation auxiliary data */
  int Rx_pol;           /* Rx polarization setting: 0:H, 1:V, 2:H+V, from aux. data */
  int Tx_pol;           /* Tx polarization setting: 0:V, 1:H, from aux. data (opposite of Rx_pol definition of H and V!!) */
  int lnum;		/* line number */
  int ns;		/* SAR echo samples */
  int rfill;		/* right fill samples */
  int year;		/* acquisition year */
  int doy;		/* day of year */		
  int msod;		/* millisec of day */
  int channel;		/* SAR channel (1 --> 4) */
  int txpol;		/* 0: H, 1: V */
  int rxpol;		/* 0: H  1: V */
  int SS_mode;		/* SCANSAR mode 1-5 expect Wide Observation mode = always 0 */
  int prf;		/* milli-Hz */
  int tdur;		/* chirp duration nano-sec */
  int chirp_rate;	/* chirp rate Hz/micro-sec */
  int rx_gain;		/* receiver gain dB */
  int zero_flg;		/* missing line flag 0: data correct  1: missing line */	
  int rho;		/* distance to first sample (m) */
  int tdel;		/* data record window (s) */ 
  int Rx_width2;	/* receive gate width 2 */
  int Rx_start2;	/* receiver gate start 2 */
  int pri;		/* radar PRI in microsec */
}PALSAR_LINE_HDR;

typedef struct{		/* parameters from the 720 byte header at the start of the IMG Level 1 (raw data) */
  int nlines;		/* number of lines in the SAR data set */
  int reclen;		/* SAR data record length (including fill and hdar) */
  int nbsar;		/* number of bytes of SAR data  (reclen - 412) */
}PALSAR_DATA_SET_HDR;

