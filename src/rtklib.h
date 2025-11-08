/*------------------------------------------------------------------------------
 * rtklib.h : RTKLIB constants, types and function prototypes
 *
 * author  : sujinglan
 *-----------------------------------------------------------------------------*/
#ifndef RTKLIB_H
#define RTKLIB_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdint.h>

#ifdef WIN32
#include <windows.h>
#define thread_t    HANDLE
#define lock_t      CRITICAL_SECTION
#define initlock(f) InitializeCriticalSection(f)
#define lock(f)     EnterCriticalSection(f)
#define unlock(f)   LeaveCriticalSection(f)
#else
#include <pthread.h>
#define thread_t    pthread_t
#define lock_t      pthread_mutex_t
#define initlock(f) pthread_mutex_init(f,NULL)
#define lock(f)     pthread_mutex_lock(f)
#define unlock(f)   pthread_mutex_unlock(f)
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* constants -----------------------------------------------------------------*/
#define PI          3.1415926535897932  /* pi */
#define D2R         (PI/180.0)          /* deg to rad */
#define R2D         (180.0/PI)          /* rad to deg */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define SC2RAD      3.1415926535898     /* semi-circle to radian (IS-GPS) */
#define AU          149597870691.0      /* 1 AU (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */
#define OMGE        7.2921151467E-5     /* earth angular velocity (IS-GPS) (rad/s) */

#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#define E1_WGS84    0.0066943799901413156
#define E_SQR       0.00669437999014    /* sqr of linear eccentricity of the ellipsoid */
#define HION        350000.0            /* ionosphere height (m) */
#define AS2R        (D2R/3600.0)        /* arc sec to radian */

#define REARRANGEOBS_BYFREQ
#define TRACE
#define ENAGLO
#define ENAGAL
#define ENACMP
#define ENAQZS

#define MAXFREQ     7                   /* max NFREQ */
#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */
#define FREQ2       1.22760E9           /* L2         frequency (Hz) */
#define FREQ5       1.17645E9           /* L5/E5a/B2a frequency (Hz) */
#define FREQ6       1.27875E9           /* E6/L6  frequency (Hz) */
#define FREQ7       1.20714E9           /* E5b    frequency (Hz) */
#define FREQ8       1.191795E9          /* E5a+b  frequency (Hz) */
#define FREQ9       2.492028E9          /* S      frequency (Hz) */
#define FREQ1_GLO   1.60200E9           /* GLONASS G1 base frequency (Hz) */
#define DFRQ1_GLO   0.56250E6           /* GLONASS G1 bias frequency (Hz/n) */
#define FREQ2_GLO   1.24600E9           /* GLONASS G2 base frequency (Hz) */
#define DFRQ2_GLO   0.43750E6           /* GLONASS G2 bias frequency (Hz/n) */
#define FREQ3_GLO   1.202025E9          /* GLONASS G3 frequency (Hz) */
#define FREQ1a_GLO  1.600995E9          /* GLONASS G1a frequency (Hz) */
#define FREQ2a_GLO  1.248060E9          /* GLONASS G2a frequency (Hz) */
#define FREQ1_CMP   1.561098E9          /* BDS B1I     frequency (Hz) */
#define FREQ2_CMP   1.20714E9           /* BDS B2I/B2b frequency (Hz) */
#define FREQ3_CMP   1.26852E9           /* BDS B3      frequency (Hz) */

#define SYS_NONE    0x00                /* navigation system: none */
#define SYS_GPS     0x01                /* navigation system: GPS */
#define SYS_SBS     0x02                /* navigation system: SBAS */
#define SYS_GLO     0x04                /* navigation system: GLONASS */
#define SYS_GAL     0x08                /* navigation system: Galileo */
#define SYS_QZS     0x10                /* navigation system: QZSS */
#define SYS_CMP     0x20                /* navigation system: BeiDou */
#define SYS_IRN     0x40                /* navigation system: IRNS */
#define SYS_LEO     0x80                /* navigation system: LEO */
#define SYS_ALL     0xFF                /* navigation system: all */

#define TSYS_GPS    0                   /* time system: GPS time */
#define TSYS_UTC    1                   /* time system: UTC */
#define TSYS_GLO    2                   /* time system: GLONASS time */
#define TSYS_GAL    3                   /* time system: Galileo time */
#define TSYS_QZS    4                   /* time system: QZSS time */
#define TSYS_CMP    5                   /* time system: BeiDou time */
#define TSYS_IRN    6                   /* time system: IRNSS time */

#ifndef NFREQ
#define NFREQ       3                   /* number of carrier frequencies */
#endif
#define NFREQGLO    2                   /* number of carrier frequencies of GLONASS */

#ifndef NEXOBS
#define NEXOBS      3                   /* number of extended obs codes */
#endif

#define SNR_UNIT    0.001               /* SNR unit (dBHz) */

#define MINPRNGPS   1                   /* min satellite PRN number of GPS */
#define MAXPRNGPS   32                  /* max satellite PRN number of GPS */
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) /* number of GPS satellites */
#define NSYSGPS     1

#ifdef ANDROID
#ifdef  ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO1  27
#define MAXPRNGLO   (MAXPRNGLO1+15)     /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif
#else
#ifdef  ENAGLO
#define MINPRNGLO   1                   /* min satellite slot number of GLONASS */
#define MAXPRNGLO   27                  /* max satellite slot number of GLONASS */
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) /* number of GLONASS satellites */
#define NSYSGLO     1
#else
#define MINPRNGLO   0
#define MAXPRNGLO   0
#define NSATGLO     0
#define NSYSGLO     0
#endif
#endif

#ifdef ENAGAL
#define MINPRNGAL   1                   /* min satellite PRN number of Galileo */
#define MAXPRNGAL   36                  /* max satellite PRN number of Galileo */
#define NSATGAL    (MAXPRNGAL-MINPRNGAL+1) /* number of Galileo satellites */
#define NSYSGAL     1
#else
#define MINPRNGAL   0
#define MAXPRNGAL   0
#define NSATGAL     0
#define NSYSGAL     0
#endif
#ifdef ENAQZS
#define MINPRNQZS   193                 /* min satellite PRN number of QZSS */
#define MAXPRNQZS   202                 /* max satellite PRN number of QZSS */
#define MINPRNQZS_S 183                 /* min satellite PRN number of QZSS L1S */
#define MAXPRNQZS_S 191                 /* max satellite PRN number of QZSS L1S */
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) /* number of QZSS satellites */
#define NSYSQZS     1
#else
#define MINPRNQZS   0
#define MAXPRNQZS   0
#define MINPRNQZS_S 0
#define MAXPRNQZS_S 0
#define NSATQZS     0
#define NSYSQZS     0
#endif
#ifdef ENACMP
#define MINPRNCMP   1                   /* min satellite sat number of BeiDou */
#define MAXPRNCMP   63                  /* max satellite sat number of BeiDou */
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) /* number of BeiDou satellites */
#define NSYSCMP     1
#else
#define MINPRNCMP   0
#define MAXPRNCMP   0
#define NSATCMP     0
#define NSYSCMP     0
#endif
#ifdef ENAIRN
#define MINPRNIRN   1                   /* min satellite sat number of IRNSS */
#define MAXPRNIRN   14                  /* max satellite sat number of IRNSS */
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1) /* number of IRNSS satellites */
#define NSYSIRN     1
#else
#define MINPRNIRN   0
#define MAXPRNIRN   0
#define NSATIRN     0
#define NSYSIRN     0
#endif
#ifdef ENALEO
#define MINPRNLEO   1                   /* min satellite sat number of LEO */
#define MAXPRNLEO   10                  /* max satellite sat number of LEO */
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) /* number of LEO satellites */
#define NSYSLEO     1
#else
#define MINPRNLEO   0
#define MAXPRNLEO   0
#define NSATLEO     0
#define NSYSLEO     0
#endif

#define MINPRNSBS   120                 /* min satellite PRN number of SBAS */
#define MAXPRNSBS   158                 /* max satellite PRN number of SBAS */
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) /* number of SBAS satellites */

#define MAXSAT      (NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATIRN+NSATSBS+NSATLEO)
                                        /* max satellite number (1 to MAXSAT) */
#ifndef MAXOBS
#define MAXOBS      96                  /* max number of obs in an epoch */
#endif
#define MAXRCV      64                  /* max receiver number (1 to MAXRCV) */
#define MAXOBSTYPE  64                  /* max number of obs type in RINEX */
#ifdef OBS_100HZ
#define DTTOL       0.005               /* tolerance of time difference (s) */
#else
#define DTTOL       0.025               /* tolerance of time difference (s) */
#endif
#define MAXDTOE     7200.0              /* max time difference to GPS Toe (s) */
#define MAXDTOE_QZS 7200.0              /* max time difference to QZSS Toe (s) */
#define MAXDTOE_GAL 14400.0             /* max time difference to Galileo Toe (s) */
#define MAXDTOE_CMP 21600.0             /* max time difference to BeiDou Toe (s) */
#define MAXDTOE_GLO 1800.0              /* max time difference to GLONASS Toe (s) */
#define MAXDTOE_IRN 7200.0              /* max time difference to IRNSS Toe (s) */
#define MAXDTOE_SBS 360.0               /* max time difference to SBAS Toe (s) */
#define MAXDTOE_S   86400.0             /* max time difference to ephem toe (s) for other */

#define INT_SWAP_TRAC 86400.0           /* swap interval of trace file (s) */

#define MAXCOMMENT  10                  /* max number of RINEX comments */
#define MAXEXFILE   1024                /* max number of expanded files */
#define MAXSBSAGEF  30.0                /* max age of SBAS fast correction (s) */
#define MAXSBSAGEL  1800.0              /* max age of SBAS long term corr (s) */
#define MAXSBSURA   8                   /* max URA of SBAS satellite */
#define MAXBAND     10                  /* max SBAS band of IGP */
#define MAXNIGP     201                 /* max number of IGP in SBAS band */
#define MAXNGEO     4                   /* max number of GEO satellites */
#define MAXSOLMSG   8191                /* max length of solution message */
#define MAXRAWLEN   8191                /* max length of receiver raw message */
#define MAXANT      64                  /* max length of station name/antenna type */
#define MAXLEAPS    64                  /* max number of leap seconds table */

#define OBSTYPE_PR  0x01                /* observation type: pseudorange */
#define OBSTYPE_CP  0x02                /* observation type: carrier-phase */
#define OBSTYPE_DOP 0x04                /* observation type: doppler-freq */
#define OBSTYPE_SNR 0x08                /* observation type: SNR */
#define OBSTYPE_ALL 0xFF                /* observation type: all */

#define FREQTYPE_L1 0x01                /* frequency type: L1/E1 */
#define FREQTYPE_L2 0x02                /* frequency type: L2/B1 */
#define FREQTYPE_L5 0x04                /* frequency type: L5/E5a/L3 */
#define FREQTYPE_L6 0x08                /* frequency type: E6/LEX/B3 */
#define FREQTYPE_L7 0x10                /* frequency type: E5b/B2 */
#define FREQTYPE_L8 0x20                /* frequency type: E5(a+b) */
#define FREQTYPE_L9 0x40                /* frequency type: S */
#define FREQTYPE_ALL 0xFF               /* frequency type: all */

#define CODE_NONE   0                   /* obs code: none or unknown */
#define CODE_L1C    1                   /* obs code: L1C/A,G1C/A,E1C (GPS,GLO,GAL,QZS,SBS) */
#define CODE_L1P    2                   /* obs code: L1P,G1P,B1P (GPS,GLO,BDS) */
#define CODE_L1W    3                   /* obs code: L1 Z-track (GPS) */
#define CODE_L1Y    4                   /* obs code: L1Y        (GPS) */
#define CODE_L1M    5                   /* obs code: L1M        (GPS) */
#define CODE_L1N    6                   /* obs code: L1codeless,B1codeless (GPS,BDS) */
#define CODE_L1S    7                   /* obs code: L1C(D)     (GPS,QZS) */
#define CODE_L1L    8                   /* obs code: L1C(P)     (GPS,QZS) */
#define CODE_L1E    9                   /* (not used) */
#define CODE_L1A    10                  /* obs code: E1A,B1A    (GAL,BDS) */
#define CODE_L1B    11                  /* obs code: E1B        (GAL) */
#define CODE_L1X    12                  /* obs code: E1B+C,L1C(D+P),B1D+P (GAL,QZS,BDS) */
#define CODE_L1Z    13                  /* obs code: E1A+B+C,L1S (GAL,QZS) */
#define CODE_L2C    14                  /* obs code: L2C/A,G1C/A (GPS,GLO) */
#define CODE_L2D    15                  /* obs code: L2 L1C/A-(P2-P1) (GPS) */
#define CODE_L2S    16                  /* obs code: L2C(M)     (GPS,QZS) */
#define CODE_L2L    17                  /* obs code: L2C(L)     (GPS,QZS) */
#define CODE_L2X    18                  /* obs code: L2C(M+L),B1_2I+Q (GPS,QZS,BDS) */
#define CODE_L2P    19                  /* obs code: L2P,G2P    (GPS,GLO) */
#define CODE_L2W    20                  /* obs code: L2 Z-track (GPS) */
#define CODE_L2Y    21                  /* obs code: L2Y        (GPS) */
#define CODE_L2M    22                  /* obs code: L2M        (GPS) */
#define CODE_L2N    23                  /* obs code: L2codeless (GPS) */
#define CODE_L5I    24                  /* obs code: L5I,E5aI   (GPS,GAL,QZS,SBS) */
#define CODE_L5Q    25                  /* obs code: L5Q,E5aQ   (GPS,GAL,QZS,SBS) */
#define CODE_L5X    26                  /* obs code: L5I+Q,E5aI+Q,L5B+C,B2aD+P (GPS,GAL,QZS,IRN,SBS,BDS) */
#define CODE_L7I    27                  /* obs code: E5bI,B2bI  (GAL,BDS) */
#define CODE_L7Q    28                  /* obs code: E5bQ,B2bQ  (GAL,BDS) */
#define CODE_L7X    29                  /* obs code: E5bI+Q,B2bI+Q (GAL,BDS) */
#define CODE_L6A    30                  /* obs code: E6A,B3A    (GAL,BDS) */
#define CODE_L6B    31                  /* obs code: E6B        (GAL) */
#define CODE_L6C    32                  /* obs code: E6C        (GAL) */
#define CODE_L6X    33                  /* obs code: E6B+C,LEXS+L,B3I+Q (GAL,QZS,BDS) */
#define CODE_L6Z    34                  /* obs code: E6A+B+C,L6D+E (GAL,QZS) */
#define CODE_L6S    35                  /* obs code: L6S        (QZS) */
#define CODE_L6L    36                  /* obs code: L6L        (QZS) */
#define CODE_L8I    37                  /* obs code: E5abI      (GAL) */
#define CODE_L8Q    38                  /* obs code: E5abQ      (GAL) */
#define CODE_L8X    39                  /* obs code: E5abI+Q,B2abD+P (GAL,BDS) */
#define CODE_L2I    40                  /* obs code: B1_2I      (BDS) */
#define CODE_L2Q    41                  /* obs code: B1_2Q      (BDS) */
#define CODE_L6I    42                  /* obs code: B3I        (BDS) */
#define CODE_L6Q    43                  /* obs code: B3Q        (BDS) */
#define CODE_L3I    44                  /* obs code: G3I        (GLO) */
#define CODE_L3Q    45                  /* obs code: G3Q        (GLO) */
#define CODE_L3X    46                  /* obs code: G3I+Q      (GLO) */
#define CODE_L1I    47                  /* obs code: B1I        (BDS) (obsolute) */
#define CODE_L1Q    48                  /* obs code: B1Q        (BDS) (obsolute) */
#define CODE_L5A    49                  /* obs code: L5A SPS    (IRN) */
#define CODE_L5B    50                  /* obs code: L5B RS(D)  (IRN) */
#define CODE_L5C    51                  /* obs code: L5C RS(P)  (IRN) */
#define CODE_L9A    52                  /* obs code: SA SPS     (IRN) */
#define CODE_L9B    53                  /* obs code: SB RS(D)   (IRN) */
#define CODE_L9C    54                  /* obs code: SC RS(P)   (IRN) */
#define CODE_L9X    55                  /* obs code: SB+C       (IRN) */
#define CODE_L1D    56                  /* obs code: B1D        (BDS) */
#define CODE_L5D    57                  /* obs code: L5D(L5S),B2aD (QZS,BDS) */
#define CODE_L5P    58                  /* obs code: L5P(L5S),B2aP (QZS,BDS) */
#define CODE_L5Z    59                  /* obs code: L5D+P(L5S) (QZS) */
#define CODE_L6E    60                  /* obs code: L6E        (QZS) */
#define CODE_L7D    61                  /* obs code: B2bD       (BDS) */
#define CODE_L7P    62                  /* obs code: B2bP       (BDS) */
#define CODE_L7Z    63                  /* obs code: B2bD+P     (BDS) */
#define CODE_L8D    64                  /* obs code: B2abD      (BDS) */
#define CODE_L8P    65                  /* obs code: B2abP      (BDS) */
#define CODE_L4A    66                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4B    67                  /* obs code: G1aL1OCd   (GLO) */
#define CODE_L4X    68                  /* obs code: G1al1OCd+p (GLO) */
#define MAXCODE     68                  /* max number of obs code */

#define PMODE_SINGLE 0                  /* positioning mode: single */
#define PMODE_DGPS   1                  /* positioning mode: DGPS/DGNSS */
#define PMODE_KINEMA 2                  /* positioning mode: kinematic */

#define SOLF_LLH    0                   /* solution format: lat/lon/height */
#define SOLF_XYZ    1                   /* solution format: x/y/z-ecef */
#define SOLF_ENU    2                   /* solution format: e/n/u-baseline */
#define SOLF_NMEA   3                   /* solution format: NMEA-183 */
#define SOLF_STAT   4                   /* solution format: solution status */
#define SOLF_GSIF   5                   /* solution format: GSI F1/F2 */

#define SOLQ_NONE   0                   /* solution status: no solution */
#define SOLQ_FIX    1                   /* solution status: fix */
#define SOLQ_FLOAT  2                   /* solution status: float */
#define SOLQ_SBAS   3                   /* solution status: SBAS */
#define SOLQ_DGPS   4                   /* solution status: DGPS/DGNSS */
#define SOLQ_SINGLE 5                   /* solution status: single */
#define SOLQ_PPP    6                   /* solution status: PPP */
#define SOLQ_DR     7                   /* solution status: dead reconing */
#define SOLQ_DEGN   8                   /* solution status: degenerate */
#define SOLQ_CHIP   9                   /* solution status: Android phone chip solution */
#define SOLQ_NLP   10                   /* solution status: Android phone solution from NLP */
#define SOLQ_ASPFT 11                   /* solution status: Android single point position using EKF filter */
#define MAXSOLQ    11                   /* max number of solution status */

#define TIMES_GPST  0                   /* time system: gps time */
#define TIMES_UTC   1                   /* time system: utc */
#define TIMES_JST   2                   /* time system: jst */

#define IONOOPT_OFF  0                  /* ionosphere option: correction off */
#define IONOOPT_BRDC 1                  /* ionosphere option: broadcast model */
#define IONOOPT_SBAS 2                  /* ionosphere option: SBAS model */
#define IONOOPT_IFLC 3                  /* ionosphere option: L1/L2 iono-free LC */
#define IONOOPT_EST  4                  /* ionosphere option: estimation */
#define IONOOPT_TEC  5                  /* ionosphere option: IONEX TEC model */
#define IONOOPT_QZS  6                  /* ionosphere option: QZSS broadcast model */
#define IONOOPT_STEC 8                  /* ionosphere option: SLANT TEC model */

#define TROPOPT_OFF  0                  /* troposphere option: correction off */
#define TROPOPT_SAAS 1                  /* troposphere option: Saastamoinen model */
#define TROPOPT_SBAS 2                  /* troposphere option: SBAS model */
#define TROPOPT_EST  3                  /* troposphere option: ZTD estimation */
#define TROPOPT_ESTG 4                  /* troposphere option: ZTD+grad estimation */
#define TROPOPT_ZTD  5                  /* troposphere option: ZTD correction */

#define EPHOPT_BRDC 0                   /* ephemeris option: broadcast ephemeris */
#define EPHOPT_PREC 1                   /* ephemeris option: precise ephemeris */
#define EPHOPT_SBAS 2                   /* ephemeris option: broadcast + SBAS */
#define EPHOPT_SSRAPC 3                 /* ephemeris option: broadcast + SSR_APC */
#define EPHOPT_SSRCOM 4                 /* ephemeris option: broadcast + SSR_COM */

#define GEOID_EMBEDDED    0             /* geoid model: embedded geoid */
#define GEOID_EGM96_M150  1             /* geoid model: EGM96 15x15" */
#define GEOID_EGM2008_M25 2             /* geoid model: EGM2008 2.5x2.5" */
#define GEOID_EGM2008_M10 3             /* geoid model: EGM2008 1.0x1.0" */
#define GEOID_GSI2000_M15 4             /* geoid model: GSI geoid 2000 1.0x1.5" */
#define GEOID_RAF09       5             /* geoid model: IGN RAF09 for France 1.5"x2" */

#define LLI_SLIP    0x01                /* LLI: cycle-slip */
#define LLI_HALFC   0x02                /* LLI: half-cycle not resovled */
#define LLI_BOCTRK  0x04                /* LLI: boc tracking of mboc signal */
#define LLI_HALFA   0x40                /* LLI: half-cycle added */
#define LLI_HALFS   0x80                /* LLI: half-cycle subtracted */

#define ARMODE_CONT 0                   /* AR mode: continuous */
#define ARMODE_INST 1                   /* AR mode: instantaneous */
#define ARMODE_FIXHOLD 2                /* AR mode: fix and hold */

#define P2_5        0.03125             /* 2^-5 */
#define P2_6        0.015625            /* 2^-6 */
#define P2_11       4.882812500000000E-04 /* 2^-11 */
#define P2_15       3.051757812500000E-05 /* 2^-15 */
#define P2_17       7.629394531250000E-06 /* 2^-17 */
#define P2_19       1.907348632812500E-06 /* 2^-19 */
#define P2_20       9.536743164062500E-07 /* 2^-20 */
#define P2_21       4.768371582031250E-07 /* 2^-21 */
#define P2_23       1.192092895507810E-07 /* 2^-23 */
#define P2_24       5.960464477539063E-08 /* 2^-24 */
#define P2_27       7.450580596923828E-09 /* 2^-27 */
#define P2_29       1.862645149230957E-09 /* 2^-29 */
#define P2_30       9.313225746154785E-10 /* 2^-30 */
#define P2_31       4.656612873077393E-10 /* 2^-31 */
#define P2_32       2.328306436538696E-10 /* 2^-32 */
#define P2_33       1.164153218269348E-10 /* 2^-33 */
#define P2_35       2.910383045673370E-11 /* 2^-35 */
#define P2_38       3.637978807091710E-12 /* 2^-38 */
#define P2_39       1.818989403545856E-12 /* 2^-39 */
#define P2_40       9.094947017729280E-13 /* 2^-40 */
#define P2_43       1.136868377216160E-13 /* 2^-43 */
#define P2_48       3.552713678800501E-15 /* 2^-48 */
#define P2_50       8.881784197001252E-16 /* 2^-50 */
#define P2_55       2.775557561562891E-17 /* 2^-55 */

#define SQR(x)       ((x)*(x))
#define SQRT(x)      ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define ROUND(x)     ((int)floor((x)+0.5))
#define SQV(x,y)     ((x)>(y)?((double)(x)/(double)(y)):0.0)
#define MIN(x,y)     ((x)<=(y)?(x):(y))
#define MAX(x,y)     ((x)>(y)?(x):(y))

#define VAR_POS      SQR(60.0)
#define VAR_VEL      SQR(60.0)
#define VAR_ACC      SQR( 5.0)
#define VAR_CLK      SQR(60.0)
#define VAR_AMB      SQR(60.0)
#define VAR_ZTD      SQR( 0.6)
#define VAR_GRA      SQR(0.01)
#define VAR_DCB      SQR(0.01)
#define VAR_GLO_IFB  SQR(0.6)

#define NP(opt)      (3)
#define NV(opt)      (3)
#define NA(opt)      (3)
#define NC(opt)      (4*NFREQ)
#define ND(opt)      (1)
#define NT(opt)      ((opt)->udppp?1:0)
#define NK(opt)      (1)
#define NB(opt)      ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NFREQ)
#define NR(opt)      (NP(opt)+NV(opt)+NC(opt)+ND(opt))

#define UDFLT_PR     (1<<1)
#define UDFLT_DP     (1<<2)
#define UDFLT_BAD    (1<<3)
#define UDFLT_FAIL   (1<<4)

/* type definitions ----------------------------------------------------------*/
typedef struct {        /* time struct */
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
} gtime_t;

typedef struct {        /* observation data record */
    gtime_t time;       /* receiver sampling time (GPST) */
    uint8_t sat,rcv;    /* satellite/receiver number */
    uint16_t SNR[NFREQ+NEXOBS]; /* signal strength (0.001 dBHz) */
    uint8_t  LLI[NFREQ+NEXOBS]; /* loss of lock indicator */
    uint8_t code[NFREQ+NEXOBS]; /* code indicator (CODE_???) */
    double L[NFREQ+NEXOBS]; /* observation data carrier-phase (cycle) */
    double P[NFREQ+NEXOBS]; /* observation data pseudorange (m) */
    double D[NFREQ+NEXOBS]; /* observation data doppler frequency (Hz) */
} obsd_t;

typedef struct {        /* observation data */
    int n,nmax;         /* number of observation data/allocated */
    obsd_t *data;       /* observation data records */
} obs_t;

typedef struct {        /* earth rotation parameter data type */
    double mjd;         /* mjd (days) */
    double xp,yp;       /* pole offset (rad) */
    double xpr,ypr;     /* pole offset rate (rad/day) */
    double ut1_utc;     /* ut1-utc (s) */
    double lod;         /* length of day (s/day) */
} erpd_t;

typedef struct {        /* earth rotation parameter type */
    int n,nmax;         /* number and max number of data */
    erpd_t *data;       /* earth rotation parameter data */
} erp_t;

typedef struct {        /* antenna parameter type */
    int sat;            /* satellite number (0:receiver) */
    char type[MAXANT];  /* antenna type */
    char code[MAXANT];  /* serial number or satellite code */
    gtime_t ts,te;      /* valid time start and end */
    double off[NFREQ][ 3]; /* phase center offset e/n/u or x/y/z (m) */
    double var[NFREQ][19]; /* phase center variation (m) */
                           /* el=90,85,...,0 or nadir=0,1,2,3,... (deg) */
} pcv_t;

typedef struct {        /* antenna parameters type */
    int n,nmax;         /* number of data/allocated */
    pcv_t *pcv;         /* antenna parameters data */
} pcvs_t;

typedef struct {        /* GPS/QZS/GAL broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode,iodc;      /* IODE,IODC */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    int week;           /* GPS/QZS: gps week, GAL: galileo week */
    int code;           /* GPS/QZS: code on L2 */
                        /* GAL: data source defined as rinex 3.03 */
                        /* BDS: data source (0:unknown,1:B1I,2:B1Q,3:B2I,4:B2Q,5:B3I,6:B3Q) */
    int flag;           /* GPS/QZS: L2 P data flag */
                        /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    int IFnav;          /* for Galileo, IFnav=1: F/nav, IFnav=0: I/Nav */
    gtime_t toe,toc,ttr; /* Toe,Toc,T_trans */
                        /* SV orbit parameters */
    double A,e,i0,OMG0,omg,M0,deln,OMGd,idot;
    double crc,crs,cuc,cus,cic,cis;
    double toes;        /* Toe (s) in week */
    double fit;         /* fit interval (h) */
    double f0,f1,f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[6];      /* group delay parameters */
                        /* GPS/QZS:tgd[0]=TGD */
                        /* GAL:tgd[0]=BGD_E1E5a,tgd[1]=BGD_E1E5b */
                        /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
                        /*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
    double Adot,ndot;   /* Adot,ndot for CNAV */
} eph_t;

typedef struct {        /* GLONASS broadcast ephemeris type */
    int sat;            /* satellite number */
    int iode;           /* IODE (0-6 bit of tb field) */
    int frq;            /* satellite frequency number */
    int svh,sva,age;    /* satellite health, accuracy, age of operation */
    gtime_t toe;        /* epoch of epherides (gpst) */
    gtime_t tof;        /* message frame time (gpst) */
    double pos[3];      /* satellite position (ecef) (m) */
    double vel[3];      /* satellite velocity (ecef) (m/s) */
    double acc[3];      /* satellite acceleration (ecef) (m/s^2) */
    double taun,gamn;   /* SV clock bias (s)/relative freq bias */
    double dtaun;       /* delay between L1 and L2 (s) */
} geph_t;

typedef struct {        /* DGPS/GNSS correction type */
    gtime_t t0;         /* correction time */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range rate correction (RRC) (m/s) */
    int iod;            /* issue of data (IOD) */
    double udre;        /* UDRE */
} dgps_t;

typedef struct {        /* SSR correction type */
    gtime_t t0[6];      /* epoch time (GPST) {eph,clk,hrclk,ura,bias,pbias} */
    double udi[6];      /* SSR update interval (s) */
    int iod[6];         /* iod ssr {eph,clk,hrclk,ura,bias,pbias} */
    int iode;           /* issue of data */
    int iodcrc;         /* issue of data crc for beidou/sbas */
    int ura;            /* URA indicator */
    int refd;           /* sat ref datum (0:ITRF,1:regional) */
    double deph [3];    /* delta orbit {radial,along,cross} (m) */
    double ddeph[3];    /* dot delta orbit {radial,along,cross} (m/s) */
    double dclk [3];    /* delta clock {c0,c1,c2} (m,m/s,m/s^2) */
    double hrclk;       /* high-rate clock corection (m) */
    float  cbias[MAXCODE]; /* code biases (m) */
    double pbias[MAXCODE]; /* phase biases (m) */
    float  stdpb[MAXCODE]; /* std-dev of phase biases (m) */
    double yaw_ang,yaw_rate; /* yaw angle and yaw rate (deg,deg/s) */
    uint8_t update;     /* update flag (0:no update,1:update) */
} ssr_t;

typedef struct {        /* SBAS fast correction type */
    gtime_t t0;         /* time of applicability (TOF) */
    double prc;         /* pseudorange correction (PRC) (m) */
    double rrc;         /* range-rate correction (RRC) (m/s) */
    double dt;          /* range-rate correction delta-time (s) */
    int iodf;           /* IODF (issue of date fast corr) */
    int16_t udre;       /* UDRE+1 */
    int16_t ai;         /* degradation factor indicator */
} sbsfcorr_t;

typedef struct {        /* SBAS long term satellite error correction type */
    gtime_t t0;         /* correction time */
    int iode;           /* IODE (issue of date ephemeris) */
    double dpos[3];     /* delta position (m) (ecef) */
    double dvel[3];     /* delta velocity (m/s) (ecef) */
    double daf0,daf1;   /* delta clock-offset/drift (s,s/s) */
} sbslcorr_t;

typedef struct {        /* SBAS satellite correction type */
    int sat;            /* satellite number */
    sbsfcorr_t fcorr;   /* fast correction */
    sbslcorr_t lcorr;   /* long term correction */
} sbssatp_t;

typedef struct {        /* SBAS satellite corrections type */
    int iodp;           /* IODP (issue of date mask) */
    int nsat;           /* number of satellites */
    int tlat;           /* system latency (s) */
    sbssatp_t sat[MAXSAT]; /* satellite correction */
} sbssat_t;

typedef struct {        /* SBAS ionospheric correction type */
    gtime_t t0;         /* correction time */
    int16_t lat,lon;    /* latitude/longitude (deg) */
    int16_t give;       /* GIVI+1 */
    float delay;        /* vertical delay estimate (m) */
} sbsigp_t;

typedef struct {        /* IGP band type */
    int16_t x;          /* longitude/latitude (deg) */
    const int16_t *y;   /* latitudes/longitudes (deg) */
    uint8_t bits;       /* IGP mask start bit */
    uint8_t bite;       /* IGP mask end bit */
} sbsigpband_t;

typedef struct {        /* SBAS ionospheric corrections type */
    int iodi;           /* IODI (issue of date ionos corr) */
    int nigp;           /* number of igps */
    sbsigp_t igp[MAXNIGP]; /* ionospheric correction */
} sbsion_t;

typedef struct {        /* SBAS ephemeris type */
    int sat;            /* satellite number */
    gtime_t t0;         /* reference epoch time (GPST) */
    gtime_t tof;        /* time of message frame (GPST) */
    int sva;            /* SV accuracy (URA index) */
    int svh;            /* SV health (0:ok) */
    double pos[3];      /* satellite position (m) (ecef) */
    double vel[3];      /* satellite velocity (m/s) (ecef) */
    double acc[3];      /* satellite acceleration (m/s^2) (ecef) */
    double af0,af1;     /* satellite clock-offset/drift (s,s/s) */
} seph_t;

typedef struct {        /* SBAS message type */
    int week,tow;       /* reception time */
    uint8_t prn,rcv;    /* SBAS satellite PRN,receiver number */
    uint8_t msg[29];    /* SBAS message (226bit) padded by 0 */
} sbsmsg_t;

typedef struct {        /* SBAS messages type */
    int n,nmax;         /* number of SBAS messages/allocated */
    sbsmsg_t *msgs;     /* SBAS messages */
} sbs_t;

typedef struct {        /* navigation data type */
    int n,nmax;         /* number of broadcast ephemeris */
    int ng,ngmax;       /* number of glonass ephemeris */
    int ns,nsmax;       /* number of sbas ephemeris */
    int ne,nemax;       /* number of precise ephemeris */
    int nc,ncmax;       /* number of precise clock */
    int na,namax;       /* number of almanac data */
    int nt,ntmax;       /* number of tec grid data */
    eph_t *eph;         /* GPS/QZS/GAL/BDS/IRN ephemeris */
    geph_t *geph;       /* GLONASS ephemeris */
    seph_t *seph;       /* SBAS ephemeris */
    double utc_gps[8];  /* GPS delta-UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
    double utc_glo[8];  /* GLONASS UTC time parameters {tau_C,tau_GPS} */
    double utc_gal[8];  /* Galileo UTC parameters */
    double utc_qzs[8];  /* QZS UTC parameters */
    double utc_cmp[8];  /* BeiDou UTC parameters */
    double utc_irn[9];  /* IRNSS UTC parameters {A0,A1,Tot,...,dt_LSF,A2} */
    double utc_sbs[4];  /* SBAS UTC parameters */
    double ion_gps[8];  /* GPS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_gal[4];  /* Galileo iono model parameters {ai0,ai1,ai2,0} */
    double ion_qzs[8];  /* QZSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    double ion_irn[8];  /* IRNSS iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
    int glo_fcn[32];    /* GLONASS FCN + 8 */
    double cbias[MAXSAT][3];    /* satellite DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
    double rbias[MAXRCV][2][3]; /* receiver DCB (0:P1-P2,1:P1-C1,2:P2-C2) (m) */
    double ISC[MAXSAT][MAXCODE];/* satellite ISC */
    sbssat_t sbssat;            /* SBAS satellite corrections */
    sbsion_t sbsion[MAXBAND+1]; /* SBAS ionosphere corrections */
    pcv_t pcvs[MAXSAT];  /* satellite antenna pcv */
    dgps_t dgps[MAXSAT]; /* DGPS corrections */
    ssr_t ssr[MAXSAT];   /* SSR corrections */
    erp_t erp;           /* earth rotation parameters */
} nav_t;

typedef struct sta {      /* station parameter type */
    char name   [MAXANT]; /* marker name */
    char marker [MAXANT]; /* marker number */
    char antdes [MAXANT]; /* antenna descriptor */
    char antsno [MAXANT]; /* antenna serial number */
    char rectype[MAXANT]; /* receiver type descriptor */
    char recver [MAXANT]; /* receiver firmware version */
    char recsno [MAXANT]; /* receiver serial number */
    int antsetup;       /* antenna setup id */
    int itrf;           /* ITRF realization year */
    int deltype;        /* antenna delta type (0:enu,1:xyz) */
    double pos[3];      /* station position (ecef) (m) */
    double del[3];      /* antenna position delta (e/n/u or x/y/z) (m) */
    double hgt;         /* antenna height (m) */
    int glo_cp_align;   /* GLONASS code-phase alignment (0:no,1:yes) */
    double glo_cp_bias[4]; /* GLONASS code-phase biases {1C,1P,2C,2P} (m) */
    gtime_t next_time;     /* new site observation time or first new site observation time (for base station) */
    struct sta* next_sta;  /* next station parameter (for base station) */
} sta_t;

typedef struct {        /* solution type */
    gtime_t time;       /* time (GPST) */
    double rr[6];       /* position/velocity (m|m/s) {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    double qr[6];       /* position variance/covariance (m^2) */
                        /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
                        /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    double qv[6];       /* velocity variance/covariance (m^2/s^2) */
    double qa[6];       /* attitude variance/covariance (rad^2)*/
    double dtr[6*NFREQ];/* receiver clock bias to time systems (s) */
    double ddrt;        /* receiver clock bias rate (m/s) */
    double dops[4];     /* DOP (dilution of precision) */
    double att[3];      /* attitude (rad) {roll,pitch,yaw} */

    uint8_t type;       /* type (0:xyz-ecef,1:enu-baseline) */
    uint8_t stat;       /* solution status (SOLQ_???) */
    uint8_t ns;         /* number of valid satellites */
    uint8_t nd;         /* number of valid satellites for doppler velocity solution */
    uint8_t motion;     /* receiver motion status (1:static,2:vehicle travels in straight line,3:low dynamic,4:middle dynamic,5:high dynamic)*/
    uint8_t inherit;    /* flag of inherit ambiguity */

    double age;         /* age of differential (s) */
    double ratio;       /* AR ratio factor for validation */
    double Ps;          /* Ambiguity Success */
    double ADOP;        /* Ambiguity dilution of precision */
    double CDOP;        /* carrier-phase dilution of precision */
    double thres;       /* AR ratio threshold for validation */
    double sigr[3][4];  /* ratio of float/fix solution DD residuals sigma/RMS/MAD/MEAN ([0]:all,[1]:phase,[2]:code) */
    double accl[3];     /* IMU accl measurement */
    double gyro[3];     /* IMU gyro measurement */
} sol_t;

typedef struct {        /* solution buffer type */
    int n,nmax;         /* number of solution/max number of buffer */
    int cyclic;         /* cyclic buffer flag */
    int start,end;      /* start/end index */
    gtime_t time;       /* current solution time */
    sol_t *data;        /* solution data */
    double rb[3];       /* reference position {x,y,z} (ecef) (m) */
    uint8_t buff[MAXSOLMSG+1]; /* message buffer */
    int nb;             /* number of byte in message buffer */
} solbuf_t;

typedef struct {        /* RTCM control struct type */
    int staid;          /* station id */
    int stah;           /* station health */
    int seqno;          /* sequence number for rtcm 2 or iods msm */
    int outtype;        /* output message type */
    gtime_t time;       /* message time */
    gtime_t time_s;     /* message start time */
    obs_t obs;          /* observation data (uncorrected) */
    nav_t nav;          /* satellite ephemerides */
    sta_t sta;          /* station parameters */
    dgps_t *dgps;       /* output of dgps corrections */
    ssr_t ssr[MAXSAT];  /* output of ssr corrections */
    char msg[128];      /* special message */
    char msgtype[256];  /* last message type */
    char msmtype[7][128]; /* msm signal types */
    int obsflag;        /* obs data complete flag (1:ok,0:not complete) */
    int ephsat;         /* input ephemeris satellite number */
    int ephset;         /* input ephemeris set (0-1) */
    double cp[MAXSAT][NFREQ+NEXOBS];     /* carrier-phase measurement */
    uint16_t lock[MAXSAT][NFREQ+NEXOBS]; /* lock time */
    uint16_t loss[MAXSAT][NFREQ+NEXOBS]; /* loss of lock count */
    gtime_t lltime[MAXSAT][NFREQ+NEXOBS];/* last lock time */
    int nbyte;          /* number of bytes in message buffer */ 
    int nbit;           /* number of bits in word buffer */ 
    int len;            /* message length (bytes) */
    uint8_t buff[1200]; /* message buffer */
    uint32_t word;      /* word buffer for rtcm 2 */
    uint32_t nmsg2[100];/* message count of RTCM 2 (1-99:1-99,0:other) */
    uint32_t nmsg3[400];/* message count of RTCM 3 (1-299:1001-1299,300-329:4070-4099,0:other) */
    char opt[256];      /* RTCM dependent options */
} rtcm_t;

typedef struct {        /* SNR mask type */
    int ena[2];         /* enable flag {rover,base} */
    double mask[NFREQ][9]; /* mask (dBHz) at 5,10,...85 deg */
} snrmask_t;

typedef struct {        /* processing options type */
    int mode;           /* positioning mode (PMODE_???) */
    int ar;             /* AR mode (0:off,1:on) */
    int navsys;         /* navigation system */
    snrmask_t snrmask;  /* SNR mask */
    double elmincode;   /* elevation mask angle (rad) for pseudorange */
    double elmindopp;   /* elevation mask angle (rad) for doppler */
    double elmin;       /* elevation mask angle (rad) for RTD/RTK */
    int sateph;         /* satellite ephemeris/clock (EPHOPT_???) */
    int minlock;        /* min lock count to fix ambiguity */
    int minfix;         /* min fix count to hold ambiguity */
    int ionoopt;        /* ionosphere option (IONOOPT_???) */
    int tropopt;        /* troposphere option (TROPOPT_???) */
    int dynamics;       /* dynamics model (0:off,1:on) */
    double maxgdop;     /* reject threshold of gdop */
    double thresar[4];  /* AR validation threshold (ratio,Ps,ADOP,CDOP) */
    double elmaskar;    /* elevation mask of AR for rising satellite (rad) */
    double elmaskinherit;
                        /* elevation mask to inherit ambiguity (deg) */
    double thresslip;   /* slip threshold of geometry-free phase (m) */
    double maxtdiff;    /* max difference of time (sec) */
    double psd[7];      /* [0-2]: acceleration PSD per axis (m^2/s^3) (H/V)
                         * [3]  : receiver clock phase-drift PSD (m^2/s)
                         * [4]  : receiver clock frequency-drift PSD (m^2/s^3)
                         * [5]  : trop random walk PSD
                         * [6]  : ambiguity bias process-noise PSD
                         * */
    double min_cn0_code;     /* CN0 mask (dBHz) for pseudorange measurement */
    double min_cn0_dopp;     /* CN0 mask (dBHz) for doppler measurement */
    double min_cn0_phas;     /* CN0 mask (dBHz) for phase measurement */
    char anttype[2][MAXANT]; /* antenna types {rover,base} */
    double antdel[2][3];     /* antenna delta {{rov_e,rov_n,rov_u},{ref_e,ref_n,ref_u}} */
    pcv_t pcvr[2];           /* receiver antenna parameters {rov,base} */
    int gloar;               /* fix GLONASS DD ambiguity */
    double rb[3];            /* base station position (ECEF/m) */
    int udvft;               /* velocity filter (0:off,1:on) */
    int udppp;               /* PPP option (0:off,1:on) */
    int udsig;               /* update RTK solution by adjust sigma (available for GNSS process) */
    int udrmm;               /* update STD of DD residuals */
    int uddir;               /* GNSS process direction (0:forward, 1:backward) */
    int maxiter;             /* max iteration for RTK filter */
    int inherit;             /* inherit DD ambiguity options (0:off,1:on) */
    double odisp[2][6*11];   /* ocean tide loading parameters {rov,base} */
    double sclkstab;         /* satellite clock stability (sec/sec) */
    int uddynamics;
    int adjsig0;
    int modear;
    double elholdamb;        /* elevation mask of AR hold (deg) */
    int enablertkflt;
    int enablertkwl;
} prcopt_t;

typedef struct {        /* solution options type */
    int posf;           /* solution format (SOLF_???) */
    int times;          /* time system (TIMES_???) */
    int timef;          /* time format (0:sssss.s,1:yyyy/mm/dd hh:mm:ss.s) */
    int timeu;          /* time digits under decimal point */
    int degf;           /* latitude/longitude format (0:ddd.ddd,1:ddd mm ss) */
    int outhead;        /* output header (0:no,1:yes) */
    int outopt;         /* output processing options (0:no,1:yes) */
    int outvel;         /* output velocity options (0:no,1:yes) */
    int datum;          /* datum (0:WGS84,1:Tokyo) */
    int height;         /* height (0:ellipsoidal,1:geodetic) */
    int geoid;          /* geoid model (0:EGM96,1:JGD2000) */
    int solstatic;      /* solution of static mode (0:all,1:single) */
    int sstat;          /* solution statistics level (0:off,1:states,2:residuals) */
    int trace;          /* debug trace level (0:off,1-5:debug) */
    double nmeaintv[2]; /* nmea output interval (s) (<0:no,0:all) */
                        /* nmeaintv[0]:gprmc,gpgga,nmeaintv[1]:gpgsv */
    char sep[64];       /* field separator */
    char prog[64];      /* program name */
    double maxsolstd;   /* max std-dev for solution output (m) (0:all) */
} solopt_t;

typedef struct {        /* RINEX options type */
    gtime_t ts,te;      /* time start/end */
    double tint;        /* time interval (s) */
    double ttol;        /* time tolerance (s) */
    double tunit;       /* time unit for multiple-session (s) */
    double rnxver;      /* RINEX version */
    int navsys;         /* navigation system */
    int obstype;        /* observation type */
    int freqtype;       /* frequency type */
    char mask[7][64];   /* code mask {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    char staid [32];    /* station id for rinex file name */
    char prog  [32];    /* program */
    char runby [32];    /* run-by */
    char marker[64];    /* marker name */
    char markerno[32];  /* marker number */
    char markertype[32];/* marker type (ver.3) */
    char name[2][32];   /* observer/agency */
    char rec [3][32];   /* receiver #/type/vers */
    char ant [3][32];   /* antenna #/type */
    double apppos[3];   /* approx position x/y/z */
    double antdel[3];   /* antenna delta h/e/n */
    char comment[MAXCOMMENT][64]; /* comments */
    char rcvopt[256];   /* receiver dependent options */
    unsigned char exsats[MAXSAT]; /* excluded satellites */
    int scanobs;        /* scan obs types */
    int outiono;        /* output iono correction */
    int outtime;        /* output time system correction */
    int outleaps;       /* output leap seconds */
    int autopos;        /* auto approx position */
    int halfcyc;        /* half cycle correction */
    int sep_nav;        /* separated nav files */
    gtime_t tstart;     /* first obs time */
    gtime_t tend;       /* last obs time */
    gtime_t trtcm;      /* approx log start time for rtcm */
    char tobs[7][MAXOBSTYPE][4]; /* obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
    int nobs[7];        /* number of obs types {GPS,GLO,GAL,QZS,SBS,CMP,IRN} */
} rnxopt_t;

typedef struct {           /* satellite status type */
    double azel[2];        /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ];    /* residuals of pseudorange (m) */
    double resd[NFREQ];    /* residuals of doppler (m) */
    double resc[NFREQ];    /* residuals of carrier-phase (m) */
    double respflt[NFREQ]; /* residuals of pseudorange (m) of RTK filter*/
    double rescflt[NFREQ]; /* residuals of carrier-phase (m) of RTK filter */
    double rs[6],dts[2];   /* satellite positions and clocks */
    double var_s;          /* satellite position variance for ephemeris*/
    uint8_t svh;           /* satellite health flag */
    uint8_t vs_spp_lsq[NFREQ]; /* valid satellite flag for pseudorange in SPP LSQ */
    uint8_t vd_spp_lsq[NFREQ]; /* valid satellite flag for doppler in SPP LSQ */
    uint8_t vs_spp_flt[NFREQ]; /* valid satellite flag for pseudorange in SPP filter */
    uint8_t vd_spp_flt[NFREQ]; /* valid satellite flag for doppler in SPP filter */
    uint8_t vs_rtd_lsq[NFREQ]; /* valid satellite flag for pseudorange in RTD LSQ */
    uint8_t vs_rtd_flt[NFREQ]; /* valid satellite flag for pseudorange in RTD filter */
    uint8_t vd_rtd_flt[NFREQ]; /* valid satellite flag for doppler in RTD filter */
    uint8_t vs_ppp_flt[NFREQ]; /* valid satellite flag for PPP filter */
    uint8_t vsat[NFREQ];   /* valid satellite flag for RTK */
    uint8_t vsatflt[NFREQ];   /* valid satellite flag for RTK filter */
    uint32_t SNR[NFREQ];   /* signal strength (0.001 dBHz) */
    uint8_t  fix[NFREQ];   /* ambiguity fix flag (0:init,1:float,2:fix,3:hold) */
    uint8_t ufix[NFREQ];   /* ambiguity fix flag (0:init,1:float,2:fix,3:hold) */
    int     slip[NFREQ];   /* cycle-slip flag (-1:unknown,0:no slip,1:slip) */
    int      dSNR[NFREQ];  /* increments of signal strength (0.001 dBHz) between current and precious epoch */
    uint32_t outc[NFREQ];  /* obs outage counter of phase */
    uint32_t outcflt[NFREQ];/* obs outage counter of phase of RTK filter */
    uint8_t half[NFREQ];   /* half-cycle valid flag */
    int     lock[NFREQ];   /* ambiguity lock flag */
    int     lockflt[NFREQ];/* ambiguity lock flag of RTK filter */
    uint32_t rejc[NFREQ];  /* reject flag (0:reject,1:accept) */
    uint8_t  pdcs[NFREQ];  /* phase/doppler consistency flag */
    uint8_t  cdcs[NFREQ];  /* code/doppler consistency flag */
    uint8_t  cdif[NFREQ];  /* consistent of rcv clock drift from inter-freq (based on L1,0:inconsistent,1:consistent) */
    double gf[NFREQ-1];    /* geometry-free phase (m) for RTK */
    double mw[NFREQ-1];    /* MW-LC (m) for RTK */
    double gf2[NFREQ-1];   /* geometry-free phase (m) for PPP */
    double mw2[NFREQ-1];   /* MW-LC (m) for PPP */
    double phw;            /* phase windup (cycle) */
    double BE[NFREQ];      /* L1/L2/L5 double-difference ambiguity from RTK-WL */
    obsd_t obs[2];         /* rover observation data of previous/current epoch */
    int ifreq[NFREQ][3];   /* indexs of RTK linear combination observation */
    int cfreq[NFREQ][2];   /* coefficient of WL linear combination observation */
    int nifreq;            /* number of RTK linear combination observation */
    int BEfixcnt[NFREQ];   /* fix counts of L1/L2/L5 double-difference ambiguity from RTK-WL */
    int BEfixflg[NFREQ];   /* fix flag of L1/L2/L5 double-difference ambiguity from RTK-WL */
    uint8_t refsat[NFREQ];
    uint8_t lflg[NFREQ];   /* large carrier-phase residual flag */
    double dd[NFREQ];      /* double-difference ambiguity */
    double ddm[NFREQ];     /* double-difference ambiguity DD-geometric-range/troposphere delay */
    double eh[3];
    double freq[NFREQ];
} ssat_t;

/* global variables ----------------------------------------------------------*/
extern const double chisqr[];             /* chi-sqr(n) table (alpha=0.001) */
extern const prcopt_t prcopt_default;     /* default positioning options */
extern const solopt_t solopt_default;     /* default solution output options */
extern const sbsigpband_t igpband1[9][8]; /* SBAS IGP band 0-8 */
extern const sbsigpband_t igpband2[2][5]; /* SBAS IGP band 9-10 */

/* satellites, systems, codes functions --------------------------------------*/
extern int  satno   (int sys, int prn);
extern int  satsys  (int sat, int *prn);
extern int  satid2no(const char *id);
extern void satno2id(int sat, char *id);
extern uint8_t obs2code(const char *obs, int *freq);
extern char *code2obs(uint8_t code, int *freq);
extern double code2freq(int sys, uint8_t code, int fcn);
extern double sat2freq(int sat, uint8_t code, const nav_t *nav);
extern int  code2idx(int sys, uint8_t code);
extern int  satexclude(int sat, double var, int svh, const prcopt_t *opt);
extern int  testsnr(int base, int freq, double el, double snr, const snrmask_t *mask);
extern void setcodepri(int sys, int idx, const char *pri);
extern int  getcodepri(int sys, uint8_t code, const char *opt);

/* matrix and vector functions -----------------------------------------------*/
extern double *mat  (int n, int m);
extern int    *imat (int n, int m);
extern double *zeros(int n, int m);
extern double *eye  (int n);
extern double dot (const double *a, const double *b, int n);
extern double norm(const double *a, int n);
extern double sum(const double *a, int n);
extern void cross3(const double *a, const double *b, double *c);
extern int  normv3(const double *a, double *b);
extern void matcpy(double *A, const double *B, int n, int m);
extern void matmul(const char *tr, int n, int k, int m, double alpha,
                   const double *A, const double *B, double beta, double *C);
extern void matmul33(const char *tr,const double *A,const double *B,const double *C,
                     int n,int p,int q,int m,double *D);
extern int  matinv(double *A, int n);
extern int mateigenvalue(const double* A, int n, double *u, double *v);
extern int matdet(const double *A, int n, double *det);
extern double mattrace(const double *A, int n);
extern int  solve (const char *tr, const double *A, const double *Y, int n, int m, double *X);
extern int  lsq   (const double *A, const double *y, int n, int m, double *x, double *Q);
extern int  wlsq  (const double *H, const double *R, const double *v, int nx, int nv, double *dx, double *Q);
extern int  filter(double *x, double *P, const double *H, const double *v,
                   const double *R, int n, int m, const int *ix_, int nx_);
extern int ifilter(double *x, double *P, const double *H, const double *v,
                   const double *R, int n, int m, const int *ix_, int nx, const double *xp);
extern int adrbfilter(double *x, double *P, const double *H, const double *v,
                      const double *R, int n, int m, const int *ix_, int nx_);
extern int smoother(const double *xf, const double *Qf, const double *xb,
                    const double *Qb, int n, double *xs, double *Qs);
extern void matfprint(const double *A, int n, int m, int p, int q, FILE *fp);
extern void vecadd(const double *a, const double *b, int n, double *c);
extern void vecsub(const double *a, const double *b, int n, double *c);

extern double median(const double *vec, int n);
extern double mad(const double *vec, int n);
extern double avg(const double *val, int n);
extern double stdv(const double *val, int n);
extern double norm_distri(double u);
extern double re_norm(double p);
extern double B(int n1, int n2, double x, double  *Ux);
extern double F(int n1, int n2, double x, double *f);
extern double t(int nn, double tv, double *f);
extern double re_t(int n,double p);
extern double re_F(int n1, int n2, double p);
extern void expm(const double *A, int n, double *Z);

/* time and string functions -------------------------------------------------*/
extern double  str2num(const char *s, int i, int n);
extern int     str2time(const char *s, int i, int n, gtime_t *t);
extern void    time2str(gtime_t t, char *str, int n);
extern gtime_t epoch2time(const double *ep);
extern void    time2epoch(gtime_t t, double *ep);
extern gtime_t gpst2time(int week, double sec);
extern double  time2gpst(gtime_t t, int *week);
extern gtime_t gst2time(int week, double sec);
extern double  time2gst(gtime_t t, int *week);
extern uint64_t time2u64(gtime_t t);
extern gtime_t bdt2time(int week, double sec);
extern double  time2bdt(gtime_t t, int *week);
extern char    *time_str(gtime_t t, int n);

extern gtime_t timeadd (gtime_t t, double sec);
extern double  timediff(gtime_t t1, gtime_t t2);
extern gtime_t gpst2utc(gtime_t t);
extern gtime_t utc2gpst(gtime_t t);
extern gtime_t gpst2bdt(gtime_t t);
extern gtime_t bdt2gpst(gtime_t t);
extern gtime_t timeget (void);
extern double  time2doy(gtime_t t);
extern double  utc2gmst(gtime_t t, double ut1_utc);

extern int adjgpsweek(int week);
extern uint32_t tickget(void);

extern int reppath(const char *path, char *rpath, gtime_t time, const char *rov, const char *base);
extern int reppaths(const char *path, char *rpaths[], int nmax, gtime_t ts, gtime_t te,
                    const char *rov, const char *base);

/* coordinates transformation ------------------------------------------------*/
extern void ecef2pos(const double *r, double *pos);
extern void pos2ecef(const double *pos, double *r);
extern void ecef2enu(const double *pos, const double *r, double *e);
extern void enu2ecef(const double *pos, const double *e, double *r);
extern void covenu  (const double *pos, const double *P, double *Q);
extern void covecef (const double *pos, const double *Q, double *P);
extern void xyz2enu (const double *pos, double *E);
extern void eci2ecef(gtime_t tutc, const double *erpv, double *U, double *gmst);
extern void deg2dms (double deg, double *dms, int ndec);
extern double dms2deg(const double *dms);

/* input and output functions ------------------------------------------------*/
extern int  sortobs(obs_t *obs);
extern void uniqnav(nav_t *nav);
extern int  screent(gtime_t time, gtime_t ts, gtime_t te, double tint);
extern void freeobs(obs_t *obs);
extern void freenav(nav_t *nav, int opt);
extern void freesta(sta_t *sta);

/* debug trace functions -----------------------------------------------------*/
extern void traceopen(const char *file);
extern void traceclose(void);
extern void tracelevel(int level);
extern void trace    (int level, const char *format, ...);
extern void tracet   (int level, const char *format, ...);
extern void tracemat (int level, const double *A, int n, int m, int p, int q);
extern void traceobs (int level, const obsd_t *obs, int n);
extern void tracenav (int level, const nav_t *nav);
extern void tracegnav(int level, const nav_t *nav);
extern void traceb   (int level, const uint8_t *p, int n);

/* platform dependent functions ----------------------------------------------*/
extern int execcmd(const char *cmd);
extern int expath(const char *path, char *paths[], int nmax);

/* positioning models --------------------------------------------------------*/
extern double satazel(const double *pos, const double *e, double *azel);
extern double geodist(const double *rs, const double *rr, double *e);
extern void dops(int ns, const double *azel, double elmin, double *dop);

/* atmosphere models ---------------------------------------------------------*/
extern double ionmodel(gtime_t t, const double *ion, const double *pos, const double *azel);
extern double ionmapf(const double *pos, const double *azel);
extern double ionppp(const double *pos, const double *azel, double re, double hion, double *pppos);
extern double tropmodel(gtime_t time, const double *pos, const double *azel, double humi);
extern double tropmapf(gtime_t time, const double *pos, const double *azel, double *mapfw);
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var);
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var);

/* antenna models ------------------------------------------------------------*/
extern void antmodel(const pcv_t *pcv, const double *del, const double *azel, int opt, double *dant);
extern void antmodel_s(const pcv_t *pcv, double nadir, double *dant);

/* earth tide models ---------------------------------------------------------*/
extern void sunmoonpos(gtime_t tutc, const double *erpv, double *rsun, double *rmoon, double *gmst);
extern void tidedisp(gtime_t tutc, const double *rr, int opt, const erp_t *erp,
                     const double *odisp, double *dr);

/* rinex functions -----------------------------------------------------------*/
extern int readrnx (const char *file, int rcv, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta);
extern int readrnxt(const char *file, int rcv, gtime_t ts, gtime_t te,
                    double tint, const char *opt, obs_t *obs, nav_t *nav, sta_t *sta);
extern int outrnxobsh(FILE *fp, const rnxopt_t *opt, const nav_t *nav);
extern int outrnxobsb(FILE *fp, const rnxopt_t *opt, const obsd_t *obs, int n,
                      int epflag);
extern int rtk_uncompress(const char *file, char *uncfile);
extern int readpcv(const char *file, pcvs_t *pcvs);
extern pcv_t *searchpcv(int sat, const char *type, gtime_t time, const pcvs_t *pcvs);

/* ephemeris and clock functions ---------------------------------------------*/
extern double eph2clk (gtime_t time, const eph_t  *eph);
extern double geph2clk(gtime_t time, const geph_t *geph);
extern double seph2clk(gtime_t time, const seph_t *seph);
extern void eph2pos (gtime_t time, const eph_t  *eph,  double *rs, double *dts, double *var);
extern void geph2pos(gtime_t time, const geph_t *geph, double *rs, double *dts, double *var);
extern void seph2pos(gtime_t time, const seph_t *seph, double *rs, double *dts, double *var);
extern int  satpos(gtime_t time, gtime_t teph, int sat, int ephopt,
                   const nav_t *nav, double *rs, double *dts, double *var, int *svh);
extern void satposs(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                    int sateph, double *rs, double *dts, double *var, int *svh);
extern void setseleph(int sys, int sel);
extern int  getseleph(int sys);
extern double baseline(const double *ru, const double *rb, double *dr);

/* receiver raw data functions -----------------------------------------------*/
extern uint32_t getbitu(const uint8_t *buff, int pos, int len);
extern int32_t  getbits(const uint8_t *buff, int pos, int len);
extern void setbitu(uint8_t *buff, int pos, int len, uint32_t data);
extern void setbits(uint8_t *buff, int pos, int len, int32_t  data);
extern uint32_t rtk_crc32 (const uint8_t *buff, int len);
extern uint32_t rtk_crc24q(const uint8_t *buff, int len);
extern uint16_t rtk_crc16 (const uint8_t *buff, int len);
extern int decode_word(uint32_t word, uint8_t *data);

/* solution functions --------------------------------------------------------*/
extern void initsolbuf(solbuf_t *solbuf, int cyclic, int nmax);
extern void freesolbuf(solbuf_t *solbuf);
extern int addsol(solbuf_t *solbuf, const sol_t *sol);
extern int readsol(char *files[], int nfile, solbuf_t *sol);
extern int readsolt(char *files[], int nfile, gtime_t ts, gtime_t te, double tint, int qflag, solbuf_t *sol);
extern int inputsol(uint8_t data, gtime_t ts, gtime_t te, double tint, int qflag, const solopt_t *opt, solbuf_t *solbuf);
extern int outsols(uint8_t *buff, const sol_t *sol, const double *rb, const solopt_t *opt);
extern int outsolexs(uint8_t *buff, const sol_t *sol, const ssat_t *ssat, const solopt_t *opt);
extern void outsol(FILE *fp, const sol_t *sol, const double *rb, const solopt_t *opt);
extern int outnmea_rmc(uint8_t *buff, const sol_t *sol);
extern int outnmea_gga(uint8_t *buff, const sol_t *sol);
extern int outnmea_gsa(uint8_t *buff, const sol_t *sol, const ssat_t *ssat);
extern int outnmea_gsv(uint8_t *buff, const sol_t *sol, const ssat_t *ssat);

/* GNSS solutions convert to KML----------------------------------------------*/
extern int sol2kml(const char *file, const solbuf_t *solbuf, int tcolor, int pcolor,
                   int outalt, int outtime);

/* observation error model----------------------------------------------------*/
extern double doppvarerr(const prcopt_t *opt, double cn0, double el, int sat);
extern double codevarerr(const prcopt_t *opt, double cn0, double el, int sat);
extern double phasevarerr(const prcopt_t *opt, double cn0, double el, int sat);


#ifdef __cplusplus
}
#endif
#endif
