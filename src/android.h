/*------------------------------------------------------------------------------
 * android.h : android gnss constants, types and function prototypes
 *
 * author   : sujinglan
 *-----------------------------------------------------------------------------*/
#ifndef ANDROID_H
#define ANDROID_H

#include "rtklib.h"

#ifdef __cplusplus
extern "C" {
#endif

/* constants -----------------------------------------------------------------*/
#define OBSQ_RAW      (1<<0)          /* observation type: raw */
#define OBSQ_PRCS_DP  (1<<1)          /* observation type: pseudorange smooth using doppler */
#define OBSQ_PRPD_DP  (1<<2)          /* observation type: pseudorange predict using precious epoch doppler */
#define OBSQ_DPPD_DP  (1<<3)          /* observation type: doppler using precious epoch doppler */
#define OBSQ_PRPD_FIT (1<<4)          /* observation type: predict pseudorange using polynomial fitting  */
#define OBSQ_PR_BAD   (1<<5)          /* observation type: bad */
#define OBSQ_PR_FIXAMB (1<<6)         /* observation type: need to fix 1ms ambiguity */
#define OBSQ_MULTP_UNKNOWN   (1<<7)   /* observation type: not available or the presence or absence of multipath is unknown */
#define OBSQ_MULTP_DETECT    (1<<8)   /* observation type: the measurement shows signs of multi-path */
#define OBSQ_MULTP_NOTDETECT (1<<9)   /* observation type: the measurement shows no signs of multi-path */
#define OBSQ_GLO_FCN         (1<<10)  /* observation type: glonass frequency channel number(FCN). since two
                                       * GLONASS satellites transmit navigation signals on the same carrier frequency,
                                       * assisted data is needed to identify the correct OSN.
                                       * */
#define OBSQ_PR_FIXAMB_FRAC  (1<<11)  /* observation type: fix pseudorange using fractional part */
#define OBSQ_PR_FIXAMB_RES   (1<<12)  /* observation type: fix 1ms ambiguity using pseudorange residual */
#define OBSQ_COARSE_TIME     (1<<13)  /* observation type: coarse observation time */

#define OBST_CODE      0x01
#define OBST_POS       0x02
#define OBST_ISB       0x04

/* type definitions ----------------------------------------------------------*/
    typedef struct {                  /* other information from android raw data type */
    double time[NFREQ];               /* time uncertainty (m) */
    double bias[NFREQ];               /* bias uncertainty (m) */
    double drift[NFREQ];              /* drift uncertainty per second (m) */
    double rt[NFREQ];                 /* received Sv time uncertainty (m) */
    double pr[NFREQ];                 /* pseudorange rate uncertainty per second (m/s) */
    double adr[NFREQ];                /* accumulated delta range uncertainty (m) */
    double cp[NFREQ];                 /* carrier phase uncertainty */
    long long tt[NFREQ];              /* signal transmitted time (ns) */
    double frq[NFREQ];                /* frequency (Hz) */
    double snr[NFREQ];                /* Signal-to-Noise ratio (SNR) in dB (post-correlation & integration) */
    double ambf[NFREQ];
    unsigned int stat[NFREQ];         /* satellite tracking state */
    unsigned int adr_stat[NFREQ];     /* status of accumulated delta range */
    unsigned int svid;
} aunc_t;

typedef struct {                      /* Android observation data record */
    gtime_t time;                     /* receiver sampling time (GPST) */
    aunc_t unc;                       /* observation uncertainty from android raw gnss data */
    unsigned char sat,rcv;            /* satellite/receiver number */
    unsigned int stat[NFREQ+NEXOBS];  /* observation status */
    unsigned char LLI [NFREQ+NEXOBS]; /* loss of lock indicator */
    unsigned char code[NFREQ+NEXOBS]; /* code indicator (CODE_???) */
    unsigned char type;               /* observation type (OBST_???) */
    unsigned char SNR[NFREQ+NEXOBS];  /* signal strength (0.25 dBHz) */
    double Pr [NFREQ+NEXOBS];         /* observation data pseudorange rate (m/s) */
    double L  [NFREQ+NEXOBS];         /* observation data carrier-phase (cycle) */
    double P  [NFREQ+NEXOBS];         /* observation data pseudorange (m) */
    double CD [NFREQ+NEXOBS];         /* receiver clock drift per second (m/s) */
    double ADR[NFREQ+NEXOBS];         /* accumulated delta range */
    double   D[NFREQ+NEXOBS];         /* observation data doppler frequency (Hz) */
    double AGC[NFREQ+NEXOBS];         /* AGC acts as a variable gain amplifier adjusting the power of the incoming signal */
    int hwclkdc[NFREQ+NEXOBS];        /* hardware clock discontinuity count */
    long long TimeNanos,FullBiasNanos;
} aobsd_t;

typedef struct {                      /* android GNSS observation data */
    int n,nmax;                       /* number of observation data/allocated */
    aobsd_t *data;                    /* observation data records */
} aobs_t;

typedef struct {                      /* receiver raw data control type */
    gtime_t time;                     /* message time */
    aobs_t obs;                       /* Android GNSS observation data */
    aobs_t obuf;                      /* Android GNSS observation data buffer */
    nav_t nav;                        /* satellite ephemerides */
    sol_t sol;                        /* Android GNSS solutions data */
    int ephsat;                       /* sat number of update ephemeris (0:no satellite) */
    int week;                         /* GPS week */
    unsigned char obsflag;            /* 1: update observation data,0: no update observation data */
    long long fullbias0;              /* fullbias value at first epoch for Android raw data decode
                                       * compute rxtime using fullbias so that includes rx clock drift since the first epoch
                                       * here as reference fullbias
                                       * */
    long long bias0;                  /* clock's sub-nanosecond bias as reference */
    long long timenanos;              /* GNSS receiver internal hardware clock value in nanoseconds. */
    int hwclk0;                       /* count of hardware clock discontinuities */
    int nbyte;                        /* number of bytes in message buffer */
    unsigned char buff[4096];         /* message buffer */
    double ptime;
} araw_t;

/* function definitions-------------------------------------------------------*/
extern int input_android(araw_t *raw, unsigned char data);
extern int input_androidf(araw_t  *raw, FILE *fp);

extern int fixglofcn(const nav_t *nav, aobsd_t *obs, int n, const double *rr, int f,
                     const prcopt_t *opt);
extern int fixobs(gtime_t teph, const nav_t *nav, aobsd_t *obs, int n, const double *rr,
                  int f, const prcopt_t *opt);
extern void udepos(const double *rr);
extern void udsatrcvdis(double pr, int sat);
extern void aobs2obs(const aobsd_t *aobs, obsd_t *obs);
extern void obs2aobs(const obsd_t *obs, aobsd_t *aobs);

extern void initaraw(araw_t *raw);
extern void freearaw(araw_t *raw);

#ifdef __cplusplus
}
#endif
#endif