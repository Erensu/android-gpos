/*------------------------------------------------------------------------------
 * agpos.h  : Android GNSS constants, types and function prototypes
 *
 * author   : sujinglan
 *-----------------------------------------------------------------------------*/
#ifndef AGPOS_H
#define AGPOS_H

#include "rtklib.h"
#include "android.h"

#ifdef __cplusplus
extern "C" {
#endif

/* type definitions ----------------------------------------------------------*/
typedef struct {                 /* SPP least square estimation solution type */
    gtime_t time;                /* time (GPST) */
    double rr[6];                /* position/velocity (m|m/s) {x,y,z,vx,vy,vz} or {e,n,u,ve,vn,vu} */
    double qr[6];                /* position variance/covariance (m^2) */
    double qv[6];                /* velocity variance/covariance (m^2/s^2) */
    double qt[5];                /* receiver clock bias to time systems variance (m)/clock drift variance (m/s) */
    double dtr[4*NFREQ];         /* receiver clock bias to time systems (s) */
    double ddrt;                 /* receiver clock bias rate (m/s) */

    double dop_pos[6];           /* DOP (dilution of precision) for position */
    double dop_vel[6];           /* DOP (dilution of precision) for velocity */

    double sig0_p0,sig0_v0;      /* position/velocity sig0=v'*P*v/(nv-nx) */
    double sig0_p1,sig0_v1;      /* position/velocity sig0=v'*v/(nv-nx) */
    double min_cn0[2];           /* cut-off CN0 of current LSQ solution (0: code, 1: doppler) */
    double dt;                   /* correction of priori observation coarse-time  */

    unsigned char psta;          /* position solution status */
    unsigned char vsta;          /* doppler velocity solution status */
    unsigned char pdeg;          /* degraded position solution status */
    unsigned char vdeg;          /* degraded position solution status */

    unsigned char ns;            /* number of valid satellites */
    unsigned char nd;            /* number of valid satellites for doppler velocity solution */
} aslsq_t;

typedef struct {                 /* SPP EKF filter solution type */
    gtime_t time;                /* last time (GPST) of EKF filter */
    int nx;                      /* number of estimated states */
    unsigned char stat;          /* EKF filter status */
    unsigned char ns;            /* number of valid satellites */
    unsigned char nd;            /* number of valid satellites for doppler velocity solution */
    double tt;                   /* time difference between current and previous (s) */
    double *xp,*Pp;              /* predict estimated states and their covariance */
    double  *x, *P;              /* measurement updated estimated states and their covariance */
} asflt_t;

typedef struct {                 /* satellite status type */
    unsigned char sys;           /* navigation system */
    double azel[2];              /* azimuth/elevation angles {az,el} (rad) */
    double resp[NFREQ];          /* residuals of pseudorange (m) */
    double resd[NFREQ];          /* residuals of doppler (m) */
    double rs[9],dts[3];         /* satellite positions and clocks */
    double var_s;                /* satellite position variance for ephemeris*/
    unsigned char vs_lsq[NFREQ]; /* valid satellite flag for pseudorange in wlsq */
    unsigned char vd_lsq[NFREQ]; /* valid satellite flag for doppler in wlsq */
    unsigned char vs_flt[NFREQ]; /* valid satellite flag for pseudorange in filter */
    unsigned char vd_flt[NFREQ]; /* valid satellite flag for doppler in filter */
    unsigned char SNR[NFREQ];    /* signal strength (0.25 dBHz) */
    unsigned int svh;            /* satellite health */
    gtime_t time[NFREQ];         /* current epoch time */
} assat_t;

typedef struct {                 /* Android GNSS single point position type */
    aslsq_t lsq;                 /* LSQ position for Android GNSS solution */
    asflt_t flt;                 /* EKF filter for Android GNSS solution */
} aspp_t;

typedef struct {                 /* Android GNSS position type */
    sol_t  sol;                  /* GNSS positioning solution for output */
    sol_t  sol_chip;             /* solution from gps chip */
    sol_t  sol_nlp;              /* solution from NLP */
    aspp_t spp;                  /* Android GNSS single positioning */
    double rb[6];                /* base position/velocity (ecef) (m|m/s) */
    assat_t ssat[MAXSAT];        /* satellite status */
    prcopt_t opt;                /* processing options */
} agpos_t;

/* function definitions-------------------------------------------------------*/
extern void agpinit(agpos_t *agp, const prcopt_t *opt);
extern void agpfree(agpos_t *agp);
extern int agpos(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav);

#ifdef __cplusplus
}
#endif
#endif