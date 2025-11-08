/* ----------------------------------------------------------------------------
 * errmdl.c: extend observation error model
 *
 * reference :
 *     [1] Groves P D , Jiang Z . Height Aiding, C/N0 Weighting and Consistency
 *         Checking for GNSS NLOS and Multipath Mitigation in Urban Areas.
 *     [2] RTCA (2006). Minimum Operational Performance Standards for Global
 *         Positioning System/Wide Area Augmentation System Airborne Equipment,
 *         DO-229D.
 *
 * author  : sujinglan
 * version : $Revision:$ $Date:$
 * history : 2018/08/31 1.0    new
 * ---------------------------------------------------------------------------*/
#include "rtklib.h"

#define SNR_EL_MDL         1

/* doppler measurement error variance ------------------------------------
 * args:    prcopt_t *opt  I  process options
 *          double cn0     I  C/N0 of observation data
 *          int sat        I  satellite no.
 *          double el      I  elevation of satellite
 * return: square of measurement error
 * ---------------------------------------------------------------------------*/
extern double doppvarerr(const prcopt_t *opt, double cn0, double el, int sat)
{
    const double cn0_zth=45.0;
    double var,sig=0.05;

    var=SQR(sig)*pow(10.0,(cn0_zth-cn0)/10.0);
    return var/SQR(sin(el));
}
/* pseudorange measurement error variance ------------------------------------
 * args:    prcopt_t *opt  I  process options
 *          double cn0     I  C/N0 of observation data
 *          double el      I  elevation of satellite
 *          int sat        I  satellite no.
 * return: square of measurement error
 * ---------------------------------------------------------------------------*/
extern double codevarerr(const prcopt_t *opt, double cn0, double el, int sat)
{
    const double cn0_zth=45.0;
    double var,sig=1.0;

    var=SQR(sig)*pow(10.0,(cn0_zth-cn0)/10.0);
    return var;
}
/* carrier phase measurements error variance ---------------------------------
 * args:    prcopt_t *opt  I  process options
 *          double cn0     I  C/N0 of observation data
 *          double el      I  elevation of satellite
 *          int sat        I  satellite no.
 * return: square of measurement error
 * ---------------------------------------------------------------------------*/
extern double phasevarerr(const prcopt_t *opt, double cn0, double el, int sat)
{
#if SNR_EL_MDL
    double cn0_zth=45.0,sig=0.003;
    return SQR(sig)*pow(10.0,(cn0_zth-cn0)/10.0)/SQR(sin(el));
#else
    double cn0_zth=45.0,cn0_thres=30.0;
    double sig0=0.001,sig1=0.003;

    if (cn0<cn0_thres) return SQR(sig1)*pow(10.0,(cn0_zth-cn0)/10.0)/SQR(sin(el));
    return SQR(sig0)+SQR(sig0)/SQR(sin(el));
#endif
}