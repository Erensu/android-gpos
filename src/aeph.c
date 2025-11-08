/*------------------------------------------------------------------------------
 * ephemeris.c : satellite ephemeris and clock functions
 *
 * author  : sujinglan
 * version : $Revision:$ $Date:$
 * history : 2019/06/16 1.1  new
*-----------------------------------------------------------------------------*/
#include "agpos.h"

/* satellite clock with broadcast ephemeris ----------------------------------*/
extern int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  double *dts);

/* satellite positions and clocks for Android GNSS positioning------------------*/
extern void satpossa(gtime_t teph, const aobsd_t *obs, int n, const nav_t *nav,
                     int ephopt, double *rs, double *dts, double *var, int *svh)
{
    static double STD_BRDCCLK=30.0;
    static aobsd_t obsp[MAXSAT]={0};
    gtime_t time[2*MAXOBS]={{0}};
    char prn[8];
    double dt,pr;
    int i,j;

    trace(3,"satpossa: teph=%s n=%d ephopt=%d\n",time_str(teph,3),n,ephopt);

    for (i=0;i<n&&i<2*MAXOBS;i++) {
        for (j=0;j<9;j++) rs [j+i*9]=0.0;
        for (j=0;j<3;j++) dts[j+i*3]=0.0;
        var[i]=0.0; svh[i]=0;

        satno2id(obs[i].sat,prn);

        /* search any pseudorange */
        for (j=0,pr=0.0;j<NFREQ;j++) if ((pr=obs[i].P[j])!=0.0) break;

        if (j>=NFREQ) {
            trace(2,"no pseudorange %s sat=%3s\n",time_str(obs[i].time,3),prn);
            if ((dt=timediff(obs[i].time,obsp[obs[i].sat-1].time))>5.0) {
                continue;
            }
            for (j=0;j<NFREQ;j++) {
                if (obsp[obs[i].sat-1].P[j]&&obsp[obs[i].sat-1].Pr[j]) {
                    pr=obsp[obs[i].sat-1].P[j]+obsp[obs[i].sat-1].Pr[j]*dt;
                    trace(2,"pseudorange through PR-rate predict ok sat=%3s dt=%6.3lf\n",prn,dt);
                    break;
                }
            }
            if (j>=NFREQ) continue;
        }
        else {
            /* save observation of precious epoch */
            obsp[obs[i].sat-1]=obs[i];
        }
        /* transmission time by satellite clock */
        time[i]=timeadd(obs[i].time,-pr/CLIGHT);

        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time[i],teph,obs[i].sat,nav,&dt)) {
            trace(3,"no broadcast clock %s sat=%3s\n",time_str(time[i],3),prn);
            continue;
        }
        time[i]=timeadd(time[i],-dt);

        /* satellite position and clock at transmission time */
        if (!satpos(time[i],teph,obs[i].sat,ephopt,nav,rs+i*9,dts+i*3,var+i,svh+i)) {
            trace(3,"no ephemeris %s sat=%3s\n",time_str(time[i],3),prn);
            continue;
        }
        /* if no precise clock available, use broadcast clock instead */
        if (dts[i*3]==0.0) {
            if (!ephclk(time[i],teph,obs[i].sat,nav,dts+i*3)) continue;
            dts[1+i*3]=0.0;
            *var=SQR(STD_BRDCCLK);
        }
    }
    for (i=0;i<n&&i<2*MAXOBS;i++) {
        trace(4,"%s sat=%3s rs=%13.3f %13.3f %13.3f vs=%13.3f %13.3f %13.3f as=%13.3f %13.3f %13.3f "
                "dts=%12.3f %12.3f %12.3f var=%7.3f svh=%02X\n",
              time_str(time[i],6),prn,rs[0+i*9],rs[1+i*9],rs[2+i*9],
              rs[3+i*9],rs[4+i*9],rs[5+i*9],
              rs[6+i*9],rs[7+i*9],rs[8+i*9],
              dts[i*3]*1E9,dts[i*3+1]*1E9,dts[i*3+2]*1E9,var[i],svh[i]);
    }
}