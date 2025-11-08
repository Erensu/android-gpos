/* ----------------------------------------------------------------------------
 * fixobsamb.c: solving for millisecond integer ambiguity
 *
 * references :
 *     [1] Digglen F V. A-GPS: Assisted GPS, GNSS, and SBAS[J]. Navigational
 *         Satellites, 2009.
 *
 * author   : sujinglan
 * version    : $Revision:$ $Date:$
 * history    : 2022/10/05 1.0    new
 * ---------------------------------------------------------------------------*/
#include "android.h"

/* constants------------------------------------------------------------------*/
#define SQR(x)                ((x)*(x))
#define MAX_VAR_EPH           SQR(300.0) /* max variance eph to reject satellite (m^2) */
#define STD_BRDCCLK           30.0       /* error of broadcast clock (m) */
#define CODE_LOCK             (1<<0)     /* tracking state: has code lock */
#define BIT_SYNC              (1<<1)     /* tracking state: has bit sync */
#define SUBFRAME_SYNC         (1<<2)     /* tracking state: has sub-frame sync. */
#define TOW_DECODED           (1<<3)     /* tracking state: has time-of-week decoded */
#define GLO_STRING_SYNC       (1<<6)     /* tracking state: has string sync */
#define GLO_TOD_DECODED       (1<<7)     /* tracking state: has time-of-day decoded */
#define TOW_KNOWN             (1<<14)    /* tracking state: has time-of-week known */
#define FIX_OBS_AMB           1          /* fix observation ambiguity */
#define FIX_OBS_GLO_FCN       1          /* fix glonass FCN ambiguity */
#define BDST_TO_GPST          14.0       /* leap seconds difference between bdst and gpst */

/* global variables-----------------------------------------------------------*/
static double epos[3]={                  /* apr. receiver position (ecef/m) */
        -2693832.0262,                   /* ecef-X */
        -4297606.6990,                   /* ecef-Y */
         3854094.4139                    /* ecef-Z */
};
static double satpr[MAXSAT]={0.0};       /* geometric distance bwt sat and rcv */
static unsigned int fcntbl[2*14]={       /* glonass FCN table: -7-+6 */
        14,10,26,26,27,27, 2, 6,         /* -7--4 */
        18,22, 9,13,12,16,11,15,         /* -3-+0 */
         1, 5,20,24,19,23,17,21,         /* +0-+4 */
         3, 7, 4, 8,                     /* +4-+6 */
};

/* satellite clock with broadcast ephemeris ----------------------------------*/
extern int ephclk(gtime_t time, gtime_t teph, int sat, const nav_t *nav,
                  double *dts);

/* update apr. receiver position----------------------------------------------*/
extern void udepos(const double *rr)
{
    matcpy(epos,rr,1,3);
}
/* update geometric distance bwt sat and rcv----------------------------------*/
extern void udsatrcvdis(double pr, int sat)
{
    if (sat<=0||sat>MAXSAT) return;
    satpr[sat-1]=pr;
}
/* satellite position/velocity/clock------------------------------------------*/
static void satpvt(gtime_t teph, const nav_t *nav, const aobsd_t *obs, int n, double *rs,
                   double *dts, int *svh, double *var, const prcopt_t *opt)
{
    gtime_t time;
    double dt=0.0;
    int i,j;

    trace(3,"satpvt: teph=%s n=%d ephopt=%d\n",time_str(teph,3),n,opt->sateph);

    for (i=0;i<n&&i<MAXOBS;i++) {
        var[i]=svh[i]=0;
        for (j=0;j<6;j++) rs[9*i+j]=0.0;
        for (j=0;j<3;j++) dts[3*i+j]=0.0;

        /* transmission time by satellite clock */
        time=timeadd(obs[i].time,-satpr[obs[i].sat-1]/CLIGHT);

        /* satellite clock bias by broadcast ephemeris */
        if (!ephclk(time,teph,obs[i].sat,nav,&dt)) {
            continue;
        }
        time=timeadd(time,-dt);

        /* satellite position and clock at transmission time */
        if (!satpos(time,teph,obs[i].sat,opt->sateph,nav,rs+i*9,dts+i*3,var+i,svh+i)) {
            continue;
        }
    }
    for (i=0;i<n&&i<2*MAXOBS;i++) {

        char prn[8];
        satno2id(obs[i].sat,prn);
        trace(4,"%s sat=%3s rs=%13.3f %13.3f %13.3f dts=%12.3f var=%7.3f svh=%02X\n",
                time_str(teph,6),prn,rs[i*9],rs[1+i*9],rs[2+i*9],
                dts[i*3]*1E9,var[i],svh[i]);
    }
}
/* search reference satellite--------------------------------------------------*/
static int selrefsat(const aobsd_t *obs, const int *tag, int n, int rf, int flag)
{
    int i,j;

    trace(3,"selrefsat: rf=%d flag=%d\n",rf,flag);

    /* select reference satellite */
    for (j=-1,i=0;i<n;i++) {
        if (!obs[i].P[rf]) continue;

        /* satellite already have excluded */
        if (tag[i]) continue;

        /* disable galileo/glonass/bds */
        if (satsys(obs[i].sat,NULL)!=SYS_GPS) {
            continue;
        }
        /* based on lock and tow status */
        if (flag) {
            if (((obs[i].unc.stat[rf]&CODE_LOCK)&&(obs[i].unc.stat[rf]&TOW_DECODED))||
                 (obs[i].unc.stat[rf]&TOW_KNOWN)) {
                j=i; break;
            }
        }
        else {
            /* based on tow and valid pseudorange */
            if ((obs[i].unc.stat[rf]&TOW_DECODED)||(obs[i].unc.stat[rf]&TOW_KNOWN)) {
                j=i; break;
            }
        }
    }
    /* index of ref sat */
    return j;
}
/* check observation ambiguity based on prng/range----------------------------*/
static void chk_obsamb(const nav_t *nav, aobsd_t *obs, int n, int f, int j, int rf, const double *r0,
                       const double *dts)
{
    double N,ri,rj,dr;
    char prn[8];
    int i;

    trace(3,"chk_obsamb: n=%d f=%d j=%d\n",n,f,j);

    for (i=0;i<n&&j>=0;i++) {
        if (i==j||r0[i]<=0.0||!obs[i].P[f]) continue;

        ri=obs[i].P[ f]-(r0[i]-CLIGHT*dts[i*3]);
        rj=obs[j].P[rf]-(r0[j]-CLIGHT*dts[j*3]);
        dr=ri-rj;

        if (!(N=round(dr/CLIGHT*1E3))) continue;
        if (fabs(dr/CLIGHT*1E3-N)<1E-1) {
            obs[i].P[f]-=CLIGHT*N/1E3;
        }
        else {
            continue;
        }
        satno2id(obs[i].sat,prn);

        trace(3,"check obs amb: sat=%3d(%3s) f=%d dt=%8.3lfms N=%3.1lf dr=%10.4lf pr=%10.4lf\n",
              obs[i].sat,prn,f,(obs[i].P[f]/CLIGHT)*1E3,N,dr,obs[i].P[f]);
    }
}
/* fix observation ambiguity by other method-----------------------------------*/
static int fixambtt_(const nav_t *nav, aobsd_t *obs, int n, const double *rrp,
                     double *rs, double *dts, double *vare, int *svh,
                     double *r0, int *ref_sat, int *ref_frq, int f, const prcopt_t *opt)
{
    double *N,*z,e[3];
    double ambf,amb,P,dr,ri,rj;
    int i,j,k=-1,week,nfix,rf,frq[MAXSAT]={0};
    char prn[8],prn_ref[8];
    gtime_t gtt;

    trace(3,"fixambtt_: n=%d f=%d time=%8.4lf\n",n,f,time2gpst(obs[0].time,NULL));

    N=zeros(1,n);
    z=zeros(1,n);

    if (f==0) {
        /* satellite position */
        satpvt(obs[0].time,nav,obs,n,rs,dts,svh,vare,opt);

        /* geometric distance/reference sat */
        for (i=0;i<n;i++) {
            r0[i]=geodist(rs+9*i,rrp,e);
        }
    }
    /* get gps week */
    time2gpst(obs[0].time,&week);
    if (week<=0) time2gpst(timeget(),&week);

    /* Nms of all satellites */
    for (i=0;i<n;i++) {
        if (r0[i]<=0.0) continue;
        rf=f;

        if (!obs[i].unc.tt[f]) {
            for (rf=-1,j=0;j<NFREQ;j++) {
                if (j==f||!obs[i].unc.tt[j]) continue;
                rf=j; break;
            }
            if (rf==-1) continue;
        }
        if (!(obs[i].stat[rf]&OBSQ_PR_FIXAMB)) continue;

        N[i]=(int)(r0[i]/CLIGHT*1E3);
        gtt=gpst2time(week,obs[i].unc.tt[rf]*1E-9);
        z[i]=timediff(obs[i].time,gtt)*1E3-(int)(timediff(obs[i].time,gtt)*1E3);
        frq[i]=rf;
    }
    /* select reference sat */
    for (j=-1,i=0;i<n;i++) {
        if ((obs[i].stat[frq[i]]&OBSQ_GLO_FCN)||!obs[i].sat) continue;
        if (N[i]<=0.0) continue;
        if (j<0||obs[i].SNR[frq[i]]>obs[j].SNR[frq[j]]) {
            j=i;
        }
        else if (k<0||obs[i].SNR[frq[i]]>obs[k].SNR[frq[k]]) {
            k=i;
        }
    }
    if (j<0||k<0) {
        free(z); free(N);
        return 0;
    }
    if (ref_sat) *ref_sat=j;
    if (ref_frq) *ref_frq=frq[j];

    /* solve ambiguity of 1ms's integers */
    for (nfix=0,i=0;i<n;i++) {
        if (!(obs[i].stat[f]&OBSQ_PR_FIXAMB)) continue;
        if (i==j||N[i]<=0.0||z[i]==0.0) continue;
        ambf=N[j]+z[j]-z[i]+(r0[i]-r0[j])/CLIGHT*1E3+(dts[3*j]-dts[3*i])*1E3;

        /* fractional part */
        obs[i].unc.ambf[f]=ambf-round(ambf);

        /* check ambiguity of 1ms's integers */
        if (fabs(ambf-round(ambf))>1E-1) {
            satno2id(obs[i].sat,prn);
            trace(2,"ambf's value fail: sat=%s ambf=%.3lf\n",prn,ambf);
            continue;
        }
        amb=round(ambf);

        /* fix full pseudorange */
        P=(amb+z[i])*1E-3*CLIGHT;
        if (!obs[i].P[f]) {
            obs[i].P[f]=P; /* fix */
        }
        else if (obs[i].P[f]>1E7&&obs[i].P[f]<4E7) {
            if (fabs(obs[i].P[f]-P)>1E3) {
                trace(2,"fix amb fail: sat=%3d f=%d dP=%.3lf\n",obs[i].sat,f,obs[i].P[f]-P);
                continue;
            }
            obs[i].P[f]=P;
        }
        else {
            obs[i].P[f]=P;
        }
        /* fix uncertainty ambiguity */
        if (obs[i].unc.rt[f]>1E3) {
            obs[i].unc.rt[f]=10.0;
        }
        obs[i].stat[f]|=OBSQ_PR_FIXAMB_FRAC;

        satno2id(obs[j].sat,prn_ref);
        satno2id(obs[i].sat,prn);
        trace(3,"fix amb: sat=%3d(%3s) f=%d dt=%8.3lfms ambf=%15.4lf pr=%10.4lf ref=%3s\n",
              obs[i].sat,prn,f,amb+z[i],ambf,obs[i].P[f],prn_ref);
        nfix++;
    }
    /* solve ambiguity of reference satellite */
    if (k>=0&&j>=0&&!(obs[j].stat[frq[j]]&OBSQ_PR_FIXAMB_FRAC)) {
        ambf=N[k]+z[k]-z[j]+(r0[j]-r0[k])/CLIGHT*1E3+(dts[3*k]-dts[3*j])*1E3;
        amb=round(ambf);

        /* fix full pseudorange of ref sat */
        P=(amb+z[j])*1E-3*CLIGHT;
        obs[j].P[frq[j]]=P;

        obs[j].stat[frq[j]]|=OBSQ_PR_FIXAMB_FRAC;

        satno2id(obs[k].sat,prn_ref);
        satno2id(obs[j].sat,prn);

        trace(3,"fix amb: sat=%3d(%3s) f=%d dt=%8.3lfms ambf=%15.4lf pr=%10.4lf ref=%3s\n",
              obs[j].sat,prn,frq[j],amb+z[j],ambf,obs[j].P[frq[j]],prn_ref);
        nfix++;
    }
    /* check ambiguity of 1ms's integers through sd */
    for (i=0;i<n;i++) {
        if (!(obs[i].stat[f]&OBSQ_PR_FIXAMB_FRAC)) continue;
        if (i==j||!r0[i]) continue;

        ri=obs[i].P[f]-(r0[i]-CLIGHT*dts[i*3]);
        rj=obs[j].P[frq[j]]-(r0[j]-CLIGHT*dts[j*3]);
        dr=ri-rj;

        if (!(N[i]=round(dr/CLIGHT*1E3))) continue;
        if (fabs(dr/CLIGHT*1E3-N[i])<1E-1) {
            obs[i].P[f]-=CLIGHT*N[i]/1E3;
            obs[i].stat[f]|=OBSQ_PR_FIXAMB_RES;
        }
        else {
            satno2id(obs[i].sat,prn);
            trace(2,"check 1ms's ambiguity fail: sat=%s ambf=%.3lf\n",
                  prn,dr/CLIGHT*1E3);
            continue;
        }
        satno2id(obs[i].sat,prn);

        trace(3,"fix amb res: sat=%3d(%3s) f=%d dt=%8.3lfms N=%3.1lf dr=%10.4lf pr=%10.4lf\n",
              obs[i].sat,prn,f,(obs[i].P[f]/CLIGHT)*1E3,N[i],dr,obs[i].P[f]);
    }
    free(z); free(N);
    return nfix;
}
/* align f2/f3 1ms's ambiguity based on f1------------------------------------*/
static void fixf2f3amb(const nav_t *nav, aobsd_t *obs, int n, int f)
{
    double ambf,amb;
    int i;

    trace(3,"fixf2f3amb: n=%d f=%d\n",n,f);

    /* fix f2/f3 1ms ambiguity based on f1 */
    for (i=0;i<n;i++) {
        if (obs[i].P[f]<=0.0||obs[i].P[0]<=0.0) continue;

        ambf=(obs[i].P[f]-obs[i].P[0])/CLIGHT*1E3;
        if (!(amb=round(ambf))) continue;

        char prn[8];
        satno2id(obs[i].sat,prn);

        obs[i].P[f]-=CLIGHT*amb/1E3;
        trace(3,"fix f2/f3 amb: sat=%3s ambf=%.3lf\n",prn,ambf);
    }
}
/* recover and check sv signal transmit time/fix GLONASS FCN ambiguity ---------
 * args:  gtime_t teph  I   time to select ephemeris (gpst)
 *        obsd_t *obs   I  observation data
 *        int n         I  number of observation data
 *        nav_t *nav    I  navigation data
 *        double *rr    I  receriver position (ecef)
 *        int f         I  frequency no.
 *        prcopt_t *opt I  process options
 * return: numbers of fix ambiguity
 * ---------------------------------------------------------------------------*/
extern int fixobs(gtime_t teph, const nav_t *nav, aobsd_t *obs, int n, const double *rr,
                  int f, const prcopt_t *opt)
{
    static aobsd_t pobs[MAXOBS];
    static double rs[6*MAXSAT],dts[2*MAXSAT],vare[MAXSAT],r0[MAXSAT];
    static int svh[MAXSAT];
    double pos[3],e[3];
    double ens,amb,ambf,rns,P,min_cn0=30.0;
    const double *rrp,nc=0.299792458;
    char prn[8],prn_ref[8];
    int i,j,week,nfix=0,rf=f,flag=1,tag[MAXOBS]={0};
    gtime_t gtt,gtt_;

    trace(3,"fixambtt: n=%d f=%d time=%8.4lf\n",n,f,time2gpst(obs[0].time,NULL));

#if !FIX_OBS_AMB
    return 0;
#else
    ecef2pos((rr==NULL||rr[0]==0.0)?(rrp=epos):(rrp=rr),pos);
    trace(3,"pos=%.6lf %.6lf %.3lf\n",pos[0]*R2D,pos[1]*R2D,pos[2]);

    if (!nav->n&&!nav->ng) return 0;
    if (n<=1) return 0;

    memcpy(pobs,obs,sizeof(obsd_t)*n);

#if FIX_OBS_GLO_FCN
    /* fix glonass FCN ambiguity */
    fixglofcn(nav,pobs,n,rrp,f,opt);
#endif

re_selrefsat:
    /* search fix amb ref satellite */
    for (j=-1,i=0;i<3;i++) {
        j=selrefsat(pobs,tag,n,rf,flag);

        /* search fail */
        if (j<0) {
            if      (i==0) {rf=0; flag=1; continue;}
            else if (i==1) {rf=f; flag=0; continue;}
            else if (i==2) {rf=0; flag=0; continue;}
        }
        else {
            /* search ok */
            break;
        }
    }
    /* no valid ref satellite */
    if (j<0) {
        trace(2,"no reference sat to fix obs amb\n");

        nfix=fixambtt_(nav,pobs,n,rrp,rs,dts,vare,svh,r0,&j,&rf,f,opt);
        goto fixamb_ret;
    }
    if (f==0) {
        /* satellite position */
        satpvt(pobs[j].time,nav,pobs,n,rs,dts,svh,vare,opt);

        /* geometric distance/reference sat*/
        for (i=0;i<n;i++) {
            if ((r0[i]=geodist(rs+9*i,rrp,e))<0.0) {
                continue;
            }
        }
    }
    /* check satellite status */
    if (svh[j]||!dts[3*j]||vare[j]>MAX_VAR_EPH||pobs[j].SNR[rf]*0.25<min_cn0) {
        tag[j]=1;

        /* re-search ref-satellite */
        goto re_selrefsat;
    }
    time2gpst(pobs[j].time,&week);

    satno2id(pobs[j].sat,prn);
    trace(3,"fix amb: ref sat=%s SNR=%.3lf rf=%d\n",prn,pobs[j].SNR[rf]*0.25,rf);

    /* fix obs ambiguity */
    for (i=0;i<n;i++) {
        if (pobs[i].unc.tt[f]==0.0||r0[i]<=0.0||i==j) continue;
        if (!(pobs[i].stat[f]&OBSQ_PR_FIXAMB)) continue;

        /* solve 1ms's integers */
        ens=pobs[j].unc.tt[rf]+(r0[j]-r0[i])/nc;
        ambf=(ens-pobs[i].unc.tt[f])/1E6-(dts[3*j]-dts[3*i])*1E3;
        amb=round(ambf);
        rns=amb*1E6+pobs[i].unc.tt[f];

        /* fractional part */
        pobs[i].unc.ambf[f]=round(ambf)-ambf;

        if (fabs(ambf-round(ambf))>1E-1) {
            satno2id(pobs[i].sat,prn);
            trace(2,"ambf's value fail: sat=%s ambf=%.3lf\n",prn,ambf);
            continue;
        }
        gtt_=gpst2time(week,pobs[i].unc.tt[f]*1E-9);
        gtt=gpst2time(week,rns*1E-9);

        /* full pseudorange */
        if ((P=timediff(pobs[i].time,gtt)*CLIGHT)<=0.0) continue;

        if (!pobs[i].P[f]) {
            pobs[i].P[f]=P; /* fix */
        }
        else if (obs[i].P[f]>1E7&&obs[i].P[f]<4E7) {
            if (fabs(pobs[i].P[f]-P)>1E2) {
                trace(2,"fix amb fail: sat=%3d f=%d dP=%.3lf\n",pobs[i].sat,f,pobs[i].P[f]-P);
                continue;
            }
            pobs[i].P[f]=P;
        }
        else {
            pobs[i].P[f]=P;
        }
        /* fix uncertainty ambiguity */
        if (pobs[i].unc.rt[f]>1E3) {
            pobs[i].unc.rt[f]=10.0;
        }
        nfix++;

        satno2id(pobs[j].sat,prn_ref);
        satno2id(pobs[i].sat,prn);
        trace(3,"fix amb: sat=%3d(%3s) f=%d dt=%8.3lfms dt_=%15.3lf "
                "ambf=%15.4lf pr=%10.4lf ref=%3s\n",
                pobs[i].sat,prn,f,timediff(pobs[i].time,gtt)*1E3,timediff(pobs[i].time,gtt_)*1E3,
                ambf,pobs[i].P[f],prn_ref);
    }
fixamb_ret:
    /* restore observation */
    memcpy(obs,pobs,sizeof(obsd_t)*n);

    /* align f2/f3 1ms's ambiguity based on f1 */
    if (f>0) fixf2f3amb(nav,obs,n,f);

    /* check observation ambiguity */
    chk_obsamb(nav,obs,n,f,j,rf,r0,dts);
    return nfix;
#endif
}
/* select glonass ephememeris ------------------------------------------------*/
static geph_t *selgeph(gtime_t time, int sat, int iode, const nav_t *nav)
{
    double t,tmax=MAXDTOE_GLO,tmin=tmax+1.0;
    int i,j=-1;

    for (i=0;i<nav->ng;i++) {
        if (nav->geph[i].sat!=sat) continue;
        if (iode>=0&&nav->geph[i].iode!=iode) continue;
        if ((t=fabs(timediff(nav->geph[i].toe,time)))>tmax) continue;
        if (iode>=0) return nav->geph+i;
        if (t<=tmin) {j=i; tmin=t;} /* toe closest to time */
    }
    if (iode>=0||j<0) return NULL;
    return nav->geph+j;
}
/* fix glonass fcn ambiguity--------------------------------------------------
 * args:  obsd_t *obs   I  observation data
 *        int n         I  number of observation data
 *        nav_t *nav    I  navigation data
 *        double *rr    I  receriver position (ecef)
 *        int f         I  frequency no.
 *        prcopt_t *opt I  process options
 * note: since two glonass satellites transmit navigation signals on the same
 * carrier frequency, assisted data is needed to identify the correct osn.
 * ---------------------------------------------------------------------------*/
extern int fixglofcn(const nav_t *nav, aobsd_t *obs, int n, const double *rr, int f,
                     const prcopt_t *opt)
{
    double rs1[3],rs2[3],dr1,dr2,e[3];
    int i,j,k,prn1,prn2;
    geph_t *geph1,*geph2;

    trace(3,"fixglofcnamb: n=%d f=%d\n",n,f);

    for (i=j=0;i<n;i++) {
        if (!(obs[i].stat[f]&OBSQ_GLO_FCN)) continue;

        /* get glonass fcn */
        k=obs[i].unc.svid-100;
        if (k<-7||k>6) continue;

        /* two glonass prn corresponding to fcn:k */
        prn1=fcntbl[(k+7)*2+0];
        prn2=fcntbl[(k+7)*2+1];
        if (prn1==prn2) {
            obs[i].sat=satno(SYS_GLO,prn1);
            continue;
        }
        /* select glonass ephememeris */
        geph1=selgeph(obs[i].time,satno(SYS_GLO,prn1),-1,nav);
        geph2=selgeph(obs[i].time,satno(SYS_GLO,prn2),-1,nav);

        if (!geph1&&!geph2) continue;
        if (geph1&&geph2) goto fixglofcn;

        if (geph1&&!geph2) {
            obs[i].sat=satno(SYS_GLO,prn1);
            continue;
        }
        if (geph2&&!geph1) {
            obs[i].sat=satno(SYS_GLO,prn2);
            continue;
        }
fixglofcn:
        /* fix glonass fcn ambiguity */
        geph2pos(timeadd(obs[i].time,-MIN(MAX(obs[i].P[f],18E6),40E6)/CLIGHT),geph1,rs1,NULL,NULL);
        geph2pos(timeadd(obs[i].time,-MIN(MAX(obs[i].P[f],18E6),40E6)/CLIGHT),geph2,rs2,NULL,NULL);

        dr1=geodist(rs1,rr,e)-obs[i].P[f];
        dr2=geodist(rs2,rr,e)-obs[i].P[f];
        if (fabs(dr1)<fabs(dr2)) {
            obs[i].sat=satno(SYS_GLO,prn1);
        }
        else {
            obs[i].sat=satno(SYS_GLO,prn2);
        }
        /* # of glonass fcn fix */
        char prn[8];
        satno2id(obs[i].sat,prn);
        trace(3,"fix glonass FCN: %s\n",prn);
        j++;
    }
    return j;
}
