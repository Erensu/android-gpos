/*------------------------------------------------------------------------------
 * spplsq.c : extend standard positioning using robust wlsq
 *
 * author   : sujinglan
 * version  : $Revision:$ $Date:$
 * history  : 2022/03/05 1.0  new
 *-----------------------------------------------------------------------------*/
#include "agpos.h"

/* constants -----------------------------------------------------------------*/
#define SQR(x)       ((x)*(x))

#define CORR_COARSE_TIME 1              /* correct effects of coarse observation time */

#if CORR_COARSE_TIME
#define NPOS          3
#define NCLK          (NFREQ*4)
#define NVEL          3
#define NCLKD         1
#define NDT           1
#else
#define NPOS          3
#define NCLK          (NFREQ*4)
#define NVEL          3
#define NCLKD         1
#define NDT           0
#endif
#define NX_P          (NPOS+NCLK+NDT)   /* # of estimated parameters in pseudorange LSQ */
#define NX_D          (NVEL+NCLKD)      /* # of estimated parameters in doppler LSQ */
#define MAXITR        10                /* max number of iteration for point pos */
#define ERR_CBIAS     0.3               /* code bias error std (m) */
#define ADJ_SIG0      1                 /* adjust position/velocity sigma0 */
#define VAR_FIXALT    SQR(3.0)          /* constraint of altitude (m^2) */
#define MIN_SIG0      0.1               /* min position/velocity sigma0 */
#define MIN_SIG1_POS  100.0             /* min sigma0 for position solution */
#define MIN_SIG1_VEL  15.0              /* min sigma0 for velocity solution */
#define MIN_DOPS_POS  6.0               /* min DOPS for position solution */
#define MIN_DOPS_VEL  6.0               /* min DOPS for velocity solution */

#define AOBS_I(val)  ((val>>16)&0xFF)   /* observation index */
#define AOBS_F(val)  ((val>> 8)&0xFF)   /* observation frq */
#define AOBS_T(val)  ((val>> 4)&0xF)    /* observation type:
                                         * 1:code, 2:doppler, 3:position, 5: clock constraint
                                         * 7: constraint to geodetic height
                                         * */
#define AEXC_I(i,f)  (NFREQ*(i)+(f))    /* measurement index */

#define AIC(s,f)     (NPOS+(s)*NFREQ+(f))/* receiver clock index  */
#define AIC_(s,f)    ((s)*NFREQ+(f))
#define AIT          (NPOS+NCLK)        /* coarse observation time index */

static double mpos[3]={                 /* default receiver position for fix geodetic height */
        40.0500502920*D2R,              /* lat (rad) */
        116.277235789*D2R,              /* lon (rad) */
        30.0                            /* height (m) */
};

/* satellite positions and clocks for Android GNSS positioning------------------*/
extern void satpossa(gtime_t teph, const aobsd_t *obs, int n, const nav_t *nav,
                     int ephopt, double *rs, double *dts, double *var, int *svh);

/* update default receiver position for fix geodetic height-------------------*/
static void updmpos(const sol_t *sol)
{
    if (sol->stat!=SOLQ_SINGLE) return;
    double pos[3];

    ecef2pos(sol->rr,pos);
    mpos[0]=pos[0];
    mpos[1]=pos[1];
}
/* get group delay parameter (m) ---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);

    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) if (nav->geph[i].sat==sat) break;
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) if (nav->eph[i].sat==sat) break;
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* psendorange with code bias correction -------------------------------------*/
static double prange(const aobsd_t *obs, const nav_t *nav, const prcopt_t *opt, double *var, int f)
{
    double P,PC,gamma,beta,b1;
    int sat,sys;

    sat=obs->sat; sys=satsys(sat,NULL);
    PC=P=obs->P[f];

    if (!P) return 0.0;
    *var=SQR(ERR_CBIAS);

    if (sys==SYS_GPS||sys==SYS_QZS) {
        if (obs->code[f]==CODE_L1C) {
            P+=nav->cbias[sat-1][1]; /* C1->P1 */
        }
        else if (obs->code[f==CODE_L2C]) {
            P+=nav->cbias[sat-1][1]; /* C2->P2 */
        }
        gamma=SQR(FREQ1/sat2freq(obs->sat,obs->code[f],nav));
        PC=P-gamma*gettgd(obs->sat,nav,0);
    }
    else if (sys==SYS_GLO) {
        if (obs->code[f]==CODE_L1C) {
            P+=nav->cbias[sat-1][1]; /* C1->P1 */
        }
        else if (obs->code[f]==CODE_L2C) {
            P+=nav->cbias[sat-1][1]; /* C2->P2 */
        }
        beta=SQR(FREQ1_GLO)/sat2freq(obs->sat,obs->code[f],nav);
        gamma=SQR(FREQ1_GLO/FREQ2_GLO);
        b1=gettgd(sat,nav,0); /* -dtaun (m) */
        return P-b1*beta/(gamma-1.0);
    }
    else if (sys==SYS_GAL) { /* E1 */
        if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
        else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */

        gamma=SQR(FREQ1/sat2freq(obs->sat,obs->code[f],nav));
        PC=P-gamma*b1;
    }
    else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
        if      (obs->code[f]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
        else if (obs->code[f]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
        else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
        PC=P-b1;
    }
    return PC;
}
/* pseudorange residuals -----------------------------------------------------*/
static int prng_res(const aobsd_t *obs, int n, const double *rs, const double *dts,
                    const double *vare, const int *svh, const nav_t *nav, const prcopt_t *opt,
                    const double *x, double *v, double *H, double *var, double *azel,
                    int *vsat, double *resp, const int *exc, int *ns, int *ind)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],dtr,e[3],P,lam_L1;
    double sat_pr,freq;
    int i,j,f,nv=0,sys;

    trace(3,"prng_res :n=%d\n",n);

    for (i=0;i<3;i++) rr[i]=x[i]; dtr=x[3];
    ecef2pos(rr,pos);

    for (f=0;f<NFREQ;f++) {
        for (*ns=i=0;i<n&&i<MAXOBS;i++) {
            azel[i*2]=azel[1+i*2]=0.0;
            resp[AEXC_I(i,f)]=vsat[AEXC_I(i,f)]=0;

            if (obs[i].type!=OBST_CODE||!(sys=satsys(obs[i].sat,NULL))) continue;
            if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;

            /* geometric distance/az/el angle */
            if ((r=geodist(rs+i*9,rr,e))<=0.0) continue;
            satazel(pos,e,azel+i*2);

            /* code bias correction */
            if ((P=prange(obs+i,nav,opt,&vmeas,f))<=RE_WGS84) continue;

            /* ionospheric corrections */
            ionocorr(obs[i].time,nav,obs[i].sat,pos,azel+i*2,opt->ionoopt,&dion,&vion);

            if (!(freq=sat2freq(obs[i].sat,obs[i].code[f],nav))) continue;
            dion*=SQR(FREQ1/freq);
            vion*=SQR(FREQ1/freq);

            /* tropospheric corrections */
            tropcorr(obs[i].time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp);

            /* pseudorange residual */
            v[nv]=P-(r-CLIGHT*dts[i*3]+dion+dtrp);

            if (CORR_COARSE_TIME) {

                /* geometric range rate and satellite-clock rate only */
                matmul("NT",1,1,3,1.0,e,rs+i*9+3,0.0,&sat_pr);
                sat_pr-=dts[3*i+1]*CLIGHT;
                v[nv]-=sat_pr*x[AIT];
            }
            /* design matrix */
            for (j=0;j<NX_P;j++) H[j+nv*NX_P]=0.0;
            for (j=0;j<3;j++) H[j+nv*NX_P]=-e[j];

            /* receiver clock bias offset correction */
            if      (sys==SYS_GPS) {v[nv]-=x[AIC(0,f)]; H[AIC(0,f)+nv*NX_P]=1.0;}
            else if (sys==SYS_QZS) {v[nv]-=x[AIC(0,f)]; H[AIC(0,f)+nv*NX_P]=1.0;}
            else if (sys==SYS_GLO) {v[nv]-=x[AIC(1,f)]; H[AIC(1,f)+nv*NX_P]=1.0;}
            else if (sys==SYS_GAL) {v[nv]-=x[AIC(2,f)]; H[AIC(2,f)+nv*NX_P]=1.0;}
            else if (sys==SYS_CMP) {v[nv]-=x[AIC(3,f)]; H[AIC(3,f)+nv*NX_P]=1.0;}
            else continue;

            if (CORR_COARSE_TIME) {
                H[AIT+nv*NX_P]=sat_pr;
            }
            resp[AEXC_I(i,f)]=v[nv];

            /* check observation CN0 */
            if (obs[i].SNR[f]*0.25<opt->min_cn0_code) continue;

            /* check satellite cut-off elevation */
            if (azel[2*i+1]<opt->elmincode) continue;

            /* excluded satellite? */
            if (exc[AEXC_I(i,f)]) continue;

            vsat[AEXC_I(i,f)]=1;

            /* error variance */
            var[nv]=codevarerr(opt,obs[i].SNR[f]*0.25,azel[2*i+1],obs[i].sat)+vare[i]+vmeas+vion+vtrp;
            ind[nv]=((i<<16)|(f<<8)|(1<<4));
            nv++; (*ns)++;

            char prn[8];
            satno2id(obs[i].sat,prn);

            trace(3,"%4s sat=%3d azel=%5.1f %5.1f res=%10.3f sig=%8.3f cn0=%8.3lf f=%d stat=%4d\n",
                  prn,obs[i].sat,azel[i*2]*R2D,azel[1+i*2]*R2D,
                  v[nv-1],sqrt(var[nv-1]),obs[i].SNR[f]*0.25,f,obs[i].unc.stat[f]);
        }
    }
    return nv;
}
/* constraint to avoid rank-deficient----------------------------------------*/
static int ctr_res(const aobsd_t *obs, const double *x, double *v, double *H, double *var,
                   int *ind, int nv)
{
    int i,j,f,sys,mask[NX_P]={0};

    for (i=0;i<nv;i++) {
        sys=satsys(obs[AOBS_I(ind[i])].sat,NULL);
        f=AOBS_F(ind[i]);
        if      (sys==SYS_GPS) mask[AIC(0,f)]=1;
        else if (sys==SYS_QZS) mask[AIC(0,f)]=1;
        else if (sys==SYS_GLO) mask[AIC(1,f)]=1;
        else if (sys==SYS_GAL) mask[AIC(2,f)]=1;
        else if (sys==SYS_CMP) mask[AIC(3,f)]=1;
        else continue;
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) {
            if (mask[AIC(i,f)]) continue;
            for (j=0;j<NX_P;j++) H[j+nv*NX_P]=j==AIC(i,f)?1.0:0.0;
            var[nv]=SQR(1E-8);
            v[nv]=0.0;
            ind[nv++]=((i<<16)|(0<<8)|(5<<4));
            trace(3,"clock constraint: s=%d f=%d var=%.2lf\n",i,f,sqrt(var[nv-1]));
        }
    }
    return nv;
}
/* constraint to coarse observation time-------------------------------------*/
static int ctr_cobst(const aobsd_t *obs, const double *x, double *v, double *H, double *var,
                     int *ind, int nv)
{
#if CORR_COARSE_TIME
    if (obs[0].stat[0]&OBSQ_COARSE_TIME) {
        return nv;
    }
    else {
        for (int i=0;i<NX_P;i++) H[i+nv*NX_P]=i==AIT?1.0:0.0;
        v[nv]=0.0;
        var[nv]=SQR(1E-12);
        ind[nv]=((nv<<16)|(0<<8)|(6<<4));
        nv++;
        trace(3,"coarse obs time constraint: nv=%d\n",nv);
        return nv;
    }
#else
    return nv;
#endif
}
/* constraint to geodetic height---------------------------------------------*/
static int ctr_rcvhgt(const double *x, double *v, double *H, double *var,
                      int *ind, int nv)
{
    if (nv>NX_P) return nv;

    double s[3],pos[3];
    int i;

    trace(3,"constraint rcv height: hgt=%.3lf\n",mpos[2]);

    ecef2pos(x,pos);

    s[0]=cos(mpos[0])*cos(mpos[1]);
    s[1]=cos(mpos[0])*sin(mpos[1]);
    s[2]=sin(mpos[0]);
    for (i=0;i<3;i++) {
        H[i+nv*NX_P]=-s[i];
    }
    v[nv]=pos[2]-mpos[2];

    ind[nv]=((nv<<16)|(0<<8)|(7<<4));
    var[nv++]=VAR_FIXALT;
    return nv;
}
/* measurement residuals-----------------------------------------------------*/
static int meas_res(const aobsd_t *obs, int n, const double *rs, const double *dts,
                    const double *vare, const int *svh, const nav_t *nav, const prcopt_t *opt,
                    const double *x, double *v, double *H, double *var,
                    double *azel, int *vsat, double *resp, const int *exc, int *ns, int *ind)
{
    /* psendorange measurement */
    int nv=prng_res(obs,n,rs,dts,vare,svh,nav,opt,x,v,H,var,azel,vsat,resp,exc,ns,ind);

    /* clock constraint */
    nv=ctr_res(obs,x,v,H,var,ind,nv);

    /* coarse obs time constraint */
    nv=ctr_cobst(obs,x,v,H,var,ind,nv);

    /* constraint to geodetic height */
    nv=ctr_rcvhgt(x,v,H,var,ind,nv);
    return nv;
}
/* initial estimated states--------------------------------------------------*/
static void init_state(const aslsq_t *sol, const prcopt_t *opt, double *x)
{
    int i,f;
    for (i=0;i<3;i++) x[i]=sol->rr[i];
    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) x[AIC(i,f)]=sol->dtr[i*NFREQ+f];
    }
}
/* compute covariance matrix of post-residuals-------------------------------*/
static int post_res_covar(const double *H, const double *var, int nv, int nx, double *Rp)
{
    double *R=zeros(nv,nv),*HRHT;
    int i,j;

    for (i=0;i<nv;i++) R[i+i*nv]=var[i];
    if (matinv(R,nv)) {
        matcpy(Rp,R,nv,nv);
        free(R); return 0;
    }
    HRHT=mat(nx,nx);
    matmul33("NNT",H,R,H,nx,nv,nv,nx,HRHT);

    if (matinv(HRHT,nx)) {
        matcpy(Rp,R,nv,nv);
        free(R); free(HRHT);
        return 0;
    }
    matmul33("TNN",H,HRHT,H,nv,nx,nx,nv,R);
    for (i=0;i<nv;i++) {
        for (j=0;j<nv;j++) {
            Rp[i+j*nv]=i==j?MAX(var[i]-R[i+i*nv],1E-6):-R[i+j*nv];
        }
    }
    free(R); free(HRHT);
    return 1;
}
/* check redundancy matrix-----------------------------------------------------*/
static int redund_mat(const double *Rp, const double *var, int nv, double *Ru,
                      int *ruflg)
{
    double thres=1.0;
    int i,j;

    /* redundancy matrix */
    for (i=0;i<nv;i++) {
        for (j=0;j<nv;j++) Ru[i+j*nv]=Rp[i+j*nv]/var[j];
    }
    /* redundancy flag */
    for (i=0;i<nv;i++) {
        for (ruflg[i]=0,j=0;j<nv;j++) {
            if (fabs(Ru[i+i*nv])*thres<fabs(Ru[j+i*nv])) {
                ruflg[i]=1; break;
            }
        }
    }
    return 1;
}
/* index navi system (m=0:gps/sbs/qzss,1:glo,2:gal,3:bds) --------------------*/
static int ind_sys(int sys)
{
    switch (sys) {
        case SYS_GPS: return 0;
        case SYS_GLO: return 1;
        case SYS_GAL: return 2;
        case SYS_CMP: return 3;
        case SYS_QZS: return 0;
        default:
            return -1;
    }
}
/* compare pseudorange residuals----------------------------------------------*/
static int cmpres(const void *p1, const void *p2)
{
    double *q1=(double*)p1,*q2=(double*)p2;
    return ((*q1-*q2)<0)?-1:1;
}
/* use box plot to detect outliers--------------------------------------------*/
static int boxplot_detoutl(const double *v, int nv, int *outl_ind, double thres)
{
    double Q1,Q3,f,Qmin,Qmax,*vv;
    int i1,i2,i,outl_n=0;

    if (nv<5) return 0;
    vv=mat(nv,1); matcpy(vv,v,nv,1);

    qsort(vv,nv,sizeof(double),cmpres);
    i1=(nv+1)/4;
    if ((i2=i1+1)>=nv) i2=i1;
    f=(nv+1.0)/4.0-i1;
    Q1=(1.0-f)*vv[i1-1]+f*vv[i2-1]; /* lower quartile */

    i1=3*(nv+1)/4;
    if ((i2=i1+1)>=nv) i2=i1;
    f=3.0*(nv+1.0)/4.0-i1;
    Q3=(1.0-f)*vv[i1-1]+f*vv[i2-1]; /* upper quartile */

    Qmin=Q1-thres*(Q3-Q1); /* lower limit */
    Qmax=Q3+thres*(Q3-Q1); /* upper limit */

    trace(3,"boxplot: Q1=%6.3lf Q3=%6.3lf Qmin=%6.3lf Qmax=%6.3lf\n",Q1,Q3,Qmin,Qmax);

    for (i=0;i<nv;i++) {
        if (v[i]<Qmin||v[i]>Qmax) {
            trace(2,"exclude sat (boxplot): index=%d\n",i);

            /* index of outlier */
            if (outl_ind) outl_ind[outl_n++]=i;
        }
    }
    free(vv); return outl_n;
}
/* detect outliers through mad------------------------------------------------*/
static int mad_detoutl(const double *v, int nv, int *outl_ind, double thres)
{
    int i,outl_n=0;
    double m,med;

    if (nv<3) return 0;

    med=median(v,nv); m=mad(v,nv);
    for (i=0;i<nv;i++) {
        if (fabs((v[i]-med)/(1.4826*m))<thres) continue;
        if (outl_ind) outl_ind[outl_n++]=i;
        trace(2,"exclude sat (mad): index=%d\n",i);
    }
    return outl_n;
}
/* detect outliers method----------------------------------------------------*/
static int detoutl_method(const double *v, int nv, int *outl_ind, double thres_box,
                          double thres_mad)
{
    int outl_n=0;
    if (thres_box) {
        outl_n+=boxplot_detoutl(v,nv,outl_ind,thres_box);
    }
    if (thres_mad) {
        outl_n+=mad_detoutl(v,nv,outl_ind+outl_n,thres_mad);
    }
    return outl_n;
}
/* detect measurement outliers-----------------------------------------------*/
static int detoutl(const double *v, int nv, int *outl_ind, double thres_box,
                   double thres_mad)
{
    return detoutl_method(v,nv,outl_ind,thres_box,thres_mad);
}
/* detect measurement outliers-----------------------------------------------*/
static int detoutl_prc(const double *v, int nv, int *outli, double thres_box,
                       double thres_mad)
{
    if (nv<=0) return 0;

    double *vk=mat(1,nv),med;
    int i,n,outln=0;

    med=median(v,nv);
    for (i=n=0;i<nv;i++) vk[n++]=v[i]-med;

    /* detect outliers */
    outln+=detoutl(vk,n,outli,thres_box,thres_mad);
    free(vk);
    return outln;
}
/* detect code measurement outliers------------------------------------------*/
static void detoutl_meas_code(const aobsd_t *obs, const double *v, const int *ind, int nv, int *exc,
                              const double *x, const prcopt_t *opt,
                              double thres_box, double thres_mad)
{
    int i,f,j,k,mi[4*NFREQ][MAXOBS]={0},ms[4*NFREQ]={0},*outli,outln;
    double *vm[4*NFREQ]={0};

    trace(3,"detoutl_meas_code: nv=%d\n",nv);

    outli=imat(2*nv,1);

    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) vm[AIC_(i,f)]=mat(nv,1);
    }
    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) {
            for (j=0;j<nv;j++) {
                if (AOBS_T(ind[j])!=1||AOBS_F(ind[j])!=f) continue;
                if (ind_sys(satsys(obs[AOBS_I(ind[j])].sat,NULL))!=i) continue;
                mi[AIC_(i,f)][ms[AIC_(i,f)]]=j;
                vm[AIC_(i,f)][ms[AIC_(i,f)]++]=v[j];
            }
        }
    }
    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) {
            if (!ms[AIC_(i,f)]) continue;

            /* exclude less satellite */
            if (ms[AIC_(i,f)]<=2) {
                for (outln=j=0;j<ms[AIC_(i,f)];j++) outli[outln++]=j;
            }
            else {
                /* detect measurement outliers */
                outln=detoutl_prc(vm[AIC_(i,f)],ms[AIC_(i,f)],outli,thres_box,thres_mad);
            }
            for (j=0;j<outln;j++) {
                if ((k=mi[AIC_(i,f)][outli[j]])<0) continue;
                if (exc[AEXC_I(AOBS_I(ind[k]),AOBS_F(ind[k]))]) continue;
                exc[AEXC_I(AOBS_I(ind[k]),AOBS_F(ind[k]))]=1; /* exclude satellite */

                char prn[8];
                satno2id(obs[AOBS_I(ind[k])].sat,prn);
                trace(2,"exclude sat=%s f=%d\n",prn,AOBS_F(ind[k]));
            }
        }
    }
    free(outli);
    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) free(vm[AIC_(i,f)]);
    }
}
/* process WLSQ--------------------------------------------------------------*/
static int awlsq(double *v, double *H, int nv, int nx, const double *var,
                 double *dx, double *Q)
{
    double sig;
    int i,j;

    /* 	normalization of v/H/R */
    for (i=0;i<nv;i++) {
        sig=sqrt(var[i]);
        v[i]/=sig;
        for (j=0;j<nx;j++) H[j+i*nx]/=sig;
    }
    /* least square estimation */
    if (lsq(H,v,nx,nv,dx,Q)) return -1;
    return 1;
}
/* DOPS check-----------------------------------------------------------------*/
static int chk_dops(int type, const int *ind, int nv, const double *azel, const aobsd_t *obs,
                    const prcopt_t *opt, double *dop)
{
    int sats[MAXOBS],i,j,k,ns,upd_flg=0;
    double azels[2*MAXOBS];

    /* large DOPS check */
    for (i=j=ns=0;i<nv;i++) {
        if (AOBS_T(ind[i])!=type) continue;
        azels[  ns*2]=azel[  AOBS_I(ind[i])*2];
        azels[1+ns*2]=azel[1+AOBS_I(ind[i])*2];

        for (upd_flg=1,k=0;k<j;k++) {
            if (obs[AOBS_I(ind[i])].sat==sats[k]) {
                upd_flg=0; break;
            }
        }
        if (upd_flg) {
            sats[j++]=obs[AOBS_I(ind[i])].sat;
            ns++;
        }
    }
    dops(ns,azels,type==1?opt->elmincode:opt->elmindopp*D2R,dop);
    trace(3,"DOP=%6.3lf %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf\n",
          dop[0],dop[1],dop[2],dop[3],dop[4],dop[5]);

    if (dop[0]<=0.0||dop[0]>opt->maxgdop) return 0;
    return 1;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const double *v, const double *Q, const aobsd_t *obs,
                  const int *ind, int nv, int nx,
                  const prcopt_t *opt, aslsq_t *lsq)
{
    double vv,rr;

    /* check position variance/residual */
    if (lsq->sig0_p0*sqrt(Q[0]+Q[1+NX_P]+Q[2+2*NX_P])>300.0||lsq->sig0_p1>30.0) {
        trace(2,"check position variance fail\n");
        return 0;
    }
    /* position DOPS check */
    if (!chk_dops(1,ind,nv,azel,obs,opt,lsq->dop_pos)) {
        trace(2,"position dops check fail\n");
        return 0;
    }
    /* chi-square valid of residuals */
    if (nv>nx) {
        vv=dot(v,v,nv);
        rr=vv/SQR(MIN(lsq->sig0_p0,1.0))/chisqr[nv-nx-1];
        if (rr>1.0) {
            trace(2,"chisqr error: vv=%.1f cs=%.1f\n",vv,chisqr[nv-nx-1]);
            return 0;
        }
    }
    trace(3,"valid position solution ok\n");
    return 1;
}
/* update solution------------------------------------------------------------*/
static void updpossol(aslsq_t *sol, const aobsd_t *obs, const double *x, const double *Q,
                      const int *ind, int nv)
{
    double Qxyz[9]={0},Qenu[9]={0},pos[3];
    double var_sig0;
    int i,j,f;

    /* update position solution */
    sol->time=timeadd(obs[0].time,-x[AIC(0,0)]/CLIGHT);

    for (i=0;i<4;i++) {
        for (f=0;f<NFREQ;f++) sol->dtr[AIC_(i,f)]=x[AIC(i,f)];
    }
    for (i=0;i<3;i++) sol->rr[i]=x[i];
    for (i=0;i<3;i++) sol->qr[i]=Q[i+i*NX_P];
    for (i=0;i<4;i++) sol->qt[i]=Q[i+i*NX_P];
    sol->qr[3]=Q[1+0*NX_P];
    sol->qr[4]=Q[2+1*NX_P];
    sol->qr[5]=Q[2+0*NX_P];

#if CORR_COARSE_TIME
    sol->dt=x[AIT];
#else
    sol->dt=0.0;
#endif

    sol->psta=1;
    for (sol->ns=i=0;i<nv;i++) {
        if (AOBS_T(ind[i])==1) sol->ns++;
    }
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) Qxyz[i+j*3]=SQR(sol->sig0_p0)*Q[i+j*NX_P];
    }
    ecef2pos(x,pos);
    covenu(pos,Qxyz,Qenu);
    tracemat(3,Qenu,3,3,12,5);

#if ADJ_SIG0
    /* variance of unit weight variance */
    var_sig0=2.0/(nv-NX_P)*SQR(SQR(sol->sig0_p0));

    /* adjust unit weight variance */
    if (sol->dop_pos[0]>MIN_DOPS_POS) {
        sol->sig0_p0=MAX(1.0,sol->sig0_p0);
    }
    if (sol->sig0_p1>MIN_SIG1_POS) {
        sol->sig0_p0=MAX(1.0,sol->sig0_p0);
    }
    if (sol->pdeg) {
        sol->sig0_p0=MAX(1.0,sol->sig0_p0);
    }
    if (nv-NX_P<2) {
        sol->sig0_p0=MAX(1.0,sol->sig0_p0);
    }
    if (sqrt(Qenu[0]+Qenu[4])>80.0) {
        sol->sig0_p0=MAX(1.0,sol->sig0_p0);
    }
    if (sqrt(var_sig0)>1.0) {
        sol->sig0_p0=MAX(1.0,sol->sig0_p0);
    }
    sol->sig0_p0=MAX(MIN_SIG0,sol->sig0_p0);
    sol->sig0_p0=MIN(1.0,sol->sig0_p0);
#endif
}
/* exclude bad satellite-----------------------------------------------------*/
static int exclude_sat(const double *v, const double *var, const double *Rp,
                       const int *redu, const int *ind, int nv, int *exc)
{
    int i,j;

    for (j=-1,i=0;i<nv;i++) {
        if (j<0||fabs(v[i]*sqrt(var[i])/sqrt(Rp[i+i*nv]))>fabs(v[j]*sqrt(var[j])/sqrt(Rp[j+j*nv]))) {
            j=i;
        }
    }
    if (j>=0) {
        exc[AEXC_I(AOBS_I(ind[j]),AOBS_F(ind[j]))]=1;
    }
    return j;
}
/* estimate receiver position main function----------------------------------*/
static int estpos(const aobsd_t *obs, const int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav, const prcopt_t *opt,
                  const assat_t *ssat, aslsq_t *lsq, double *azel, int *vsat, double *resp)
{
    int i,j=0,ns,nx=NX_P,try_cnt,stat=0,nv=3*n*NFREQ+NX_P,*exc,*ind,*ru;
    double *x,*v,*H,*var,*dx,*Q,*Rp,*Ru,*vp;
    aslsq_t lsq_;

    trace(3,"estpos: n=%d\n",n);

    ind=imat(nv,1); ru=imat(nv,1); vp=mat(nv,1);
    x=mat(nx,1); v=mat(1,nv); Ru=mat(nv,nv); H=zeros(nv,nx); exc=imat(nv,1); Q=mat(nx,nx);
    dx=mat(nx,1); Rp=mat(nv,nv); var=mat(nv,1);

    memset(vsat,0,sizeof(int)*NFREQ*n*2);

    /* initial estimated states */
    init_state(lsq,opt,x);
    lsq_=*lsq;

    /* measurement residuals */
    if ((nv=meas_res(obs,n,rs,dts,vare,svh,nav,opt,x,v,H,var,azel,vsat,resp,exc,&ns,ind))<nx) {
        trace(2,"no valid measurement\n");
        stat=-1;
        goto __ret_poslsq;
    }
    /* detect measurement outlier */
    detoutl_meas_code(obs,v,ind,nv,exc,x,opt,4.0,8.0);

    /* start robust WLSQ estimator */
    for (stat=0,try_cnt=0;try_cnt<n;try_cnt++) {
        init_state(&lsq_,opt,x);

        /* iteration for WLSQ */
        for (i=0;i<MAXITR;i++) {
            if ((nv=meas_res(obs,n,rs,dts,vare,svh,nav,opt,x,v,H,var,azel,vsat,resp,exc,&ns,ind))<NX_P) {
                trace(2,"no valid measurement\n");
                stat=-1;
                goto __ret_poslsq;
            }
            /* post-residuals covariance */
            post_res_covar(H,var,nv,nx,Rp);

            /* compute redundancy matrix */
            redund_mat(Rp,var,nv,Ru,ru);

            /* do WLSQ estimation */
            matcpy(vp,v,1,nv);
            if (awlsq(v,H,nv,nx,var,dx,Q)<0) {
                trace(1,"position WLSQ error\n");
                stat=-3;
                break;
            }
            /* update WLSQ solutions */
            for (j=0;j<nx;j++) x[j]+=dx[j];

            /* validation of solution */
            if (norm(dx,nx)<1E-2) {
                lsq->sig0_p1=sqrt(dot(vp,vp,nv)/(nv-nx));
                lsq->sig0_p0=sqrt(dot(v,v,nv)/(nv-nx));
                stat=valsol(azel,v,Q,obs,ind,nv,nx,opt,lsq);

                trace(3,"pos-sig1=%6.3lf\n",lsq->sig0_p1);
                trace(3,"pos-sig0=%6.3lf\n",lsq->sig0_p0);
                break;
            }
        }
        /* validate solution ok */
        if (stat==1) {
            updpossol(lsq,obs,x,Q,ind,nv);
            break;
        }
        /* do RAIM FDE */
        if ((j=exclude_sat(v,var,Rp,ru,ind,nv,exc))<0) {
            stat=-5;
            break;
        }
        char prn[8];
        satno2id(obs[AOBS_I(ind[j])].sat,prn);
        trace(2,"exclude sat=%s f=%d\n",prn,AOBS_F(ind[j]));
    }
__ret_poslsq:
    free(x); free(v); free(H); free(var); free(dx); free(Q);
    free(exc); free(Ru); free(vp);
    free(ind); free(ru); free(Rp);
    return stat;
}
/* detect code measurement outliers------------------------------------------*/
static void detoutl_meas_dopp(const aobsd_t *obs, const double *v, const int *ind, int nv, int *exc,
                              double thres_box, double thres_mad)
{
    int i,*outli,outln;
    char prn[8];

    trace(3,"detoutl_meas_dopp: nv=%d\n",nv);

    outli=imat(4*nv,1);

    outln=detoutl_prc(v,nv,outli,thres_box,thres_mad);
    for (i=0;i<outln;i++) {
        if (exc[AEXC_I(AOBS_I(ind[outli[i]]),AOBS_F(ind[outli[i]]))]) continue;
        exc[AEXC_I(AOBS_I(ind[outli[i]]),AOBS_F(ind[outli[i]]))]=1;

        satno2id(obs[AOBS_I(ind[outli[i]])].sat,prn);
        trace(2,"exclude sat=%s f=%d\n",prn,AOBS_F(ind[outli[i]]));
    }
    free(outli);
}
/* doppler residuals ---------------------------------------------------------*/
static int resdop(int post, const aobsd_t *obs, int n, const double *rs, const double *dts, const double *vare, const int *svh,
                  const prcopt_t *opt, const nav_t *nav, const assat_t *ssat,
                  const double *rr, const double *x, const double dt, const double *azel, double *v, double *H,
                  double *var, double *resd, int *vsat, const int *exc, int *ind)
{
    double rate,pos[3],E[9],a[3],e[3],vs[3],cosel,y;
    char prn[8];
    int i,j,f=0,nv=0,nx=NX_D;

    trace(3,"resdop: post=%d\n",post);

    ecef2pos(rr,pos); xyz2enu(pos,E);

    for (f=0;f<NFREQ;f++) {
        for (i=0;i<n&&i<MAXOBS;i++) {
            resd[AEXC_I(i,f)]=0.0; vsat[AEXC_I(i,f)]=0;

            /* get doppler measurement data */
            if (obs[i].type!=OBST_CODE||!(y=obs[i].Pr[f])) continue;
            if (norm(rs+3+i*9,3)<=0.0) continue;
            if (satexclude(obs[i].sat,vare[i],svh[i],opt)) continue;

            /* line-of-sight vector in ECEF */
            cosel=cos(azel[1+i*2]);
            a[0]=sin(azel[i*2])*cosel;
            a[1]=cos(azel[i*2])*cosel;
            a[2]=sin(azel[1+i*2]);
            matmul("TN",3,1,3,1.0,E,a,0.0,e);

            /* satellite velocity relative to receiver in ECEF */
            for (j=0;j<3;j++) vs[j]=rs[j+3+i*9]+rs[j+6+i*9]*dt-x[j];

            /* range rate with earth rotation correction */
            rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*9]*rr[0]+rs[1+i*9]*x[0]-rs[3+i*9]*rr[1]-rs[i*9]*x[1]);

            /* doppler measurement residual */
            v[nv]=y-(rate+x[3]-CLIGHT*(dts[1+i*3]+dts[2+i*3]*dt));

            /* check observation CN0 */
            if (obs[i].SNR[f]*0.25<opt->min_cn0_dopp) continue;

            /* check elevation mask angle */
            if (azel[1+i*2]<opt->elmindopp*D2R) continue;
            if (exc[AEXC_I(i,f)]) continue;

            /* measurement variance */
            var[nv]=doppvarerr(opt,obs[i].SNR[f]*0.25,azel[2*i+1],obs[i].sat);
            resd[AEXC_I(i,f)]=v[nv];

            satno2id(obs[i].sat,prn);

            trace(3,"%4s sat=%3d azel=%5.1f %5.1f res=%8.3f prr=%8.3lf sig=%8.3lf cn0=%8.3lf f=%d\n",
                  prn,obs[i].sat,azel[i*2]*R2D,azel[1+i*2]*R2D,v[nv],
                  y,sqrt(var[nv]),obs[i].SNR[f]*0.25,f);

            /* design matrix */
            for (j=0;j<nx;j++) H[j+nv*nx]=j<3?-e[j]:1.0;
            ind[nv++]=((i<<16)|(f<<8)|(2<<4));
            vsat[AEXC_I(i,f)]=1;
        }
    }
    return nv;
}
/* initial estimated states--------------------------------------------------*/
static void init_state_dopp(const aslsq_t *sol, double *x)
{
    int i;
    for (i=0;i<3;i++) x[i]=sol->rr[3+i];
    x[3]=sol->ddrt;
}
/* validation of doppler solution---------------------------------------------*/
static int valsol_dopp(const double *vr, const double *v, const double *azel, const int *ind,
                       const aobsd_t *obs, const prcopt_t *opt, int nv, int nx, const double *Q, aslsq_t *lsq)
{
    double vv,rr;

    trace(3,"valsol_dopp: n=%d nv=%d\n",nx,nv);

    /* check velocity variance/residual */
    if (lsq->sig0_v0*sqrt(Q[0]+Q[1+1*nx]+Q[2+2*nx])>10.0||lsq->sig0_v1>5.0) {
        trace(2,"check velocity variance fail\n");
        return 0;
    }
    /* check velocity value */
    if ((vv=norm(vr,3))>=1E2) {
        trace(2,"too large velocity: v=%.2lf\n",vv);
        return 0;
    }
    /* valid chi-square of residuals */
    if (nv>nx) {
        vv=dot(v,v,nv);
        rr=vv/SQR(MIN(lsq->sig0_v0,1.0))/chisqr[nv-nx-1];
        if (rr>1.0) {
            trace(2,"chi-square error: vv=%.1f cs=%.1f\n",vv,chisqr[nv-nx-1]);
            return 0;
        }
    }
    /* velocity DOPS check */
    if (!chk_dops(2,ind,nv,azel,obs,opt,lsq->dop_vel)) {
        trace(2,"velocity dops check fail\n");
        return 0;
    }
    trace(3,"velocity solution ok\n");
    return 1;
}
/* update velocity solution---------------------------------------------------*/
static void updvelsol(aslsq_t *sol, const double *x, const double *Q, int nv)
{
    double var_sig0;
    int i;

    /* update velocity solution */
    for (i=0;i<3;i++) sol->rr[i+3]=x[i];
    for (i=0;i<3;i++) sol->qv[i]=Q[i+i*4];
    sol->qv[3]=Q[ 1];
    sol->qv[4]=Q[ 6];
    sol->qv[5]=Q[ 2];
    sol->qt[4]=Q[15];
    sol->nd=nv; sol->vsta=1;
    sol->ddrt=x[3];

#if ADJ_SIG0
    /* variance of unit weight variance */
    var_sig0=2.0/(nv-NX_D)*SQR(SQR(sol->sig0_v0));

    /* adjust unit weight variance */
    if (sol->dop_vel[0]>MIN_DOPS_VEL) {
        sol->sig0_v0=MAX(1.0,sol->sig0_v0);
    }
    if (sol->sig0_v1>MIN_SIG1_VEL) {
        sol->sig0_v0=MAX(1.0,sol->sig0_v0);
    }
    if (sol->vdeg) {
        sol->sig0_v0=MAX(1.0,sol->sig0_v0);
    }
    if (sqrt(var_sig0)>1.0) {
        sol->sig0_v0=MAX(1.0,sol->sig0_v0);
    }
    sol->sig0_v0=MAX(MIN_SIG0,sol->sig0_v0);
    sol->sig0_v0=MIN(1.0,sol->sig0_v0);
#endif
}
/* estimate receiver velocity ------------------------------------------------*/
static int estvel(const aobsd_t *obs, int n, const double *rs, const double *dts, const double *vare,
                  const int *svh, const nav_t *nav, const prcopt_t *opt, aslsq_t *lsq,
                  const assat_t *ssat, const double *azel, int *vsatd, double *resd)
{
    int i,j=0,try_cnt,stat=0,nv=n*NFREQ,*exc,*ind,*ru,post;
    const int maxiter=n,nx=NX_D;
    double *x,*v,*H,*var,*dx,*Q,*Rp,*Ru,*vp,dt=lsq->dt;
    aslsq_t lsq_={0};

    trace(3,"estvel: n=%d\n",n);

    if (!lsq->psta) return 0;

    x=mat(nx,1); v=mat(1,nv); var=mat(nv,1); H=zeros(nv,nx);
    exc=imat(nv,1); ru=imat(nv,1);
    ind=imat(nv,1); Rp=mat(nv,nv); Q=mat(nx,nx); dx=mat(nx,1);
    Ru=mat(nv,nv); vp=mat(nv,1);

    memset(vsatd,0,sizeof(int)*NFREQ*n);

    /* initial estimated states */
    lsq_=*lsq;
    init_state_dopp(lsq,x);

    /* measurement residuals */
    if ((nv=resdop(-1,obs,n,rs,dts,vare,svh,opt,nav,ssat,lsq_.rr,x,dt,azel,v,H,var,resd,vsatd,exc,ind))<=4) {
        trace(2,"no valid measurement\n");
        stat=-1;
        goto vellsqquit;
    }
    /* detect meas outlier */
    detoutl_meas_dopp(obs,v,ind,nv,exc,3.0,6.0);

    /* start robust WLSQ estimator */
    for (stat=0,post=0,try_cnt=j;try_cnt<maxiter;try_cnt++) {
        init_state_dopp(&lsq_,x);

        for (i=0;i<MAXITR;i++) {

            /* measurement residuals */
            if ((nv=resdop(post,obs,n,rs,dts,vare,svh,opt,nav,ssat,lsq_.rr,x,dt,azel,v,H,var,resd,vsatd,exc,ind))<=NX_D) {
                trace(2,"no valid measurement\n");
                stat=-1;
                break;
            }
            /* post-residuals covariance */
            post_res_covar(H,var,nv,nx,Rp);

            /* compute redundancy matrix */
            redund_mat(Rp,var,nv,Ru,ru);

            /* do WLSQ estimation */
            matcpy(vp,v,1,nv);
            if (awlsq(v,H,nv,nx,var,dx,Q)<0) {
                trace(1,"velocity WLSQ error\n");
                stat=-3;
                break;
            }
            for (j=0;j<nx;j++) x[j]+=dx[j];
            if (norm(dx,nx)<1E-2) {

                /* validate solution */
                lsq->sig0_v1=sqrt(dot(vp,vp,nv)/(nv-nx));
                lsq->sig0_v0=sqrt(dot(v,v,nv)/(nv-nx));
                stat=valsol_dopp(x,v,azel,ind,obs,opt,nv,nx,Q,lsq);

                trace(3,"vel-sig0=%6.3lf\n",lsq->sig0_v0);
                trace(3,"vel-sig1=%6.3lf\n",lsq->sig0_v1);
                break;
            }
        }
        /* update velocity solution */
        if (stat==1) {
            updvelsol(lsq,x,Q,nv);
            break;
        }
        /* do raim fde */
        if ((j=exclude_sat(v,var,Rp,ru,ind,nv,exc))<0) {
            stat=-5;
            break;
        }
        char prn[8];
        satno2id(obs[AOBS_I(ind[j])].sat,prn);
        trace(2,"exclude sat=%s f=%d\n",prn,AOBS_F(ind[i]));
    }
vellsqquit:
    free(x); free(v); free(H); free(var); free(dx); free(vp);
    free(exc); free(Ru); free(Q);
    free(ind); free(ru); free(Rp);
    return stat;
}
/* update satellite status information --------------------------------------*/
static void udsatstat(assat_t *ssat, const double *rs, const double *dts, const double *var,
                      const double *azel, const int *vp, const int *vd, const double *resp,
                      const double *resd, const aobsd_t *obs, int n)
{
    int i,j;

    /* update satellite status information */
    for (i=0;i<MAXSAT;i++) {
        for (j=0;j<NFREQ;j++) {
            ssat[i].vs_lsq[j]=ssat[i].vd_lsq[j]=0;
            ssat[i].resp[j]=ssat[i].resd[j]=0.0;
        }
        ssat[i].azel[0]=0.0;
        ssat[i].azel[1]=0.0;
    }
    for (i=0;i<n;i++) {
        /* update satellite position/velocity/clocks */
        matcpy(ssat[obs[i].sat-1].azel,azel+2*i,1,2);
        matcpy(ssat[obs[i].sat-1].rs,rs+9*i,1,9);
        matcpy(ssat[obs[i].sat-1].dts,dts+3*i,1,3);
        ssat[obs[i].sat-1].var_s=var[i];

        /* update satellite valid flag and residuals */
        for (j=0;j<NFREQ;j++) {
            if (obs[i].SNR[j]) {
                ssat[obs[i].sat-1].SNR[j]=obs[i].SNR[j];
                ssat[obs[i].sat-1].time[j]=obs[i].time;
            }
            ssat[obs[i].sat-1].vs_lsq[j]=vp[AEXC_I(i,j)];
            ssat[obs[i].sat-1].vd_lsq[j]=vd[AEXC_I(i,j)];
            ssat[obs[i].sat-1].resp[j]=resp[AEXC_I(i,j)];
            ssat[obs[i].sat-1].resd[j]=resd[AEXC_I(i,j)];
        }
    }
}
/* update LSQ solution-------------------------------------------------------*/
static int updlsqsol(gtime_t time, aslsq_t *lsq, sol_t *sol)
{
    int i;

    trace(3,"updlsqsol: sow=%.3lf\n",time2gpst(time,NULL));

    if (!lsq->psta&&!lsq->vsta) {
        trace(2,"update lsq solution fail\n");
        return 0;
    }
    sol->type=0;
    sol->stat=SOLQ_SINGLE;
    sol->time=lsq->time;

    /* update position LSQ solution */
    if (lsq->psta) {
        for (i=0;i<3;i++) sol->rr[i]=lsq->rr[i];
        for (i=0;i<4*NFREQ;i++) sol->dtr[i]=lsq->dtr[i];
        matcpy(sol->qr,lsq->qr,1,6);
        sol->ns=lsq->ns;
    }
    else {
        for (i=0;i<3;i++) sol->rr[i]=lsq->rr[i];
        for (i=0;i<4*NFREQ;i++) sol->dtr[i]=lsq->dtr[i];
        for (i=0;i<6;i++) sol->qr[i]=0.0;
        matcpy(sol->qr,lsq->qr,1,6);
        sol->ns=0;
        lsq->sig0_p0=10.0;
        trace(3,"update lsq position solution fail\n");
    }
    /* update velocity LSQ solution */
    if (lsq->vsta) {
        for (i=3;i<6;i++) sol->rr[i]=lsq->rr[i];
        sol->ddrt=lsq->ddrt;
        matcpy(sol->qv,lsq->qv,1,6);
        sol->nd=lsq->nd;
    }
    else {
        for (i=3;i<6;i++) sol->rr[i]=lsq->rr[i];
        for (i=0;i<6;i++) sol->qv[i]=0.0;
        sol->ddrt=lsq->ddrt;
        sol->nd=0;
        lsq->sig0_v0=10.0;
        trace(3,"update lsq velocity solution fail\n");
    }
    trace(3,"update LSQ solution ok: pos=%8.3f %8.3f %8.3f vel=%6.3f %6.3f %6.3f\n",
          sol->rr[0],sol->rr[1],sol->rr[2],
          sol->rr[3],sol->rr[4],sol->rr[5]);
    return 1;
}
/* GNSS navigation LSQ using android raw gnss measurement data----------------
 * args:    agpos_t* agp   IO  Android GNSS position type
 *          obsd_t* obs    I   observation data
 *          int n          I   number of observation data
 *          nav_t* nav     I   navigation data
 * return : status(1:ok,0:error)
 * --------------------------------------------------------------------------*/
extern int aspplsq(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav)
{
    double *rs,*dts,*var,*azel,*resp,*resd;
    int i,j,nv=2*n*NFREQ,*vsatp,svh[MAXOBS],*vsatd,stat;
    aslsq_t *lsq=&agp->spp.lsq;
    assat_t *ssat=agp->ssat;
    prcopt_t opt_=agp->opt;

    trace(3,"aspplsq: sow=%6.3lf n=%d\n",time2gpst(obs[0].time,NULL),n);

    lsq->time=obs[0].time;
    lsq->psta=lsq->vsta=0;
    lsq->pdeg=lsq->vdeg=0;
    for (i=0;i<6;i++) lsq->qr[i]=0.0;
    for (i=0;i<6;i++) lsq->qv[i]=0.0;

    lsq->sig0_p0=lsq->sig0_v0=1.0;
    lsq->sig0_p1=lsq->sig0_v1=0.0;
    lsq->dt=0.0;
    agp->sol.stat=SOLQ_NONE;

    lsq->min_cn0[0]=agp->opt.min_cn0_code;
    lsq->min_cn0[1]=agp->opt.min_cn0_dopp;

    rs=mat(9,n); dts=mat(3,n); var=mat(1,n); vsatp=imat(1,nv); vsatd=imat(1,nv);
    azel=zeros(2,n); resp=mat(1,nv); resd=mat(1,nv);

    /* satellite position/velocity/clocks */
    satpossa(obs[0].time,obs,n,nav,agp->opt.sateph,rs,dts,var,svh);

    /* estimate receiver position */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&agp->opt,ssat,lsq,azel,vsatp,resp);
    if (stat!=1) {
        for (i=0;i<MAXITR;i++,lsq->pdeg++) {
            agp->opt.min_cn0_code-=3.0;
            lsq->min_cn0[0]=agp->opt.min_cn0_code;

            stat=estpos(obs,n,rs,dts,var,svh,nav,&agp->opt,ssat,lsq,azel,vsatp,resp);
            if (agp->opt.min_cn0_code<10.0) break;
            if (stat==1) break;
        }
    }
    /* estimate receiver velocity */
    stat=estvel(obs,n,rs,dts,var,svh,nav,&agp->opt,lsq,ssat,azel,vsatd,resd);
    if (stat!=1) {
        for (i=0;i<MAXITR;i++,lsq->vdeg++) {
            agp->opt.min_cn0_dopp-=3.0;
            lsq->min_cn0[1]=agp->opt.min_cn0_dopp;

            stat=estvel(obs,n,rs,dts,var,svh,nav,&agp->opt,lsq,ssat,azel,vsatd,resd);
            if (agp->opt.min_cn0_dopp<10.0) break;
            if (stat==1) break;
        }
    }
    /* update LSQ solution */
    updlsqsol(lsq->time,lsq,&agp->sol);

    /* update default receiver position */
    updmpos(&agp->sol);

    /* update satellite status information */
    udsatstat(ssat,rs,dts,var,azel,vsatp,vsatd,resp,resd,obs,n);

    agp->opt=opt_;
    free(dts); free(var); free(rs);
    free(vsatp); free(vsatd);
    free(azel); free(resp); free(resd);
    return agp->sol.stat;
}