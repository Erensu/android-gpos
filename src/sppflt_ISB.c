/*------------------------------------------------------------------------------
 * sppflt_ISB.c: GNSS navigation filter using Android raw gnss measurement data
 *
 * version    : $Revision: 1.2 $ $Date: 2008/07/14 00:05:05 $
 * history    : 2019/12/17 1.0  new
 *-----------------------------------------------------------------------------*/
#include "agpos.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define VAR_CDT     SQR(10.0)     /* init variance corrections of coarse observation time (s^2) */
#define VAR_CBIAS   SQR(0.3)      /* code bias error std (m) */

#define ERR_SAAS    0.3           /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5           /* broadcast iono model error factor */
#define ERR_CBIAS   0.3           /* code bias error std (m) */
#define REL_HUMI    0.7           /* relative humidity for saastamoinen model */
#define CORR_TIME_ISB   3600.0    /* correlation time of inter-system clock bias */
#define CORR_TIME_CLKD  120.0     /* correlation time of clock drift */

#define CORR_COARSE_TIME 0        /* correct effects of coarse observation time */
#define OUT_DBG_INFO 1            /* output debugs information to logger */
#define DETECT_OUTL  1            /* detect measurement outliers */
#define MAXTDIFF     5.0          /* max time difference to initial ekf */
#define MAXITR       30           /* max number of iteration for ekf */
#define ADRBFILTER   1            /* adaptive robust kalman filtering */

/* number and index of states */
#define ANP(opt)     (6)
#define ANC(opt)     (4+1)

#if CORR_CORSE_TIME
#define ANT(opt)     (1)
#else
#define ANT(opt)     (0)
#endif

#define ANX(opt)     (ANP(opt)+ANC(opt)+ANT(opt))
#define AIP(opt)     (0)
#define AIV(opt)     (3)
#define AIC(s,opt)   (ANP(opt)+(s))
#define AIT(opt)     (ANP(opt)+ANC(opt))

#define AOBS_I(val)  ((val>>16)&0xFF)
#define AOBS_F(val)  ((val>> 8)&0xFF)
#define AOBS_T(val)  ((val>> 4)&0xF)

#define AEXC_I(val)  (NFREQ*AOBS_I(val)+AOBS_F(val))

/* measurement residual callback function type--------------------------------*/
typedef int meas_res_func_t(const aobsd_t *obs, int n, const nav_t *nav, int f,
                            int uCN0, int uLSQ, int uSig0,
                            const double *x, agpos_t *agp, double *v, double *H,
                            double *var, const int *exc, int *ind);

extern sol_t* getphonesol(gtime_t time);
extern sol_t* getnlpsol(gtime_t time);
extern sol_t* getrtksol(gtime_t time);

static void initxP(double *x, double *P, double xi, double var, int i, int nx)
{
    int j;
    for (x[i]=xi,j=0;j<nx;j++) {
        P[i+j*nx]=P[j+i*nx]=i==j?var:0.0;
    }
}
/* solution to covariance ----------------------------------------------------*/
static void soltocov(const aslsq_t *sol, double *P, int type)
{
    double s=type?SQR(sol->sig0_p0):SQR(sol->sig0_v0);

    P[0]     =type?s*sol->qr[0]:s*sol->qv[0];
    P[4]     =type?s*sol->qr[1]:s*sol->qv[1];
    P[8]     =type?s*sol->qr[2]:s*sol->qv[2];
    P[1]=P[3]=type?s*sol->qr[3]:s*sol->qv[3];
    P[5]=P[7]=type?s*sol->qr[4]:s*sol->qv[4];
    P[2]=P[6]=type?s*sol->qr[5]:s*sol->qv[5];
}
/* determine whether initial ekf filter states--------------------------------*/
static int determine_init_state(const asflt_t *flt, const gtime_t cur_time, const prcopt_t *opt)
{
    double var;
    int i;

    trace(3,"determine_init_state: nx=%d\n",flt->nx);

    /* check whether at first epoch */
    if (flt->x&&norm(flt->x,3)<=0.0) {
        trace(2,"first epoch to initial\n");
        return 1;
    }
    /* check time diff of between current and last filter */
    if (timediff(cur_time,flt->time)>MAXTDIFF) {
        trace(2,"re-init due to large time age\n");
        return 1;
    }
    /* check variance of estimated position */
    for (var=i=0;i<3;i++) {
        var+=flt->P[i+i*flt->nx];
    }
    if (var/3.0>10*VAR_POS) {
        trace(2,"re-initial due to large var: %.3f\n",var);
        return 2;
    }
    return 0;
}
/* check initial ekf filter states conditions--------------------------------*/
static int chk_init_cond(const aslsq_t *lsq, const prcopt_t *opt, int reset_flag)
{
    double min_pos_dop=10.0;
    double min_vel_dop=10.0;
    double min_sig_pos=reset_flag?60.0:300.0;
    double min_sig_vel=10.0;
    double min_sig0_pos=reset_flag?60.0:100.0;
    double min_sig0_vel=reset_flag?5.0:10.0;
    double Qxyz[9]={0},Qpe[9],Qve[9],pos[3];

    trace(3,"chk_init_cond: sow=%6.3lf\n",time2gpst(lsq->time,NULL));

    ecef2pos(lsq->rr,pos);

    /* position/velocity cov. */
    soltocov(lsq,Qxyz,1);
    covenu(pos,Qxyz,Qpe);

    soltocov(lsq,Qxyz,0);
    covenu(pos,Qxyz,Qve);

    trace(3,"position Qenu=\n"); tracemat(3,Qpe,3,3,12,5);
    trace(3,"velocity Qenu=\n"); tracemat(3,Qve,3,3,12,5);

    /* check lsq solution status */
    if (!lsq->psta||!lsq->vsta) {
        trace(2,"invalid position or velocity to initial\n");
        return 0;
    }
    /* check lsq position solution variance */
    if (sqrt(Qpe[0]+Qpe[4])>min_sig_pos||sqrt((Qpe[0]+Qpe[4])/SQR(lsq->sig0_p0))>(reset_flag?30.0:100.0)) {
        trace(2,"large variance of position to initial\n");
        return 0;
    }
    /* check lsq velocity solution variance */
    if (sqrt(Qve[0]+Qve[4])>min_sig_vel||sqrt((Qve[0]+Qve[4])/SQR(lsq->sig0_v0))>10.0) {
        trace(2,"large variance of velocity to initial\n");
        return 0;
    }
    /* check lsq position DOPS */
    if (lsq->dop_pos[0]<=0.0||lsq->dop_pos[2]>min_pos_dop) {
        trace(2,"check position dops fail\n");
        return 0;
    }
    /* check lsq velocity DOPS */
    if (lsq->dop_vel[0]<=0.0||lsq->dop_vel[2]>min_vel_dop) {
        trace(2,"check velocity dops fail\n");
        return 0;
    }
    /* check lsq solution sig0 */
    if (lsq->sig0_p1>min_sig0_pos||lsq->sig0_v1>min_sig0_vel) {
        trace(2,"check position/velocity sig0 fail\n");
        return 0;
    }
#if !CORR_COARSE_TIME
    /* check correction of priori observation coarse-time */
    if (fabs(lsq->dt)>3000.0*1E-3) {
        trace(2,"check correction of priori observation coarse-time fail\n");
        return 0;
    }
#endif
    trace(3,"check initial conditions ok\n");
    return 1;
}
/* check whether rcv clock has been initialized -----------------------------*/
static void chk_ISB_init(asflt_t *flt, const aslsq_t *lsq)
{
    int i;

    for (i=1;i<4;i++) {
        if (!lsq->psta||!lsq->dtr[i]) continue;
        if (flt->x[AIC(i,opt)]) continue;

        initxP(flt->x,flt->P,lsq->dtr[i]-lsq->dtr[0],VAR_CLK,AIC(i,opt),flt->nx);
        trace(3,"rcv ISB init: i=%d\n",i);
    }
}
/* initial correction of priori observation coarse-time----------------------*/
static void init_corr_obst(double *x, double *P, int nx, const aslsq_t *lsq,
                           const prcopt_t *opt)
{
    if (CORR_COARSE_TIME) {
        initxP(x,P,lsq->psta?fabs(lsq->dt)>5.0?0.0:lsq->dt:0.0,VAR_CDT,AIT(opt),nx);
    }
}
/* initial filter estimated states-------------------------------------------*/
static int init_flt_state(agpos_t *agp, const aobsd_t *obs, int n, int init_flag, int reset_flg,
                          double *x0, double *P0)
{
    const aslsq_t *lsq=&agp->spp.lsq;
    int i,nx=agp->spp.flt.nx;
    double *x,*P;
    asflt_t *flt=&agp->spp.flt;
    const prcopt_t *opt=&agp->opt;

    trace(3,"initial estimated states: n=%d\n",n);

    /**
     * reset_flg: 1-position, 2-velocity, 3-position/velocity
     */
    /* check whether ISB has been initialized */
    chk_ISB_init(flt,&agp->spp.lsq);

    /* determine whether initial ekf filter states */
    if (!reset_flg&&!determine_init_state(flt,obs[0].time,&agp->opt)) {
        return 1;
    }
    /* check initial conditions */
    if (!chk_init_cond(&agp->spp.lsq,&agp->opt,reset_flg)) {
        trace(2,"check initial conditions fail\n");
        return 0;
    }
    x=zeros(1,nx); P=zeros(nx,nx);
    if (x0) matcpy(x,x0, 1,nx);
    if (P0) matcpy(P,P0,nx,nx);

    /* initialize position */
    if (init_flag||(reset_flg==1||reset_flg==3)) {
        for (i=AIP(opt);i<AIP(opt)+3;i++) {
            initxP(x,P,lsq->rr[i],reset_flg?VAR_POS:lsq->qr[i]*SQR(lsq->sig0_p0),i,nx);
        }
    }
    /* initialize velocity */
    if (init_flag||(reset_flg==2||reset_flg==3)){
        for (i=AIV(opt);i<AIV(opt)+3;i++) {
            initxP(x,P,lsq->rr[i],reset_flg?VAR_VEL:lsq->qv[i-AIV(opt)]*SQR(lsq->sig0_v0),i,nx);
        }
    }
    /* initialize clock bias/rate */
    initxP(x,P,lsq->dtr[0],VAR_CLK,AIC(0,opt),nx);
    initxP(x,P,lsq->ddrt,VAR_CLK,AIC(4,opt),nx);

    for (i=1;i<4;i++) {
        initxP(x,P,lsq->dtr[i]?lsq->dtr[i]-lsq->dtr[0]:0.0,VAR_CLK,AIC(i,opt),nx);
    }
    /* initial correction of priori obs coarse-time */
    init_corr_obst(x,P,nx,lsq,opt);

    if (reset_flg) {
        matcpy(x0,x, 1,nx);
        matcpy(P0,P,nx,nx);
    }
    if (init_flag) {
        flt->time=lsq->time;
        flt->tt=timediff(flt->time,obs[0].time);

        matcpy(flt->x,x, 1,nx);
        matcpy(flt->P,P,nx,nx);
    }
    trace(3,"initialize estimated states ok\n");

    trace(3,"x0=\n");
    tracemat(3,x,1,nx,12,5);

    trace(3,"P0=\n");
    tracemat(3,P,nx,nx,12,5);

    free(x); free(P);
    return 1;
}
/* estimate process noise ---------------------------------------------------*/
static void est_prc_noise(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav,
                          const double dt, double *cpsd)
{
    static gtime_t ptime;
    static double ppsd[5];
    double pos[3],acc[3],ae[3],Qxyz[9],Qv[9];
    const aslsq_t *lsq=&agp->spp.lsq;
    int i;

    trace(3,"est_prc_noise: dt=%6.3lf\n",dt);

    if (!ptime.time||timediff(obs[0].time,ptime)>3.0) {
        matcpy(ppsd,agp->opt.psd,1,5);
        matcpy(cpsd,agp->opt.psd,1,5);
        ptime=obs[0].time;
        return;
    }
    matcpy(cpsd,ppsd,1,5);

    if (!lsq->vsta||lsq->sig0_v1>1.0||fabs(dt)<1E-2) return;
    soltocov(lsq,Qxyz,0);

    ecef2pos(agp->sol.rr,pos);
    covenu(pos,Qxyz,Qv);
    if (sqrt(Qv[0]+Qv[4])>0.5) return;

    for (i=0;i<3;i++) {
        acc[i]=(lsq->rr[i+3]-agp->spp.flt.x[i+3])/dt;
    }
    ecef2enu(pos,acc,ae);

    cpsd[0]=SQR(1.5*ae[0])*dt;
    cpsd[1]=SQR(1.5*ae[1])*dt;

    ptime=obs[0].time;
    matcpy(ppsd,cpsd,1,5);
}
/* update ekf estimated states-----------------------------------------------*/
static int timeupdflt(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav,
                      double *x, double *P)
{
    asflt_t *flt=&agp->spp.flt;
    double *F,*FP,*Qk,pos[3],Q[9]={0},Qv[9];
    double dt=timediff(obs[0].time,flt->time),psd[5];
    int i,j,nx=flt->nx;

    trace(3,"update ekf estimated states: dt=%.2lf\n",dt);

    if (dt<-0.5) {
        trace(2,"time difference should be >0.0\n");
        return 0;
    }
    trace(3,"xk0=\n");
    tracemat(3,flt->x,1,nx,8,3);

    trace(3,"Pk0=\n");
    tracemat(3,flt->P,nx,nx,8,2);

    matcpy(P,flt->P,nx,nx);
    matcpy(x,flt->x, 1,nx);
    flt->tt=dt;

    /* estimate process noise */
    est_prc_noise(agp,obs,n,nav,dt,psd);

    trace(3,"psd=%.3lf %.3lf %.3lf\n",psd[0],psd[1],psd[2]);

    /* update receiver position/velocity */
    for (i=0;i<3;i++) {
        x[i+AIV(opt)]=flt->x[i+AIV(opt)];
        x[i+AIP(opt)]=flt->x[i+AIP(opt)]+flt->x[i+AIV(opt)]*dt;
    }
    /* update receiver clock bias */
    x[AIC(0,opt)]=flt->x[AIC(0,opt)]+flt->x[AIC(4,opt)]*dt;

    /* receiver clock bias rate */
    x[AIC(4,opt)]=flt->x[AIC(4,opt)];

    /* inter-system clock bias */
    for (i=1;i<4;i++) x[AIC(i,opt)]=flt->x[AIC(i,opt)];

    F=eye(nx); FP=mat(nx,nx);
    Qk=zeros(nx,nx);

    /* state transition of position */
    for (i=AIP(opt);i<AIP(opt)+3;i++) {
        for (j=AIV(opt);j<AIV(opt)+3;j++) if (i-AIP(opt)==j-AIV(opt)) F[i+j*nx]=dt;
    }
    /* clock error rate */
    F[AIC(4,opt)+AIC(4,opt)*nx]=exp(-dt/CORR_TIME_CLKD);
    F[AIC(0,opt)+AIC(4,opt)*nx]=CORR_TIME_CLKD-CORR_TIME_CLKD*exp(-dt/CORR_TIME_CLKD);

    /* receiver inter-system clock bias */
    for (i=1;i<4;i++) {
        if (!x[AIC(i,opt)]) continue;
        F[AIC(i,opt)+AIC(i,opt)*nx]=exp(-dt/CORR_TIME_ISB);
    }
    trace(3,"Fk_\n");
    tracemat(3,F,nx,nx,8,2);

    /* process noise added to only acceleration */
    ecef2pos(x,pos);
    Q[0]=psd[0];
    Q[4]=psd[1];
    Q[8]=psd[2];
    covecef(pos,Q,Qv);

    /* system noise matrix for position/velocity */
    for (i=AIP(opt);i<AIP(opt)+3;i++) {
        for (j=AIP(opt);j<AIP(opt)+3;j++) {
            Qk[i+j*nx]=1.0/3.0*Qv[i-AIP(opt)+(j-AIP(opt))*3]*pow(dt,3.0);
        }
        for (j=AIV(opt);j<AIV(opt)+3;j++) {
            Qk[i+j*nx]=1.0/2.0*Qv[i-AIP(opt)+(j-AIV(opt))*3]*SQR(dt);
        }
    }
    for (i=AIV(opt);i<AIV(opt)+3;i++) {
        for (j=AIP(opt);j<AIP(opt)+3;j++) {
            Qk[i+j*nx]=1.0/2.0*Qv[i-AIV(opt)+(j-AIP(opt))*3]*SQR(dt);
        }
        for (j=AIV(opt);j<AIV(opt)+3;j++) {
            Qk[i+j*nx]=Qv[i-AIV(opt)+(j-AIV(opt))*3]*dt;
        }
    }
    /* system noise matrix for clock rate */
    i=AIC(0,opt);
    j=AIC(4,opt);
    Qk[i+i*nx]=1.0/3.0*psd[3]*pow(dt,3.0)+psd[4]*dt;
    Qk[i+j*nx]=0.5*psd[3]*SQR(dt);
    Qk[j+i*nx]=0.5*psd[3]*SQR(dt);
    Qk[j+j*nx]=psd[3]*dt;

    /* system noise matrix for ISB */
    for (i=1;i<4;i++) {
        if (!x[AIC(i,opt)]) continue;
        Qk[AIC(i,opt)+AIC(i,opt)*nx]=SQR(CLIGHT*agp->opt.sclkstab*dt);
    }
    trace(3,"Qk_=\n");
    tracemat(3,Qk,nx,nx,8,2);

    /* P=F*P*F+Q */
    matmul("NN",nx,nx,nx,1.0,F,flt->P,0.0,FP);
    matmul("NT",nx,nx,nx,1.0,FP,F,0.0,P);

    for (i=0;i<nx;i++) {
        for (j=0;j<nx;j++) P[i+j*nx]+=Qk[i+j*nx];
    }
    /* initial correction of priori obs coarse-time */
    init_corr_obst(x,P,nx,&agp->spp.lsq,&agp->opt);

    trace(3,"xk_=\n");
    tracemat(3,x,1,nx,8,3);

    trace(3,"Pk_=\n");
    tracemat(3,P,nx,nx,8,2);

    matcpy(flt->Pp,P,nx,nx);
    matcpy(flt->xp,x, 1,nx);

    free(FP); free(F); free(Qk);
    return 1;
}
/* get tgd parameter (m) -----------------------------------------------------*/
static double gettgd(int sat, int code, const nav_t *nav, double *tgd1, double *tgd2)
{
    int i,sys=satsys(sat,NULL);

    for (i=0;i<nav->n;i++) {
        if (nav->eph[i].sat!=sat) continue;
        if (sys==SYS_CMP||sys==SYS_GAL) {
            if (tgd1) *tgd1=CLIGHT*nav->eph[i].tgd[0];
            if (tgd2) *tgd2=CLIGHT*nav->eph[i].tgd[1];
        }
        return CLIGHT*nav->eph[i].tgd[0];
    }
    return 0.0;
}
/* psendorange with code bias correction -------------------------------------*/
static double prange(const aobsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var, int f)
{
    double PC,P,P1_C1,P2_C2,tgd1=0.0,tgd2=0.0,tgd;
    double beta;
    int i,if_nav=0,sys;

    if (var) *var=0.0;

    if (!(sys=satsys(obs->sat,NULL))) return 0.0;
    if (!(PC=P=obs->P[f])) return 0.0;

    P1_C1=nav->cbias[obs->sat-1][1];
    P2_C2=nav->cbias[obs->sat-1][2];
    tgd=gettgd(obs->sat,obs->code[f],nav,&tgd1,&tgd2);

    beta=SQR(sat2freq(obs->sat,obs->code[0],nav))/SQR(sat2freq(obs->sat,obs->code[f],nav));

    /* GPS/QZSS/GLONASS */
    if (sys==SYS_GPS||sys==SYS_QZS||sys==SYS_GLO) {
        if (obs->code[f]==CODE_L1C) {
            P+=P1_C1; /* C1->P1 */
            PC=P-beta*tgd;
        }
        else if (obs->code[f]==CODE_L2C||obs->code[f]==CODE_L2X||obs->code[f]==CODE_L2L||
                 obs->code[f]==CODE_L2S) {
            P+=P2_C2; /* C2->P2 */
            PC=P-beta*tgd;
        }
        else if (obs->code[f]==CODE_L1P) {
            PC=P-tgd;
        }
        else if (obs->code[f]==CODE_L2P) {
            PC=P-beta*tgd;
        }
    }
    /* BDS */
    else if (sys==SYS_CMP) {
        if      (obs->code[f]==CODE_L1I||obs->code[f]==CODE_L2I) PC=P-tgd1;
        else if (obs->code[f]==CODE_L7I) PC=P-tgd2;
    }
    /* GALILEO */
    else if (sys==SYS_GAL) {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==obs->sat) {
                if_nav=nav->eph[i].IFnav;
                break;
            }
        }
        if (if_nav) {
            if (obs->code[f]==CODE_L5A) {
                PC=P-SQR(FREQ1/FREQ5)*tgd1;
            }
        }
        else {
            if (obs->code[f]==CODE_L1C||obs->code[f]==CODE_L1B) {
                PC=P-tgd2;
            }
            else if (obs->code[f]==CODE_L5B) {
                PC=P-SQR(FREQ1/FREQ5)*tgd2;
            }
        }
    }
    /* other system */
    else {
        PC=P-tgd;
    }
    if (var) *var=VAR_CBIAS;
    return PC;
}
/* code residuals -----------------------------------------------------------*/
static int prng_res(const aobsd_t *obs, int n, const nav_t *nav, int uf, int uCN0,
                    int uLSQ, int uSig0,
                    const double *x, agpos_t *agp, double *v, double *H, double *var, const int *exc,
                    int *ind)
{
    double r,dion,dtrp,vmeas,vion,vtrp,rr[3],pos[3],e[3],P,lam_L1,azel[2],sig,dtr;
    double sig0=agp->spp.lsq.sig0_p0;
    double sat_pr,corr_dt,freq;
    int i,k,f,sys,sat,nx=agp->spp.flt.nx,nv,nf=NFREQ;
    assat_t *ssat=agp->ssat;
    const prcopt_t *opt=&agp->opt;

    trace(3,"prng_res: n=%d sig0=%.3lf\n",n,agp->spp.lsq.sig0_p0);

    for (i=0;i<3;i++) rr[i]=x[i]; dtr=x[AIC(0,opt)];
    ecef2pos(rr,pos);
    agp->spp.flt.ns=0;

    if (CORR_COARSE_TIME) {
        corr_dt=x[AIT(opt)];
    }
    else {
        corr_dt=agp->spp.lsq.dt;
    }
    for (nv=f=0;f<nf;f++) {
        for (i=0;i<n&&i<MAXOBS;i++) {
            if (uf>=0&&f!=uf) continue;

            sat=obs[i].sat; sys=satsys(sat,NULL);
            ssat[sat-1].vs_flt[f]=0;

            /* geometric distance/az/el angle */
            if ((r=geodist(ssat[sat-1].rs,x,e))<=0.0||satazel(pos,e,azel)<opt->elmincode) {
                continue;
            }
            /* code bias correction */
            if ((P=prange(obs+i,nav,opt,&vmeas,f))<=RE_WGS84) continue;
            if (H) {
                for (k=0;k<nx;k++) H[k+nv*nx]=0.0;
            }
            /* ionospheric corrections */
            ionocorr(obs[i].time,nav,obs[i].sat,pos,azel,opt->ionoopt,&dion,&vion);

            if (!(freq=sat2freq(obs[i].sat,obs[i].code[f],nav))) continue;
            dion*=SQR(FREQ1/freq);
            vion*=SQR(FREQ1/freq);

            /* tropospheric corrections */
            tropcorr(obs[i].time,nav,pos,azel,opt->tropopt,&dtrp,&vtrp);

            /* pseudorange residual */
            v[nv]=P-(r+dtr-CLIGHT*ssat[sat-1].dts[0]+dion+dtrp);
            if (H) {
                for (k=AIP(opt);k<(AIP(opt)+3);k++) H[k+nx*nv]=-e[k];
            }
            /* rcv clock bias offset correction */
            if      (sys==SYS_GLO) {v[nv]-=x[AIC(1,opt)]; if (H) H[AIC(1,opt)+nv*nx]=1.0;}
            else if (sys==SYS_GAL) {v[nv]-=x[AIC(2,opt)]; if (H) H[AIC(2,opt)+nv*nx]=1.0;}
            else if (sys==SYS_CMP) {v[nv]-=x[AIC(3,opt)]; if (H) H[AIC(3,opt)+nv*nx]=1.0;}

            matmul("NT",1,1,3,1.0,e,ssat[sat-1].rs+3,0.0,&sat_pr);
            sat_pr-=ssat[sat-1].dts[1]*CLIGHT;

            /* correct effects of coarse observation time */
            v[nv]-=sat_pr*corr_dt;
            if (H) {
                if (CORR_COARSE_TIME) {
                    H[AIT(opt)+nv*nx]=sat_pr;
                }
            }
            ssat[sat-1].resp[f]=v[nv];
            matcpy(ssat[sat-1].azel,azel,1,2);

            /* check observation cut-off CN0 */
            if (uCN0&&obs[i].SNR[f]*0.25<agp->spp.lsq.min_cn0[0]) continue;

            /* excluded satellite? */
            if (uLSQ&&!ssat[sat-1].vs_lsq[f]) continue;
            if ((exc&&exc[i*NFREQ+f])) continue;

            /* error variance */
            sig=codevarerr(opt,obs[i].SNR[f]*0.25,azel[1],obs[i].sat)+ssat[sat-1].var_s+vmeas+vion+vtrp;
            if (uSig0) sig*=SQR(sig0);

            if (var) var[nv]=sig;
            if (ind) {
                ind[nv]=((i<<16)|(f<<8)|(1<<4));
            }
            nv++;
            ssat[sat-1].vs_flt[f]=1;
            agp->spp.flt.ns++;

            char prn[8];
            satno2id(obs[i].sat,prn);

            trace(3,"%4s sat=%3d azel=%5.1f %5.1f res=%10.3f cdt=%10.3f sig=%8.3f cn0=%8.3lf f=%d stat=%4d\n",
                  prn,obs[i].sat,azel[0]*R2D,azel[1]*R2D,
                  v[nv-1],sat_pr*corr_dt,sqrt(sig),obs[i].SNR[f]*0.25,
                  f,obs[i].unc.stat[f]);
        }
    }
    return nv;
}
/* doppler residuals----------------------------------------------------------*/
static int dopp_res(const aobsd_t *obs, int n, const nav_t *nav, int uf, int uCN0,
                    int uLSQ, int uSig0,
                    const double *x, agpos_t *agp, double *v, double *H, double *var,
                    const int *exc, int *ind)
{
    double rate,pos[3],E[9],a[3],e[3],vs[3],cosel,y,sig;
    double corr_dt=0.0;
    const double *azel,*rs,*dts,*rr,*vr;
    char prn[8];
    int i,j,f=0,nv=0,sat,nx=agp->spp.flt.nx;
    const prcopt_t *opt=&agp->opt;
    assat_t *ssat=agp->ssat;

    trace(3,"dopp_res: n=%d sig0=%.3lf\n",n,agp->spp.lsq.sig0_v0);

    rr=x; vr=x+AIV(opt);
    ecef2pos(x,pos); xyz2enu(pos,E);
    agp->spp.flt.nd=0;

    if (CORR_COARSE_TIME) {
        corr_dt=x[AIT(opt)];
    }
    else {
        corr_dt=agp->spp.lsq.dt;
    }
    for (f=0;f<NFREQ;f++) {
        for (i=0;i<n&&i<MAXOBS;i++) {
            if (uf>=0&&f!=uf) continue;

            ssat[obs[i].sat-1].vd_flt[f]=0;
            sat=obs[i].sat;

            rs=ssat[sat-1].rs; dts=ssat[sat-1].dts;
            azel=ssat[sat-1].azel;

            /* pseudorange rate measurement data */
            if (!(y=obs[i].Pr[f])) continue;

            /* line-of-sight vector in ecef */
            cosel=cos(azel[1]);
            a[0]=sin(azel[0])*cosel;
            a[1]=cos(azel[0])*cosel;
            a[2]=sin(azel[1]);
            matmul("TN",3,1,3,1.0,E,a,0.0,e);

            /* satellite velocity relative to receiver in ecef */
            for (j=0;j<3;j++) vs[j]=rs[j+3]+rs[j+6]*corr_dt-x[j+AIV(opt)];

            /* range rate with earth rotation correction */
            rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4]*rr[0]+rs[1]*vr[0]-rs[3]*rr[1]-rs[0]*vr[1]);

            /* doppler residual */
            v[nv]=y-(rate+x[AIC(4,opt)]-CLIGHT*(dts[1]+dts[2]*corr_dt));

            /* jacobian of doppler wrt velocity */
            if (H) {
                for (j=0;j<3;j++) H[j+AIV(opt)+nv*nx]=-e[j];
                H[AIC(4,opt)+nv*nx]=1.0;
            }
            /* measurement variance */
            sig=doppvarerr(opt,obs[i].SNR[f]*0.25,azel[1],obs[i].sat);
            if (uSig0) sig*=SQR(agp->spp.lsq.sig0_v0);

            ssat[sat-1].resd[f]=v[nv];

            /* check elevation mask angle */
            if (norm(rs+3,3)<=0.0||azel[1]<opt->elmindopp*D2R) continue;

            /* check observation cut-off CN0 */
            if (uCN0&&obs[i].SNR[f]*0.25<agp->spp.lsq.min_cn0[1]) continue;

            /* exclude satellite? */
            if (exc&&exc[i*NFREQ+f]) continue;
            if (uLSQ&&!ssat[sat-1].vd_lsq[f]) continue;

            if (var) var[nv]=sig;
            if (ind) {
                ind[nv]=((i<<16)|(f<<8)|(2<<4));
            }
            satno2id(obs[i].sat,prn);
            trace(3,"%4s sat=%3d azel=%5.1f %5.1f res=%8.3f prr=%8.3lf sig=%8.3lf cn0=%8.3lf f=%d\n",
                  prn,obs[i].sat,azel[0]*R2D,azel[1]*R2D,v[nv],
                  y,sqrt(sig),obs[i].SNR[f]*0.25,f);

            ssat[sat-1].vd_flt[f]=1;
            agp->spp.flt.nd++;
            nv++;
        }
    }
    return nv;
}
/* update filter solution----------------------------------------------------*/
static int updfltsol(gtime_t obs_time, agpos_t *agp, const double *x, const double *P,
                     int nx, int stat)
{
    double pos[3],E[9],vn[3];
    int i;
    asflt_t *flt=&agp->spp.flt;
    aslsq_t *lsq=&agp->spp.lsq;
    sol_t *sol=&agp->sol;

    trace(3,"updfltsol: sow=%6.3lf stat=%d\n",time2gpst(obs_time,NULL),stat);

    ecef2pos(flt->xp,pos);
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,flt->xp+3,0.0,vn);

    trace(3,"vn0=%8.4lf %8.4lf %8.4lf\n",vn[0],vn[1],vn[2]);

    if (stat>=1) {
        matcpy(flt->P,P,nx,nx);
        matcpy(flt->x,x, 1,nx);
    }
    else if (stat<=0) {
        matcpy(flt->x,flt->xp, 1,nx);
        matcpy(flt->P,flt->Pp,nx,nx);
    }
    sol->ns=flt->ns;
    sol->nd=flt->nd;

    /* current filter time */
    if (stat==1) {
        flt->time=sol->time=timeadd(obs_time,-x[AIC(0,opt)]/CLIGHT);
    }
    else if (lsq->psta) {
        flt->time=sol->time=lsq->time;
    }
    else flt->time=sol->time=obs_time;

    /* receiver position */
    for (i=0;i<3;i++) {
        sol->rr[i]=x[i];
        sol->qr[i]=(float)flt->P[i+i*nx];
    }
    sol->qr[3]=(float)flt->P[1+0*nx];
    sol->qr[4]=(float)flt->P[2+1*nx];
    sol->qr[5]=(float)flt->P[2+0*nx];

    /* receiver clock/clock drift error */
    for (i=0;i<4;i++) sol->dtr[i]=flt->x[AIC(i,opt)];
    sol->ddrt=flt->x[AIC(4,opt)];

    if (stat>0) {
        sol->stat=SOLQ_ASPFT;
    }
    else {
        sol->stat=SOLQ_DR;
    }
    ecef2pos(flt->x,pos);
    xyz2enu(pos,E);
    matmul("NN",3,1,3,1.0,E,flt->x+3,0.0,vn);

    trace(3,"vnk=%8.4lf %8.4lf %8.4lf\n",vn[0],vn[1],vn[2]);

    trace(3,"updfltsol: rr=%10.4lf %10.4lf %10.4lf vel=%10.4lf %10.4lf %10.4lf\n",
          flt->x[0],flt->x[1],flt->x[2],
          flt->x[3],flt->x[4],flt->x[5],
          vn[0],vn[1],vn[2]);

    trace(3,"xk=\n");
    tracemat(3,flt->x,1,flt->nx,12,5);

    trace(3,"Pk=\n");
    tracemat(3,flt->P,flt->nx,flt->nx,12,5);
    return stat;
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
/* grubbs test----------------------------------------------------------------*/
static int grb_detoutl(const double *v, int nv, int *outl_ind, double alpha)
{
    double max_val,mean,std,tcrit,zcrit;
    int i,flag[MAXOBS]={0},j,max_i,ind_[MAXOBS],outl_n=0;

    if (alpha<=0.0) return 0;
    while (1) {
        for (mean=i=j=0;i<nv;i++) {
            if (flag[i]) continue;
            mean+=v[i]; ind_[j]=i;
            j++;
        }
        if (j<=2) break;
        mean/=j;

        for (max_val=fabs(v[ind_[0]]-mean),max_i=i=0;i<j;i++) {
            if (fabs(v[ind_[i]]-mean)>fabs(v[ind_[max_i]]-mean)) {
                max_val=v[ind_[i]];
                max_i=i;
            }
        }
        for (std=i=0;i<j;i++) std+=SQR(v[ind_[i]]-mean);
        std=sqrt(std/(j-1));

        tcrit=re_t(j-2,alpha/(2*j));
        zcrit=(j-1)/sqrt(j)*(sqrt(SQR(tcrit)/(j-2+SQR(tcrit))));
        if (fabs(max_val-mean)/std>zcrit) {

            /* detect outliers */
            if (outl_ind) outl_ind[outl_n]=ind_[max_i];
            flag[ind_[max_i]]=1;
            outl_n++;

            trace(2,"exclude sat (grubbs): index=%d\n",ind_[max_i]);
            continue;
        }
        else {
            /* no detect outliers */
            break;
        }
    }
    return outl_n;
}
/* detect method-------------------------------------------------------------*/
static int detoutl_method(const double *v, int nv, int *outl_ind, int *flag, double thres_box,
                          double thres_mad, double thres_grb)
{
    int i,outl_n=0;
    for (i=0;flag&&i<nv;i++) flag[i]=1;
    if (thres_box) {
        outl_n+=boxplot_detoutl(v,nv,outl_ind,thres_box);
    }
    if (thres_mad) {
        outl_n+=mad_detoutl(v,nv,outl_ind+outl_n,thres_mad);
    }
    if (thres_grb) {
        outl_n+=grb_detoutl(v,nv,outl_ind+outl_n,thres_grb);
    }
    for (i=0;flag&&i<outl_n;i++) flag[outl_ind[i]]=0;
    return outl_n;
}
/* ekf filter----------------------------------------------------------------*/
static int filt(double *x, double *P, double *dx, double *Ko, const double *H, const double *v,
                const double *var,
                int n, int m, int upd_flag)
{
    double *xp,*Pp,*F,*Q,*K,*I,*R,lam=1.0;
    int i,info;

    xp=mat(n,1); F=mat(n,m);
    K=mat(n,m); I=eye(n);
    Pp=mat(n,n); Q=mat(m,m); R=zeros(m,m);

    for (i=0;i<m;i++) R[i+i*m]=var[i];
    matcpy(xp,x,n,1);
    matcpy(Q,R,m,m);

    matmul("NN",n,m,n,lam,P,H,0.0,F); /* Q=H*P*H'+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    matcpy(R,Q,m,m);

    /* ekf filter */
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0, K); /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp); /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I); /* Pp=(I-K*H)*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);

        /* dx=K*v */
        if (dx) {
            matmul("NN",n,1,m,1.0,K,v,0.0,dx);
        }
        if (Ko) {
            matcpy(Ko,K,n,m);
        }
        /* update x/P */
        if (upd_flag) {
            matcpy(x,xp,1,n);
            matcpy(P,Pp,n,n);
        }
    }
    free(K); free(I); free(xp);
    free(R); free(F); free(Pp); free(Q);
    return info;
}
/* validation of ekf post-fit residuals--------------------------------------*/
static int valsol(int type, agpos_t *agp, const double *v, const int *ind, int nv, const double *var,
                  const double *H, const double *xp, const double *Pp,
                  const double *xe, const double *Pe, int nx,
                  const aobsd_t *obs, double *vr,
                  int *clk_jump_ind, int *reset_ind)
{
    double *g,*R,*M,*N;
    double thres_gain=4.0,thres_res=4.0;
    int i,flag=1,reset_pos=0,reset_vel=0;
    unsigned int iflg[6]={0};

    trace(3,"valsol: nv=%d\n",nv);

    g=mat(nx,1); R=mat(nv,nv);
    N=mat(nv,nv);
    M=mat(nv,nv);

    matmul33("TNN",H,Pp,H,nv,nx,nx,nv,R);
    for (i=0;i<nv;i++) R[i+i*nv]+=var[i],M[i+i*nv]=var[i];

    if (matinv(R,nv)) {
        free(g); free(R); free(N);
        free(M);
        return 0;
    }
    matmul33("NNN",M,R,M,nv,nv,nv,nv,N);

    for (i=0;i<nx;i++) g[i]=xe[i]-xp[i];
    trace(3,"g=\n");
    tracemat(3,g,1,nx,8,3);

    /* check ekf gain */
    for (*clk_jump_ind=-1,*reset_ind=0,i=0;i<nx;i++) {
        if ((vr[i]=fabs(g[i]/sqrt(Pp[i+i*nx])))>thres_gain) {
            trace(2,"large gain: i=%d r=%.3lf\n",i,vr[i]);

            if (i>=6) {
                if (*clk_jump_ind<0||vr[i]>vr[*clk_jump_ind]) {
                    *clk_jump_ind=i;
                }
            }
            else if (i<3) reset_pos++;
            else if (i<6) reset_vel++;
        }
        else {
            vr[i]=0.0;
        }
    }
    /* check post-fit residuals */
    for (i=0;i<nv;i++) {
        if (AOBS_T(ind[i])!=type+1) continue;
        vr[nx+i]=fabs(v[i])/sqrt(N[i+i*nv]);

        if (vr[nx+i]>thres_res) {
            trace(2,"large post-fit residual: sat=%3d f=%d type=%d cn0=%6.3lf r=%.3lf\n",
                  obs[AOBS_I(ind[i])].sat,
                  AOBS_F(ind[i]),
                  AOBS_T(ind[i]),
                  obs[AOBS_I(ind[i])].SNR[AOBS_F(ind[i])]*0.25,
                  vr[nx+i]);
            flag=0;
        }
    }
    free(g); free(R); free(N);
    free(M);

    if (*clk_jump_ind>=6) return -1;

    if (reset_pos||reset_vel) {
        if (reset_pos>=2) *reset_ind=1;
        if (reset_vel>=2) *reset_ind=2;
        if (*reset_ind) return -2;
    }
    return flag;
}
/* exclude measurement outliers----------------------------------------------*/
static int excludeoutl(const int *ind, const double *vr, const double *var, int nv,
                       const aobsd_t *obs, int *exc)
{
    int i,j;

    for (j=-1,i=0;i<nv;i++) {
        if (vr[i]==0.0) continue;
        if (j<0||fabs(vr[i])>fabs(vr[j])) {
            j=i;
        }
    }
    if (j>=0) {
        exc[AEXC_I(ind[j])]=1;
        trace(2,"exclude: sat=%3d f=%d type=%d\n",
              obs[AOBS_I(ind[j])].sat,
              AOBS_F(ind[j]),
              AOBS_T(ind[j]));
    }
    return j;
}
/* filter model validation---------------------------------------------------*/
static int val_fltmdl(int type, int reset_cnt, const double *v, const int *ind, const double *H,
                      const double *var, const aobsd_t *obs, double *x, double *P,
                      int nv, int nx, int *exc)
{
    double *Rv,*vv,thres=re_norm(0.999);
    double s,thres_bad=8.0,thres_reset=4.0;
    int i,ne=0;

    trace(3,"val_fltmdl: type=%d\n",type);

    Rv=mat(nv,nv);
    vv=mat(nv,1);

    /* test of statistics */
    matmul33("TNN",H,P,H,nv,nx,nx,nv,Rv);
    for (i=0;i<nv;i++) Rv[i+i*nv]+=var[i];

    if (matinv(Rv,nv)) {
        free(Rv);
        free(vv);
        return 0;
    }
    matmul33("NNT",v,Rv,v,1,nv,nv,1,&s);
    trace(3,"s=%.3lf chispr=%.3lf\n",s,chisqr[nv-1]);

    if (s/chisqr[nv-1]>thres_bad) {
        free(Rv); free(vv);
        return -1;
    }
    if (reset_cnt) {
        if (s/chisqr[nv-1]>thres_reset) {
            free(Rv); free(vv);
            return -2;
        }
    }
    matmul("NN",nv,1,nv,1.0,Rv,v,0.0,vv);

    for (i=0;i<nv;i++) {
        if (fabs(vv[i]/sqrt(Rv[i+i*nv]))<thres) continue;
        if (exc[AEXC_I(ind[i])]) continue;
        exc[AEXC_I(ind[i])]=1;
        ne++;

        trace(2,"exclude: sat=%3d f=%d type=%d\n",
              obs[AOBS_I(ind[i])].sat,
              AOBS_F(ind[i]),
              AOBS_T(ind[i]));
    }
    free(Rv); free(vv);
    return nv-ne;
}
/* index navi system (m=0:gps/sbs/qzss/sbas,1:glo,2:gal,3:bds) -----------------*/
static int ind_sys(int sys)
{
    switch (sys) {
        case SYS_GPS: return 0;
        case SYS_SBS: return 0;
        case SYS_GLO: return 1;
        case SYS_GAL: return 2;
        case SYS_CMP: return 3;
        case SYS_QZS: return 0;
        default:
            return -1;
    }
}
/* reset rcv clock-----------------------------------------------------------*/
static void reset_clka(double *xp, double *Pp, int nx, const double clk_est, int i,
                       const prcopt_t *opt)
{
    trace(3,"reset_clk: i=%d\n",i);

    initxP(xp,Pp,xp[AIC(i,opt)]+clk_est,VAR_CLK,AIC(i,opt),nx);

    trace(3,"xp=\n"); tracemat(3,xp, 1,nx,12,5);
    trace(3,"Pp=\n"); tracemat(3,Pp,nx,nx,12,5);
}
/* reset ISB------------------------------------------------------------------*/
static void reset_ISB(double *xp, double *Pp, int nx, const double ISB_est, int i,
                      const prcopt_t *opt)
{
    trace(3,"reset_ISB: i=%d\n",i);

    initxP(xp,Pp,xp[AIC(i,opt)]+ISB_est,Pp[AIC(i,opt)+AIC(i,opt)*nx]*SQR(3.0),AIC(i,opt),nx);

    trace(3,"xp=\n"); tracemat(3,xp, 1,nx,12,5);
    trace(3,"Pp=\n"); tracemat(3,Pp,nx,nx,12,5);
}
/* reset rcv clock drift-----------------------------------------------------*/
static void reset_clka_dri(double *xp, double *Pp, int nx, const double clk_dri_est,
                           const prcopt_t *opt)
{
    trace(3,"reset_clk_dri:\n");

    initxP(xp,Pp,xp[AIC(1,op1)]+clk_dri_est,VAR_CLK,AIC(4,opt),nx);

    trace(3,"xp=\n"); tracemat(3,xp, 1,nx,12,5);
    trace(3,"Pp=\n"); tracemat(3,Pp,nx,nx,12,5);
}
/* process clock jump---------------------------------------------------------*/
static void init_rcv_clk(int clk_jump, double *xp, double *Pp, int nx, const double *xe,
                         const prcopt_t *opt)
{
    initxP(xp,Pp,xe[clk_jump],VAR_CLK,clk_jump,nx);

    trace(3,"xp=\n"); tracemat(3,xp, 1,nx,12,5);
    trace(3,"Pp=\n"); tracemat(3,Pp,nx,nx,12,5);
}
/* detect measurement outliers------------------------------------------------*/
static int detoutl_meas(const double *v, const int *ind, int nv, int *exc, const aobsd_t *obs,
                        double thres_box,
                        double thres_mad,
                        double thres_grb)
{
    int i,*outl=imat(3*nv,1),n=0,nk,*id=imat(nv,1);
    double *vk=mat(nv,1);

    for (nk=i=0;i<nv;i++) {
        if (exc[AEXC_I(ind[i])]) continue;
        vk[nk]=v[i];
        id[nk++]=ind[i];
    }
    n=detoutl_method(vk,nk,outl,NULL,thres_box,thres_mad,thres_grb);

    for (i=0;i<n;i++) {
        if (exc[AEXC_I(id[outl[i]])]) continue;
        exc[AEXC_I(id[outl[i]])]=1;

        trace(2,"exclude: sat=%3d f=%d type=%d cn0=%6.3lf\n",
              obs[AOBS_I(id[outl[i]])].sat,
              AOBS_F(id[outl[i]]),
              AOBS_T(id[outl[i]]),
              obs[AOBS_I(id[outl[i]])].SNR[AOBS_F(id[outl[i]])]*0.25);
    }
    free(outl);
    free(vk); free(id);
    return n;
}
/* detect rcv clock jump-----------------------------------------------------*/
static int detect_clk_jump(meas_res_func_t *res_func, const aobsd_t *obs, int n,
                           const nav_t *nav, agpos_t *agp,
                           double *xp, double *Pp)
{
    const double thres=re_norm(0.999);
    const double min_cn0=MAX(7.0,agp->spp.lsq.min_cn0[0]),min_el=10.0;
    double *v,*vi,med,clk_est;
    int i,m,mk,nx=agp->spp.flt.nx,nv=n*NFREQ*n+3,*ind,*flag;
    int *outl,noutl,jump_flg=0;

    trace(3,"detect_clk_jump_prc: n=%d\n",n);

    v=mat(nv,1); vi=mat(nv,1);
    ind=imat(nv,1); outl=imat(nv,1);
    flag=imat(nv,1);

    /* update measurement residuals */
    if ((nv=res_func(obs,n,nav,-1,1,1,1,xp,agp,v,NULL,NULL,NULL,ind))<=0) {
        trace(2,"no valid meas\n");

        free(v); free(vi); free(ind);
        free(outl); free(flag);
        return 0;
    }
    /* i=0:gps/sbs/qzs,1:glo,2:gal,3:bds */
    for (m=i=0;i<nv;i++) {
        if (obs[AOBS_I(ind[i])].SNR[AOBS_F(ind[i])]*0.25<min_cn0) continue;
        if (satsys(obs[AOBS_I(ind[i])].sat,NULL)!=SYS_GPS) continue;
        if (agp->ssat[obs[AOBS_I(ind[i])].sat-1].azel[1]<min_el*D2R) continue;
        vi[m++]=v[i];
    }
    /* detect meas outliers */
    noutl=detoutl_method(vi,m,outl,flag,1.5,2.5,0.001);
    for (mk=i=0;i<m;i++) {
        if (!flag[i]) continue;
        vi[mk++]=vi[i];
    }
    /* estimate rcv clock bias */
    if (mk<2||mad(vi,mk)>60.0) {
        free(outl); free(flag);
        free(v); free(vi); free(ind);
        return 0;
    }
    med=median(vi,mk);

    /* detect rcv clock jump */
    if (fabs(med)/sqrt(Pp[AIC(0,opt)+AIC(0,opt)*agp->spp.flt.nx])>thres) {
        reset_clka(xp,Pp,nx,med,0,&agp->opt);
        jump_flg++;
    }
    free(outl); free(flag);
    free(v); free(vi); free(ind);
    return jump_flg;
}
/* detect rcv clock drift jump-----------------------------------------------*/
static int detect_clk_drift_jump(meas_res_func_t *res_func, const aobsd_t *obs, int n,
                                 const nav_t *nav, agpos_t *agp,
                                 double *xp, double *Pp)
{
    const double thres=re_norm(0.999);
    double vm[MAXSAT],med,clk_est;
    static int flag[MAXSAT],outl[MAXSAT],noutl;
    int i,m,mk,nx=agp->spp.flt.nx,nv;

    trace(3,"detect_clk_drift_jump_prc: n=%d\n",n);

    /* measurement residuals */
    if ((nv=res_func(obs,n,nav,-1,1,1,1,xp,agp,vm,NULL,NULL,NULL,NULL))<=0) {
        return 0;
    }
    /* detect measurement outliers */
    noutl=detoutl_method(vm,nv,outl,flag,1.5,2.5,0.001);
    for (m=i=0;i<nv;i++) {
        if (!flag[i]) continue;
        vm[m++]=vm[i];
    }
    /* estimate rcv clock drift bias */
    if (m<2||mad(vm,m)>5.0) return 0;
    med=median(vm,m);

    /* detect rcv clock drift jump */
    if (fabs(med)/sqrt(Pp[AIC(4,opt)+AIC(4,opt)*agp->spp.flt.nx])>thres) {
        reset_clka_dri(xp,Pp,nx,med,&agp->opt);
        return 1;
    }
    return 0;
}
/* detect ISB jump-----------------------------------------------------------*/
static int detect_ISB_jump(meas_res_func_t *res_func, const aobsd_t *obs, int n,
                           const nav_t *nav, agpos_t *agp,
                           double *xp, double *Pp)
{
    const double thres=re_norm(0.999);
    const double min_cn0=20.0,min_el=15.0;
    double *v,*vi,med,clk_est;
    int i,j,m,mk,nx=agp->spp.flt.nx,nv=n*NFREQ*n+3,*ind,*flag;
    int *outl,noutl,jump_flg=0;

    trace(3,"detect_ISB_jump: n=%d\n",n);

    v=mat(nv,1); vi=mat(nv,1);
    ind=imat(nv,1); outl=imat(nv,1);
    flag=imat(nv,1);

    /* update measurement residuals */
    if ((nv=res_func(obs,n,nav,-1,1,1,1,xp,agp,v,NULL,NULL,NULL,ind))<=0) {
        trace(2,"no valid meas\n");

        free(v); free(vi); free(ind);
        free(outl); free(flag);
        return 0;
    }
    /* i=1:glo,2:gal,3:bds */
    for (i=1;i<4;i++) {
        for (m=j=0;j<nv;j++) {
            if (obs[AOBS_I(ind[i])].SNR[AOBS_F(ind[i])]*0.25<min_cn0) continue;
            if (ind_sys(satsys(obs[AOBS_I(ind[j])].sat,NULL))!=i) continue;
            if (agp->ssat[obs[AOBS_I(ind[j])].sat-1].azel[1]<min_el*D2R) continue;
            vi[m++]=v[j];
        }
        /* detect meas outliers */
        noutl=detoutl_method(vi,m,outl,flag,1.5,2.5,0.001);
        for (mk=j=0;j<m;j++) {
            if (!flag[j]) continue;
            vi[mk++]=vi[j];
        }
        /* estimate rcv clock bias */
        if (mk<2||mad(vi,mk)>10.0) continue;
        med=median(vi,mk);

        /* detect rcv clock jump */
        if (fabs(med)/sqrt(Pp[AIC(i,opt)+AIC(i,opt)*agp->spp.flt.nx])>thres) {
            reset_ISB(xp,Pp,nx,med,i,&agp->opt);
            jump_flg++;
        }
    }
    free(outl); free(flag);
    free(v); free(vi); free(ind);
    return jump_flg;
}
/* filter process------------------------------------------------------------*/
static int updfilt_prc(int type, meas_res_func_t *res_func, const aobsd_t *obs, int n,
                       const nav_t *nav, agpos_t *agp, double *xp, double *Pp)
{
    double *x,*P,*v,*H,*var,*xt,*Pt,*vr,*R,clk_est;
    int i,j,nx=agp->spp.flt.nx,nv=NFREQ*n+3,stat=0,*ind,*exc,info;
    int clk_jump_ind,reset_cnt=-1,uLSQ=1,reset_ind;
    const prcopt_t *opt=&agp->opt;
    const int max_reset_cnt=5,det_outl_cnt=1;

    trace(3,"updfilt_prc: type=%d\n",type);

    if (type==0) {
        if (!agp->spp.lsq.psta) return 0;
    }
    else if (type==1) {
        if (!agp->spp.lsq.vsta) return 0;
    }
    Pt=mat(nx,nx); P=mat(nx,nx);
    x=mat(1,nx); H=mat(nv,nx); v=mat(1,nv);
    vr=mat(nx+nv,1); R=zeros(nv,nv);
    xt=mat(nx,1); var=mat(nv,1);

    ind=imat(nv,1); exc=imat(nv,1);
    matcpy(xt,xp, 1,nx);
    matcpy(Pt,Pp,nx,nx);

sppfltrestart:
    /* reset filter status */
    memset(exc,0,sizeof(int)*NFREQ*n);
    matcpy(x,xt, 1,nx);
    matcpy(P,Pt,nx,nx);

    /* check reset filter counts */
    if (reset_cnt++>max_reset_cnt) {
        trace(2,"reset filter counts overflow: %d\n",reset_cnt);
        stat=-2;
        goto sppfltquit;
    }
    if (type) {
        /* detect rcv clock drift jump */
        if (detect_clk_drift_jump(res_func,obs,n,nav,agp,xt,Pt)) {
            goto sppfltrestart;
        }
    }
    else {
        /* detect rcv clock jump */
        if (detect_clk_jump(res_func,obs,n,nav,agp,xt,Pt)) {
            goto sppfltrestart;
        }
        /* detect ISB jump */
        if (detect_ISB_jump(res_func,obs,n,nav,agp,xt,Pt)) {
            goto sppfltrestart;
        }
    }
#if DETECT_OUTL
    /* meas residuals before exclude outliers */
    if ((nv=res_func(obs,n,nav,-1,1,uLSQ,1,x,agp,v,H,var,exc,ind))<=0) {
        trace(2,"no valid meas\n");
        stat=-1;
        goto sppfltquit;
    }
    /* detect meas outliers */
    for (i=0;i<det_outl_cnt;i++) {
        if (!detoutl_meas(v,ind,nv,exc,obs,1.5,2.5,0.001)) {
            break;
        }
    }
#endif
    /* measurement residuals */
    if ((nv=res_func(obs,n,nav,-1,1,uLSQ,1,x,agp,v,H,var,exc,ind))<=0) {
        trace(2,"no valid meas\n");
        stat=-1;
        goto sppfltquit;
    }
    /* filter model validation */
    if ((stat=val_fltmdl(type,reset_cnt,v,ind,H,var,obs,x,P,nv,nx,exc))<=0) {

        /* filter diverged and reset filter states */
        if (stat==-1) {
            if (init_flt_state(agp,obs,n,0,type?2:1,xt,Pt)) {
                goto sppfltrestart;
            }
            else {
                goto sppfltquit;
            }
        }
        /* reset filter states fail */
        else if (stat==-2) {
            trace(2,"reset filter states fail\n");
            goto sppfltquit;
        }
        /* filter model valid fail */
        else {
            trace(2,"filter model valid fail\n");
            stat=-3;
            goto sppfltquit;
        }
    }
    /* start robust ekf filter */
    for (stat=i=0;i<MAXITR;i++) {

        /* predict measurement residuals */
        if ((nv=res_func(obs,n,nav,-1,1,uLSQ,1,x,agp,v,H,var,exc,ind))<=0) {
            trace(2,"no valid meas\n");
            stat=-1;
            break;
        }
#if ADRBFILTER
        /* measurement variance */
        memset(R,0,sizeof(double)*nv*nv);
        for (j=0;j<nv;j++) R[j+j*nv]=var[j];

        /* do ekf filter */
        for (memset(R,0,sizeof(double)*nv*nv),i=0;i<nv;i++) R[i+i*nv]=var[i];
        if (adrbfilter(x,P,H,v,R,nx,nv,NULL,0)) {
            trace(2,"ekf filter error\n");
            stat=-2;
            break;
        }
#else
        /* do ekf filter */
        if (filt(x,P,NULL,NULL,H,v,var,nx,nv,1)) {
            trace(2,"ekf filter error\n");
            stat=-2;
            break;
        }
#endif
        /* post-fit residuals */
        if ((nv=res_func(obs,n,nav,-1,1,uLSQ,1,x,agp,v,H,var,exc,ind))<=0) {
            trace(2,"no valid meas\n");
            stat=-1;
            break;
        }
        /* validation of solution */
        if ((info=valsol(type,agp,v,ind,nv,var,H,xt,Pt,x,P,nx,obs,vr,&clk_jump_ind,&reset_ind))==1) {
            matcpy(xp,x, 1,nx);
            matcpy(Pp,P,nx,nx);
            stat=1;
            break;
        }
            /* solution validation fail */
        else if (info==0) {

excludeoutl:
            /* exclude outliers */
            excludeoutl(ind,vr+nx,var,nv,obs,exc);

            /* reset ekf x/P */
            matcpy(x,xt, 1,nx);
            matcpy(P,Pt,nx,nx);
            continue;
        }
        else if (info==-1) {

            /* detected clock/ISB jump */
            init_rcv_clk(clk_jump_ind,xt,Pt,nx,x,opt);
            goto sppfltrestart;
        }
        else if (info==-2) {

            /* reset ekf filter states */
            if (init_flt_state(agp,obs,n,0,reset_ind,xt,Pt)) {
                goto sppfltrestart;
            }
            else {
                goto excludeoutl;
            }
        }
    }
sppfltquit:
    free(exc); free(ind); free(var);
    free(v); free(H); free(x); free(P);
    free(Pt); free(xt); free(vr);
    return stat;
}
/* update filter using psendorange/doppler measurement-----------------------*/
static int updmeasflt(const aobsd_t *obs, int n, agpos_t *agp, const nav_t *nav, const prcopt_t *opt,
                      double *xp, double *Pp)
{
    trace(3,"updmeasflt: n=%d\n",n);
    int stat=0;

    /* update psendorange measurement */
    if (updfilt_prc(0,prng_res,obs,n,nav,agp,xp,Pp)>0) stat++;

    /* update doppler measurement */
    if (updfilt_prc(1,dopp_res,obs,n,nav,agp,xp,Pp)>0) stat++;
    return stat;
}
/* output some debug information---------------------------------------------*/
static void outputdbg(const asflt_t *flt, const aslsq_t *lsq, const agpos_t *agp)
{
#define TIME_INTTERP  1

    static FILE *fp_dbg=NULL;
    sol_t *solp=getphonesol(flt->time);
    sol_t *solr=getrtksol(flt->time);
    double dp[3],denu1[3],denu2[3],ved[3],vep[3],ver[3];
    double Q[9],pos[3];
    double denu3[3],denu4[3];
    int i;

    if (!fp_dbg) if (!(fp_dbg=fopen(".\\dbg.out","w"))) return;

    ecef2pos(lsq->rr,pos);
    xyz2enu(pos,Q);

    double gdp[3],gdv[3];
    double gep[3],gev[3];
    for (i=0;i<3;i++) gdp[i-0]=flt->x[i]-flt->xp[i];
    for (i=3;i<6;i++) gdv[i-3]=flt->x[i]-flt->xp[i];
    matmul("NN",3,1,3,1.0,Q,gdp,0.0,gep);
    matmul("NN",3,1,3,1.0,Q,gdv,0.0,gev);

    double gdp_h,gdv_h;
    gdp_h=norm(gep,2);
    gdv_h=norm(gev,2);

    if (!solp||!solr) {
        if (!solp) trace(2,"no phone solution\n");
        if (!solr) trace(2,"no agp solution\n");
        return;
    }
    double dt_r=0.0,dt_p=0.0,dt_l=0.0;

#if TIME_INTTERP
    dt_r=timediff(flt ->time,solr->time);
    dt_p=timediff(solp->time,solr->time);
    dt_l=timediff(lsq ->time,solr->time);

    if (dt_r>5.0) dt_r=0.0;
    if (dt_p>5.0) dt_p=0.0;
    if (dt_l>5.0) dt_l=0.0;
#endif

    for (i=0;i<3;i++) dp[i]=flt->x[i]-(solr->rr[i]+solr->rr[3+i]*dt_r);
    matmul("NN",3,1,3,1.0,Q,dp,0.0,denu1);

    for (i=0;i<3;i++) dp[i]=solp->rr[i]-(solr->rr[i]+solr->rr[3+i]*dt_p);
    matmul("NN",3,1,3,1.0,Q,dp,0.0,denu2);

    for (i=0;i<3;i++) dp[i]=lsq->rr[i]-(solr->rr[i]+solr->rr[3+i]*dt_l);
    matmul("NN",3,1,3,1.0,Q,dp,0.0,denu4);

    for (i=0;i<3;i++) dp[i]=flt->xp[i]-(solr->rr[i]+solr->rr[3+i]*dt_r);
    matmul("NN",3,1,3,1.0,Q,dp,0.0,denu3);

    matmul("NN",3,1,3,1.0,Q,lsq->rr+3,0.0,ved);
    matmul("NN",3,1,3,1.0,Q,solp->rr+3,0.0,vep);
    matmul("NN",3,1,3,1.0,Q,solr->rr+3,0.0,ver);

    fprintf(fp_dbg,"%10.4lf  "
                   "%8.3lf %8.3lf %8.3lf (%8.3lf %8.3lf %8.3lf, %8.3lf %8.3lf %8.3lf, %8.3lf %8.3lf, %8.3lf %8.3lf, %8.3lf %8.3lf)  %8.3lf %8.3lf %8.3lf  "
                   "%8.3lf %8.3lf %8.3lf  "
                   "%8.3lf %8.3lf %8.3lf  "
                   "%8.3lf %8.3lf %8.3lf\n",
            time2gpst(flt->time,NULL),
            denu1[0],denu1[1],denu1[2],
            denu3[0],denu3[1],denu3[2],
            denu4[0],denu4[1],denu4[2],
            lsq->sig0_p0,lsq->sig0_v0,
            lsq->sig0_p1,lsq->sig0_v1,
            gdp_h,gdv_h,
            denu2[0],denu2[1],denu2[2],
            ved[0],ved[1],ved[2],
            vep[0],vep[1],vep[2],
            ver[0],ver[1],ver[2]);
    fflush(fp_dbg);
}
/* GNSS navigation filter using android raw GNSS measurement data------------
 * args:    agpos* agp     IO  agp struct
 *          aobsd_t* obs   I   observation data
 *          int n          I   number of observation data
 *          nav_t* nav     I   navigation data
 * return : status (1: ok, 0: error)
 * --------------------------------------------------------------------------*/
extern int asppflt_ISB(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav)
{
    int nx=agp->spp.flt.nx,stat=0;
    double *xp,*Pp;

    trace(3,"asppflt_ISB: sow=%6.3lf n=%d\n",time2gpst(obs[0].time,NULL),n);

    Pp=zeros(nx,nx);
    xp=zeros(1,nx);

    /* initialize ekf estimated states */
    if (!init_flt_state(agp,obs,n,1,0,NULL,NULL)) {
        trace(2,"initialize EKF states fail\n");
        free(xp); free(Pp);
        return 0;
    }
    /* time update of ekf filter */
    if (!timeupdflt(agp,obs,n,nav,xp,Pp)) {
        trace(2,"update EKF filter fail\n");
        free(xp); free(Pp);
        return 0;
    }
    /* update measurement data */
    stat=updmeasflt(obs,n,agp,nav,&agp->opt,xp,Pp);

    /* update filter solution */
    updfltsol(obs[0].time,agp,xp,Pp,nx,stat);

#if OUT_DBG_INFO
    /* output debugs information */
    outputdbg(&agp->spp.flt,&agp->spp.lsq,agp);
#endif
    free(xp); free(Pp);
    return stat;
}
