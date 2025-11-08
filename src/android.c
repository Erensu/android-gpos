/*------------------------------------------------------------------------------
 * android.c : decode android raw gnss measurement data functions
 *
 * reference :
 *     [1] USING GNSS RAW MEASUREMENTS ON ANDROID DEVICES
 *         (https://www.gsa.europa.eu/gnss-raw-measurements-task-force.)
 *     [2] https://github.com/google/gps-measurement-tools
 *     [3] Android GNSS Raw measurement data format: see file raw/gnsslogger
 *
 * author  : sujinglan
 * version : $Revision: 1.2 $ $Date: 2008/07/14 00:05:05 $
 * history : 2022/10/09 1.0  new
 *-----------------------------------------------------------------------------*/
#include "android.h"

#ifdef ANDROID
/* constants------------------------------------------------------------------*/
#define STATE_CODE_LOCK             (1<<0)     /* tracking state: has code lock */
#define STATE_BIT_SYNC              (1<<1)     /* tracking state: has bit sync */
#define STATE_SUBFRAME_SYNC         (1<<2)     /* tracking state: has sub-frame sync. */
#define STATE_TOW_DECODED           (1<<3)     /* tracking state: has time-of-week decoded */
#define STATE_MSEC_AMBIGUOUS        (1<<4)     /* tracking state: contains millisecond ambiguity. */
#define STATE_SYMBOL_SYNC           (1<<5)     /* tracking state: has symbol sync */
#define STATE_GLO_STRING_SYNC       (1<<6)     /* tracking state: has string sync */
#define STATE_GLO_TOD_DECODED       (1<<7)     /* tracking state: has time-of-day decoded */
#define STATE_BDS_D2_BIT_SYNC       (1<<8)     /* tracking state: has D2 bit sync */
#define STATE_BDS_D2_SUBFRAME_SYNC  (1<<9)     /* tracking state: has D2 sub-frame sync */
#define STATE_GAL_E1BC_CODE_LOCK    (1<<10)    /* tracking state: has E1B/C code lock. */
#define STATE_GAL_E1C_2ND_CODE_LOCK (1<<11)    /* tracking state: has E1C secondary code lock */
#define STATE_GAL_E1B_PAGE_SYNC     (1<<12)    /* tracking state: has E1B page sync */
#define STATE_SBAS_SYNC             (1<<13)    /* tracking state: has whole second level sync */
#define STATE_TOW_KNOWN             (1<<14)    /* tracking state: has time-of-week known */
#define STATE_GLO_TOD_KNOWN         (1<<15)    /* tracking state: has time-of-day known */

#define MAXFIELD              64               /* max number of fields in a record */
#define MAXTDIFF              10.0             /* max time difference between transmited and revceive time */
#define WEEKSEC               604800           /* number of seconds in a week */
#define CONST_GPS             1                /* constellation ID: GPS */
#define CONST_BEIDOU          5                /* constellation ID: BDS */
#define CONST_GALILEO         6                /* constellation ID: GALILEO */
#define CONST_GLONASS         3                /* constellation ID: GLONASS */
#define CONST_QZSS            4                /* constellation ID: QZSS */
#define CONST_SBAS            2                /* constellation ID: SBAS */
#define BDST_TO_GPST          14.0             /* leap seconds difference between BDST and GPST */
#define GLOT_TO_UTC           10800.0          /* time difference between GLOT and UTC in seconds */
#define CHK_GPS_WEEK          0                /* check gps week */
#define END_FLAG_STR          "EPOCHEND"     /* epoch end flag for gnss observation */

/* get gps frequency index----------------------------------------------------*/
static int getfrq_gps(const double fhz, double *freq)
{
    static const double gps_frq[3]={FREQ1,FREQ2,FREQ5};
    int i;
    for (i=0;i<3;i++) {
        if (fabs(gps_frq[i]-fhz)<1E4) {
            *freq=gps_frq[i];
            return i;
        }
    }
    return 0;
}
/* get qzs frequency index----------------------------------------------------*/
static int getfrq_qzs(const double fhz, double *freq)
{
    static const double qzs_frq[5]={FREQ1,FREQ2,FREQ5,FREQ6};
    int i;
    for (i=0;i<4;i++) {
        if (fabs(qzs_frq[i]-fhz)<1E4) {
            *freq=qzs_frq[i];
            return i;
        }
    }
    return 0;
}
/* get gal frequency index----------------------------------------------------*/
static int getfrq_gal(const double fhz, double *freq)
{
    int i;
    static const double gal_frq[3]={FREQ1,FREQ7,FREQ5};

    for (i=0;i<3;i++) {
        if (fabs(gal_frq[i]-fhz)<1E4) {
            *freq=gal_frq[i];
            return i;
        }
    }
    return 0;
}
/* get bds frequency index----------------------------------------------------*/
static int getfrq_bds(const double fhz, double *freq)
{
    int i;
    static const double bds_frq[6]={FREQ1,FREQ1_CMP,FREQ2_CMP,FREQ5,FREQ3_CMP,FREQ8};

    for (i=0;i<6;i++) {
        if (fabs(bds_frq[i]-fhz)<1E4) {
            *freq=bds_frq[i];
            return i;
        }
    }
    return 0;
}
/* get glo frequency index----------------------------------------------------*/
static int getfrq_glo(const double fhz, int svid, double *freq)
{
    static const double L1[]={1.59290E9,1.61000E9};
    static const double L2[]={1.23780E9,1.25680E9};
    static const double L3[]={1.19035E9,1.21223E9};

    *freq=fhz;
    if      (fhz>=L1[0]&&fhz<=L1[1]) return 0; /* L1 */
    else if (fhz>=L2[0]&&fhz<=L2[1]) return 1; /* L2 */
    else if (fhz>=L3[0]&&fhz<=L3[1]) return 2; /* L3 */
    return 0;
}
/* gnss signal index----------------------------------------------------------*/
static int sigindex(double fhz,int ctype, int *find, int stat, int svid, double *freq)
{
    int frq=0,code=CODE_L1C;

    switch (ctype) {
        case CONST_GPS: { /* GPS */
            frq=getfrq_gps(fhz,freq);
            if (frq==0) {code=obs2code("1C",NULL); break;}
            if (frq==1) {code=obs2code("2W",NULL); break;}
            if (frq==2) {code=obs2code("5Q",NULL); break;}
        }
        case CONST_GALILEO: { /* GALILEO */
            frq=getfrq_gal(fhz,freq);
            if (frq==0) {code=obs2code("1X",NULL); break;}
            if (frq==1) {code=obs2code("7X",NULL); break;}
            if (frq==2) {code=obs2code("5X",NULL); break;}
            if (frq==3) {code=obs2code("6X",NULL); break;}
            if (frq==4) {code=obs2code("8X",NULL); break;}
        }
        case CONST_GLONASS: { /* GLONASS */
            frq=getfrq_glo(fhz,svid,freq);
            if (frq==0) {code=obs2code("1C",NULL); break;}
            if (frq==1) {code=obs2code("2C",NULL); break;}
            if (frq==2) {code=obs2code("3X",NULL); break;}
        }
        case CONST_BEIDOU: { /* BDS */
            frq=getfrq_bds(fhz,freq);
            if (frq==0) {code=obs2code("1C",NULL); break;} /* B1C */
            if (frq==1) {code=obs2code("2I",NULL); break;} /* B1I */
            if (frq==2) {code=obs2code("7I",NULL); break;} /* B2I/B2b */
            if (frq==3) {code=obs2code("5P",NULL); break;} /* B2a */
            if (frq==4) {code=obs2code("6I",NULL); break;} /* B3 */
            if (frq==5) {code=obs2code("8P",NULL); break;} /* B2ab */
        }
        case CONST_QZSS: { /* QZSS */
            frq=getfrq_qzs(fhz,freq);
            if (frq==0) {code=obs2code("1C",NULL); break;}
            if (frq==1) {code=obs2code("2W",NULL); break;}
            if (frq==2) {code=obs2code("5Q",NULL); break;}
            if (frq==3) {code=obs2code("6Q",NULL); break;}
        }
        case CONST_SBAS: { /* SBAS */
            frq=getfrq_gps(fhz,freq);
            if (frq==0) {code=obs2code("1C",NULL); break;}
            if (frq==2) {code=obs2code("5I",NULL); break;}
        }
    }
    if (find) *find=frq;
    return code;
}
/* get satellite id-----------------------------------------------------------*/
static int satid(int ctype, int svid, int *sys)
{
    int prn=0;
    *sys=SYS_NONE;

    if (ctype==CONST_GLONASS&&(svid<=106&&svid>=93)) {
        svid=MAXPRNGLO1+svid-93;
    }
    switch (ctype) {
        case CONST_GPS    : *sys=SYS_GPS; prn=svid; break;
        case CONST_GLONASS: *sys=SYS_GLO; prn=svid; break;
        case CONST_BEIDOU : *sys=SYS_CMP; prn=svid; break;
        case CONST_GALILEO: *sys=SYS_GAL; prn=svid; break;
        case CONST_QZSS   : *sys=SYS_QZS; prn=svid; break;
        case CONST_SBAS   : *sys=SYS_SBS; prn=svid; break;
        default: {
            trace(2,"no supported constellation type\n");
        }
    }
    return satno(*sys,prn);
}
/* update observation data buffer---------------------------------------------*/
static void updobsbuf(araw_t *raw)
{
    static aobsd_t obs0={0};
    int i,j;

    trace(3,"update observation data buffer: n=%d\n",raw->obuf.n);

    raw->obs.n=0;

    /* update observation data */
    for (i=0;i<raw->obuf.n&&i<raw->obs.nmax;i++) {
        raw->obs.data[i]=raw->obuf.data[i];
        raw->obs.n++;

        for (j=0;j<NFREQ;j++) {
            if (raw->obuf.data[i].code[j]==CODE_NONE) continue;
            if (raw->obuf.data[i].stat[j]&OBSQ_COARSE_TIME) {
                raw->obs.data[i].time=raw->obs.data[0].time;
            }
        }
    }
    for (i=0;i<MAXOBS;i++) {
        raw->obuf.data[i]=obs0;
    }
    /* reset counts of obs buffer */
    raw->obuf.n=0;
}
/* determine whether it is the same epoch ------------------------------------*/
static int issamepoch(araw_t *raw, double *tp, const double tc)
{
    if (fabs(*tp-tc)>10) {

        /* new epoch and update observation data buffer */
        updobsbuf(raw);
        *tp=tc;
        raw->obsflag=1;
        return 0;
    }
    raw->obsflag=0;
    return 1;
}
static int chksameepoch(araw_t *raw, const double tc)
{
    /* determine whether it is the same epoch */
    if (raw->ptime<=0) {raw->ptime=tc; return 1;}
    else {
        return issamepoch(raw,&raw->ptime,tc);
    }
}
/* check gnss measurement's tracking state------------------------------------*/
static int chktrkstat(int ctype, int stat)
{
    switch (ctype) {
        case CONST_GPS    : {
            if (!(stat&STATE_CODE_LOCK            )) return 0;
            if (!(stat&STATE_BIT_SYNC             )) return 0;
            if (!(stat&STATE_SUBFRAME_SYNC        )) return 0;
            if (!(stat&STATE_TOW_DECODED          )) return 0;
            return 1;
        }
        case CONST_GALILEO: {
            if (!(stat&STATE_CODE_LOCK            )) return 0;
            if (!(stat&STATE_GAL_E1C_2ND_CODE_LOCK)) return 0;
            if (!(stat&STATE_TOW_DECODED          )) return 0;
            return 1;
        }
        case CONST_BEIDOU:  {
            if (!(stat&STATE_CODE_LOCK            )) return 0;
            if (!(stat&STATE_BIT_SYNC             )) return 0;
            if (!(stat&STATE_SUBFRAME_SYNC        )) return 0;
            if (!(stat&STATE_TOW_DECODED          )) return 0;
            return 1;
        }
        case CONST_GLONASS: {
            if (!(stat&STATE_CODE_LOCK            )) return 0;
            if (!(stat&STATE_GLO_TOD_DECODED      )) return 0;
            return 1;
        }
        case CONST_QZSS:    {
            if (!(stat&STATE_CODE_LOCK            )) return 0;
            if (!(stat&STATE_BIT_SYNC             )) return 0;
            if (!(stat&STATE_SUBFRAME_SYNC        )) return 0;
            if (!(stat&STATE_TOW_DECODED          )) return 0;
            return 1;
        }
        case CONST_SBAS:    {
            if (!(stat&STATE_CODE_LOCK            )) return 0;
            if (!(stat&STATE_BIT_SYNC             )) return 0;
            if (!(stat&STATE_SUBFRAME_SYNC        )) return 0;
            if (!(stat&STATE_TOW_DECODED          )) return 0;
            return 1;
        }
    }
    return 0;
}
/* compute glonass satellite transmitted time---------------------------------*/
static gtime_t cmpglott(const double tt, const double tr, const int week)
{
    double tod=fmod(tt,86400.0),dow=(int)floor(tr/86400.0);

    gtime_t gt=gpst2time(week,tr);
    gt=timeadd(gt,-GLOT_TO_UTC+dow*86400.0-tr);
    gt=timeadd(gt,tod);
    return utc2gpst(gt);
}
/* adjust received sv time to gpst--------------------------------------------*/
static gtime_t adjttr(const int ctype, const double tt, const double tr, const int week)
{
    double tt_;
    if      (ctype==CONST_BEIDOU ) {tt_=tt+BDST_TO_GPST; return gpst2time(week,tt_);}
    else if (ctype==CONST_GLONASS) {
        return cmpglott(tt,tr,week);
    }
    return gpst2time(week,tt);
}
/* check week crossover-------------------------------------------------------*/
static void chkweek(gtime_t *tt, const gtime_t tr)
{
    if (fabs(timediff(tr,*tt))>WEEKSEC/2.0) *tt=timeadd(*tt,WEEKSEC);
}
/* adjust observation index --------------------------------------------------*/
static int adjobsindex(int sat, const aobs_t *obs, int *iobs)
{
    int i;

    if (sat==0) return -1;
    for (i=0;i<obs->n;i++) {
        if (obs->data[i].sat==sat) {*iobs=i; return 1;}
    }
    return -1;
}
/* update observation data uncertainty----------------------------------------*/
static void updobsunc(aobsd_t *obs, const char **val, int nval, int fi, const double fhz)
{
    obs->unc.time [fi]=atof(val[ 4])*CLIGHT*1E-9;
    obs->unc.bias [fi]=atof(val[ 7])*CLIGHT*1E-9;
    obs->unc.rt   [fi]=atof(val[15])*CLIGHT*1E-9;
    obs->unc.drift[fi]=atof(val[ 9])*CLIGHT*1E-9;
    obs->unc.pr   [fi]=atof(val[18]);
    obs->unc.cp   [fi]=atof(val[25]);
    obs->unc.adr  [fi]=atof(val[21]);
    obs->unc.tt   [fi]=atoll(val[14]);
    obs->unc.stat [fi]=atoi(val[13]);
    obs->unc.adr_stat[fi]=atoi(val[19]);
    obs->unc.snr[fi]=atof(val[27]);
    obs->unc.frq[fi]=fhz;
    obs->unc.svid=atoi(val[11]);
}
/* update status of multipath-------------------------------------------------*/
static void updmultipath(const char **val, aobsd_t *obs, int f)
{
    switch (atoi(val[26])) {
        case 0: obs->stat[f]|=OBSQ_MULTP_UNKNOWN  ; break;
        case 1: obs->stat[f]|=OBSQ_MULTP_DETECT   ; break;
        case 2: obs->stat[f]|=OBSQ_MULTP_NOTDETECT; break;
        default:
            obs->stat[f]|=OBSQ_MULTP_NOTDETECT; break;
    }
}
/* decode gps raw measurement data--------------------------------------------*/
static int decode_meas_raw(araw_t *raw, const char **val, int nval, int ctype)
{
    long long tr_ns,week_ns,fullbias,bias;
    double tr,tt,dt,fhz,freq=0.0;
    int fi=0,code,week,stat=atoi(val[13]),svid=atoi(val[11]),sys,sat,iobs=0,amb_flag=1,week_cur;
    int corse_obst=0;
    static int upd_flag=1;
    char prn[8];
    gtime_t gtr={0},gtt={0};

    trace(3,"decode_gps_raw: ctype=%2d svid=%3d stat=%4d\n",ctype,svid,stat);

    /* GPS week number */
    week=(int)floor(-atoll(val[5])*1E-9/WEEKSEC);
    if (week<=0) {
        week=raw->week;
#if CHK_GPS_WEEK
        if (week<=0) time2gpst(timeget(),&week);
#endif
        corse_obst=1;
    }
#if CHK_GPS_WEEK
    time2gpst(timeget(),&week_cur);
    if (week_cur!=week) {
        trace(2,"gps week is unreliable: week=%d week_cir=%d\n",week,week_cur);
        week=week_cur;
    }
#endif
    fhz=atof(val[22]);
    sat=satid(ctype,svid,&sys);

    if (atoll(val[2])<=0) {
        trace(2,"value of timenanos should be >0\n");
        return 0;
    }
    /* build android version is O */
    if (nval>30) fhz=atof(val[30]);

    /* check hardware clock discontinuities */
    if (raw->hwclk0) {
        if (atoll(val[10])!=raw->hwclk0) {
            raw->fullbias0=raw->bias0=0;
            raw->hwclk0=0;
        }
    }
    /* check hardware clock jump or overflow */
    if (raw->timenanos) {
        if (atoll(val[2])-raw->timenanos<0.0||fabs(atoll(val[2])-raw->timenanos)>1.5E9) {
            raw->fullbias0=raw->bias0=0;
            raw->hwclk0=0;
        }
    }
    /* check biasnanos value */
    if (fabs(atoll(val[7])*CLIGHT*1E-9)>1E6) {
        raw->fullbias0=raw->bias0=0;
        raw->hwclk0=0;
    }
#if 0
    if (!(fullbias=raw->fullbias0)) fullbias=atoll(val[5]);
    if (!(bias=raw->bias0)) bias=atoll(val[6]);
#else
    fullbias=atoll(val[5]);
    bias=atoll(val[6]);
#endif
    /* phone receive time in gpst */
    week_ns=week*WEEKSEC*1E9;
    tr_ns=atoll(val[2])-fullbias-week_ns; /* ns */
    tr=(tr_ns-atoll(val[12])-bias)*1E-9;
    gtr=gpst2time(week,tr);

    /* determine whether it is the same epoch */
    if (upd_flag) {
        chksameepoch(raw,atof(val[2]));
    }
    if (corse_obst&&nval>31) {
        if (strstr(val[nval-1],END_FLAG_STR)) {
            gtr=gpst2time(week,atof(val[nval-2]));
        }
        else {
            gtr=gpst2time(week,atof(val[nval-1]));
        }
    }
    if (!chktrkstat(ctype,stat)) {
        trace(2,"current state of satellite is in wrong track status. ctype=%2d "
                "svid=%3d stat=%4d time=%lld\n",
                ctype,svid,stat,atoll(val[2]));

        switch (ctype) {
            case CONST_GALILEO:
                if (stat&STATE_GAL_E1BC_CODE_LOCK) amb_flag=0;
                break;
            case CONST_GLONASS:
            case CONST_BEIDOU:
            case CONST_GPS:
            case CONST_QZSS:
                if (stat&STATE_CODE_LOCK) amb_flag=0;
                break;
            default:
                amb_flag=1; break;
        }
    }
    /* satellite time at the measurement time */
    tt=atoll(val[14])*1E-9;

    /* check week crossover */
    chkweek(&gtt,gtr);

    /* adjust transmitted sv time to gpst(s) */
    gtt=adjttr(ctype,tt,tr,week);

    /* transmission time interval */
    if (fabs(dt=fmod(timediff(gtr,gtt),WEEKSEC))>MAXTDIFF||dt<=50.0*1E-3) {
        trace(2,"too large transmission time interval: dt(ms)=%12.8lf\n",dt*1E3);

        /* fix 1ms ambiguity */
        dt=0.0;
        amb_flag=0;
    }
    /* update receiver observation */
    code=sigindex(fhz,ctype,&fi,stat,svid,&freq);

    /* observation index adjust for dual frequency */
    iobs=adjobsindex(sat,&raw->obuf,&iobs)<0?raw->obuf.n:iobs;

    /* update observation uncertainty */
    updobsunc(&raw->obuf.data[iobs],val,nval,fi,fhz);

    /* too larger uncertainty */
    if (raw->obuf.data[iobs].unc.rt[fi]>1E3) {
        amb_flag=0;
    }
    /* observation data pseudorange rate (m/s) */
    if (fhz>0.0) {
        raw->obuf.data[iobs].D[fi]=(float)(-atof(val[17])/CLIGHT*fhz);
    }
    /* build android version is O */
    if (nval>=30) {

        /* gain amplifier adjusting the power of the incoming signal */
        raw->obuf.data[iobs].AGC[fi]=atof(val[29]);
    }
    /* hardware clock discontinuity count */
    raw->obuf.data[iobs].hwclkdc[fi]=atoi(val[10]);

    /* receiver clock drift per second (m/s) */
    raw->obuf.data[iobs].CD[fi]=atof(val[8])*1E-9*CLIGHT;

    /* pseudorange (m) */
    raw->obuf.data[iobs].P   [fi]=CLIGHT*dt;
    raw->obuf.data[iobs].L   [fi]=atof(val[20])*fhz/CLIGHT;
    raw->obuf.data[iobs].SNR [fi]=(unsigned char)(atof(val[16])*4.0+0.5);
    raw->obuf.data[iobs].Pr  [fi]=atof(val[17]);
    raw->obuf.data[iobs].ADR [fi]=atof(val[20]);
    raw->obuf.data[iobs].code[fi]=code;
    raw->obuf.data[iobs].time=gtr;
    raw->obuf.data[iobs].sat =sat;
    raw->obuf.data[iobs].rcv =1;
    raw->obuf.data[iobs].stat[fi]=OBSQ_RAW;
    raw->obuf.data[iobs].type=OBST_CODE;

    raw->obuf.data[iobs].TimeNanos=atoll(val[2]);
    raw->obuf.data[iobs].FullBiasNanos=atoll(val[5]);

    if (!amb_flag) {
        raw->obuf.data[iobs].stat[fi]|=OBSQ_PR_FIXAMB;
        trace(2,"fix observation ms ambiguity: sat=%d\n",sat);
    }
    if (ctype==CONST_GLONASS&&(svid<=106&&svid>=93)) {
        raw->obuf.data[iobs].stat[fi]|=OBSQ_GLO_FCN;
        trace(2,"fix glonass fcn ambiguity: svid=%d\n",svid);
    }
    if (corse_obst) {
        raw->obuf.data[iobs].stat[fi]|=OBSQ_COARSE_TIME;
        trace(2,"use observation coarse time: sat=%d\n",sat);
    }
    /* multipath status */
    updmultipath(val,&raw->obuf.data[iobs],fi);
    raw->obuf.n+=(iobs==raw->obuf.n);

    /* fullbias at first epoch */
    if (!raw->fullbias0) {
        if (ctype==CONST_GPS&&chktrkstat(ctype,stat)&&amb_flag) {
            trace(3,"fullbias0=%lld bias0=%lld hwclk0=%d\n",atoll(val[5]),atoll(val[6]),atoi(val[10]));

            raw->fullbias0=atoll(val[5]);
            raw->bias0=atoll(val[6]);
            raw->hwclk0=atoi(val[10]);
        }
    }
    raw->timenanos=atoll(val[2]);
    satno2id(sat,prn);

    trace(3,"track satellite(%8.3lf): %4s code=%2d f=%d timenanos=%lld fullbias=%lld tt=%10.6f ambflag=%d\n",
          time2gpst(gtr,NULL),
          prn,code,fi,atoll(val[2]),atoll(val[5]),tt,amb_flag);

    if (strstr(val[nval-1],END_FLAG_STR)) {
        updobsbuf(raw);
        upd_flag=0;
        raw->obsflag=1;
    }
    else if (!upd_flag) {
        raw->obsflag=0;
    }
    return raw->obsflag;
}
/* decode gnss raw measurement data-------------------------------------------*/
static int decode_graw(araw_t *raw,char *pbuf)
{
    char *val[MAXFIELD],*p,*q;
    int n=0,ctype;

    trace(3,"decode_graw:\n");

    /* parse fields */
    for (p=pbuf;*p&&n<MAXFIELD;p=q+1) {
        if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {val[n++]=p; *q='\0';}
        else break;
    }
    ctype=atoi(val[28]);
    switch (ctype) {
        case CONST_GPS    : return decode_meas_raw(raw,val,n,ctype);
        case CONST_GALILEO: return decode_meas_raw(raw,val,n,ctype);
        case CONST_GLONASS: return decode_meas_raw(raw,val,n,ctype);
        case CONST_BEIDOU : return decode_meas_raw(raw,val,n,ctype);
        case CONST_QZSS   : return decode_meas_raw(raw,val,n,ctype);
        case CONST_SBAS   : return decode_meas_raw(raw,val,n,ctype);
        default:
            return 0;
    }
}
/* decode gnss navigation message data----------------------------------------*/
static int decode_gnav(araw_t *raw, const char *pbuf)
{
    trace(3,"decode_gnav:\n");
    return 0;
}
/* decode gnss solution provided by android message data----------------------*/
static int decode_gsol(araw_t *raw, char *pbuf)
{
    char *val[MAXFIELD],*p,*q;
    double pos[3],tutc;
    int n=0;

    trace(3,"decode_gsol:\n");

    /* parse fields */
    for (p=(char*)pbuf;*p&&n<MAXFIELD;p=q+1) {
        if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {val[n++]=p; *q='\0';}
        else break;
    }
    pos[0]=atof(val[2])*D2R; /* latitude (rad) */
    pos[1]=atof(val[3])*D2R; /* longitude (rad) */
    pos[2]=atof(val[4]);     /* height (m) */
    pos2ecef(pos,raw->sol.rr);

    /* covariance */
    raw->sol.qr[0]=raw->sol.qr[1]=raw->sol.qr[2]=atof(val[6])/3.0;

    /* solution time */
    tutc=atof(val[n-1])*1E-3;
    raw->sol.time.time=(int)tutc;
    raw->sol.time.sec =tutc-(int)tutc;
    raw->sol.time=utc2gpst(raw->sol.time);
    raw->sol.stat=SOLQ_CHIP;
    return 4;
}
/* decode NLP solution-------------------------------------------------------*/
static int decode_nlp(araw_t *raw, char *pbuf)
{
    char *val[MAXFIELD],*p,*q;
    double pos[3],tutc;
    int n=0;

    trace(3,"decode_nlp:\n");

    /* parse fields */
    for (p=(char*)pbuf;*p&&n<MAXFIELD;p=q+1) {
        if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {val[n++]=p; *q='\0';}
        else break;
    }
    pos[0]=atof(val[2])*D2R; /* latitude (rad) */
    pos[1]=atof(val[3])*D2R; /* longitude (rad) */
    pos[2]=atof(val[4]);     /* height (m) */
    pos2ecef(pos,raw->sol.rr);

    /* covariance */
    raw->sol.qr[0]=raw->sol.qr[1]=raw->sol.qr[2]=atof(val[6])/3.0;

    /* solution time */
    tutc=atof(val[n-1])*1E-3;
    raw->sol.time.time=(int)tutc;
    raw->sol.time.sec =tutc-(int)tutc;
    raw->sol.time=utc2gpst(raw->sol.time);
    raw->sol.stat=SOLQ_NLP;
    return 4;
}
/* sync a newline------------------------------------------------------------*/
static int syncnewline(const unsigned char *buff, int nb)
{
    if (buff[nb-1]=='\n'||(buff[nb-2]=='\r'&&buff[nb-1]=='\n')) return 1;
    return 0;
}
/* clear message buffer------------------------------------------------------*/
static void clearbuff(araw_t *raw)
{
    memset(raw->buff,'\0',sizeof(raw->buff));
    raw->nbyte=0;
}
/* input android raw message from stream --------------------------------------
 * fetch next android raw data and input a mesasge from stream
 * args   : raw_t *raw         IO  receiver raw data control struct
 *          unsigned char data I   stream data (1 byte)
 * ---------------------------------------------------------------------------*/
extern int input_android(araw_t *raw, unsigned char data)
{
    char *pbuf=NULL;
    int ret=0;

    raw->buff[raw->nbyte++]=data;
    if (raw->nbyte>=MAXRAWLEN) { /* buffer overflow and reset decoder */
        clearbuff(raw);
        return 0;
    }
    /* sync an new line if detected */
    if (!syncnewline(raw->buff,raw->nbyte)) {
        return 0;
    }
    /* decode raw data */
    if      ((pbuf=strstr((char*)raw->buff,  "#"))) ret=0;
    else if ((pbuf=strstr((char*)raw->buff,"Raw"))) ret=decode_graw(raw,pbuf);
    else if ((pbuf=strstr((char*)raw->buff,"Nav"))) ret=decode_gnav(raw,pbuf);
    else if ((pbuf=strstr((char*)raw->buff,"Fix"))) ret=decode_gsol(raw,pbuf);
    else if ((pbuf=strstr((char*)raw->buff,"NLP"))) ret=decode_nlp (raw,pbuf);
    else ret=0;

    clearbuff(raw);
    return ret;
}
/* input android raw message from file ----------------------------------------
 * fetch next android raw data and input a message from file
 * args   : raw_t  *raw   IO     receiver raw data control struct
 *          FILE   *fp    I      file pointer
 * return : status(-2: end of file, -1...9: same as above)
 *-----------------------------------------------------------------------------*/
extern int input_androidf(araw_t  *raw, FILE *fp)
{
    int i,data,ret;

    for (i=0;i<4096;i++) {
        if ((data=fgetc(fp))==EOF) return -2;
        if ((ret=input_android(raw,(unsigned char)data))) return ret;
    }
    /* return at every 4k bytes */
    return 0;
}
#else
/* input Android raw message from stream --------------------------------------
 * fetch next android raw data and input a mesasge from stream
 * args   : raw_t *raw         IO  receiver raw data control struct
 *          unsigned char data I   stream data (1 byte)
 * ---------------------------------------------------------------------------*/
extern int input_android(araw_t *raw, unsigned char data)
{
    return -2;
}
/* input Android raw message from file ----------------------------------------
 * fetch next android raw data and input a message from file
 * args   : raw_t  *raw   IO     receiver raw data control struct
 *          FILE   *fp    I      file pointer
 * return : status(-2: end of file, -1...9: same as above)
 *-----------------------------------------------------------------------------*/
extern int input_androidf(araw_t  *raw, FILE *fp)
{
    return -2;
}
#endif

/* Android observation data convert to observation data------------------------*/
extern void aobs2obs(const aobsd_t *aobs, obsd_t *obs)
{
    int i;

    memset(obs,0,sizeof(*obs));

    for (i=0;i<NFREQ+NEXOBS;i++) {
        obs->code[i]=aobs->code[i];
        obs->SNR[i]=aobs->SNR[i];
        obs->LLI[i]=aobs->LLI[i];
    }
    obs->time=aobs->time;
    obs->sat=aobs->sat;
    obs->rcv=aobs->rcv;
    matcpy(obs->P,aobs->P,1,NFREQ+NEXOBS);
    matcpy(obs->L,aobs->L,1,NFREQ+NEXOBS);
    matcpy(obs->D,aobs->D,1,NFREQ+NEXOBS);
}
/* observation data convert to Android observation data------------------------*/
extern void obs2aobs(const obsd_t *obs, aobsd_t *aobs)
{
    int i;

    memset(aobs,0,sizeof(*aobs));

    for (i=0;i<NFREQ+NEXOBS;i++) {
        aobs->code[i]=obs->code[i];
        aobs->SNR[i]=obs->SNR[i];
        aobs->LLI[i]=obs->LLI[i];
    }
    aobs->time=obs->time;
    aobs->sat=obs->sat;
    aobs->rcv=obs->rcv;
    matcpy(aobs->P,obs->P,1,NFREQ+NEXOBS);
    matcpy(aobs->L,obs->L,1,NFREQ+NEXOBS);
    matcpy(aobs->D,obs->D,1,NFREQ+NEXOBS);
}
/* initial Android raw data struct----------------------------------------------*/
extern void initaraw(araw_t *raw)
{
    aobsd_t data0={{0}};
    eph_t  eph0 ={0,-1,-1};
    geph_t geph0={0,-1};
    int i;

    memset(raw,0,sizeof(*raw));
    raw->obs.data =calloc(MAXOBS,sizeof(aobsd_t));
    raw->obuf.data=calloc(MAXOBS,sizeof(aobsd_t));
    raw->obs.n =0; raw->obs.nmax =MAXOBS;
    raw->obuf.n=0; raw->obuf.nmax=MAXOBS;

    raw->nav.eph =NULL;
    raw->nav.geph=NULL;

    if (!(raw->nav.eph =(eph_t *)malloc(sizeof(eph_t )*MAXSAT))||
        !(raw->nav.geph=(geph_t*)malloc(sizeof(geph_t)*NSATGLO))) {
        freearaw(raw);
        return;
    }
    raw->obs.n =0; raw->obs.nmax =MAXOBS;
    raw->obuf.n=0; raw->obuf.nmax=MAXOBS;
    raw->nav.n =MAXSAT;
    raw->nav.ng=NSATGLO;
    for (i=0;i<MAXOBS ;i++) raw->obs.data [i]=data0;
    for (i=0;i<MAXOBS ;i++) raw->obuf.data[i]=data0;
    for (i=0;i<MAXSAT ;i++) raw->nav.eph  [i]=eph0;
    for (i=0;i<NSATGLO;i++) raw->nav.geph [i]=geph0;

    raw->fullbias0=0;
    raw->bias0=0;
}
/* free Android raw data struct-------------------------------------------------*/
extern void freearaw(araw_t *raw)
{
    free(raw->obs.data ); raw->obs.data =NULL; raw->obs.n =0;
    free(raw->obuf.data); raw->obuf.data=NULL; raw->obuf.n=0;
    free(raw->nav.eph  ); raw->nav.eph  =NULL; raw->nav.n =0;
    free(raw->nav.geph ); raw->nav.geph =NULL; raw->nav.ng=0;
}