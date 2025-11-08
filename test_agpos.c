#include <stdlib.h>
#include "src/agpos.h"

#define SQR(x)      ((x)*(x))

#define READ_NLP_SOLS   1
#define READ_PHONE_SOLS 1
#define READ_RKT_SOLS   1
#define TIME_INTTERP    1
#define LEVER_ARM       1

#define FMT_POS     0x01
#define FMT_NMA     0x02
#define FMT_SOL     0x03

static solbuf_t solbuf_phone={0};
static solbuf_t solbuf_nlp={0};
static solbuf_t solbuf_rtk={0};

extern sol_t* getphonesol(gtime_t time);
extern sol_t* getnlpsol(gtime_t time);
extern sol_t* getrtksol(gtime_t time);
extern int sort_solbuf(solbuf_t *solbuf);

static sol_t* get_sol(gtime_t time, const solbuf_t *solbuf)
{
    sol_t *sol1=NULL,*sol2=NULL;
    int i,j;

    for (i=0;i<solbuf->n;i++) {
        if (fabs(timediff(time,solbuf->data[i].time))<0.1) {
            if (i-1>=0) {
                sol1=&solbuf->data[i  ];
                sol2=&solbuf->data[i-1];
                if (fabs(timediff(sol1->time,sol2->time))>5.0) {
                    return sol1;
                }
                for (j=0;j<3;j++) {
                    sol1->rr[j+3]=(sol2->rr[j]-sol1->rr[j])/timediff(sol2->time,sol1->time);
                }
            }
            return sol1;
        }
    }
    return NULL;
}
static void outputkml(const solbuf_t *phone, const solbuf_t *demo, const solbuf_t *rtk)
{
    const char *demokml=".\\demo.kml";
    const char *phonekml=".\\phone.kml";
    const char *rtkkml=".\\rtk.kml";

    sol2kml(demokml,demo,1,3,0,0);
    sol2kml(phonekml,phone,1,4,0,0);
    sol2kml(rtkkml,rtk,1,5,0,0);
}

static void lever_arm(const double venu[3], const double *pos, double dxyz[3],
                      double arm_horizental)
{
    double arm_vertical=0.0;
    double speed=sqrt(venu[0]*venu[0]+venu[1]*venu[1]);
    double denu[3]={0};

    if (speed<2.0) return;
    denu[0]=arm_horizental*venu[0]/speed;
    denu[1]=arm_horizental*venu[1]/speed;
    denu[2]=arm_vertical;
    enu2ecef(pos,denu,dxyz);
}
static int cmpres(const void *p1, const void *p2)
{
    double *q1=(double*)p1,*q2=(double*)p2;
    return ((*q1-*q2)<0)?-1:1;
}
static int cmpprc(const solbuf_t *phone, const solbuf_t *demo, const solbuf_t *rtk)
{
    double *dv_chip_rtk,*dh_chip_rtk;
    double *dv_demo_rtk,*dh_demo_rtk;
    double pos[3],dr[3],e1[3],e2[3];
    gtime_t time;
    int i,j,n,nmax;

    if (phone->n<=0||demo->n<=0||rtk->n<=0) return 0;
    FILE *fp=fopen(".\\stat.out","w");

    nmax=MAX(demo->n,phone->n);
    dv_chip_rtk=mat(1,nmax);
    dh_chip_rtk=mat(1,nmax);
    dv_demo_rtk=mat(1,nmax);
    dh_demo_rtk=mat(1,nmax);

    /* start compare three solutions */
    for (n=i=0;i<demo->n;i++) {
        time=demo->data[i].time;

        sol_t *sol_p=getphonesol(time);
        sol_t *sol_r=getrtksol(time);
        sol_t *sol_d=&demo->data[i];

        if (!sol_r||!sol_p) continue;
        ecef2pos(sol_r->rr,pos);

        double lever_arm_dxyz[3]={0.0};
        double dt_r=0.0,dt_p=0.0;

#if TIME_INTTERP
        dt_r=timediff(time,sol_r->time);
        dt_p=timediff(sol_p->time,sol_r->time);
#endif

#if LEVER_ARM
        double arm_horizental=2.0,venu[3]={0.0};
        ecef2enu(pos,&sol_r->rr[3],venu);
        lever_arm(venu,pos,lever_arm_dxyz,arm_horizental);
#endif
        if (dt_r>5.0) dt_r=0.0;
        if (dt_p>5.0) dt_p=0.0;

        for (j=0;j<3;j++) {
            dr[j]=sol_d->rr[j]-(sol_r->rr[j]+sol_r->rr[3+j]*dt_r+lever_arm_dxyz[j]);
        }
        ecef2enu(pos,dr,e1);
        for (j=0;j<3;j++) {
            dr[j]=sol_p->rr[j]-(sol_r->rr[j]+sol_r->rr[3+j]*dt_p+lever_arm_dxyz[j]);
        }
        ecef2enu(pos,dr,e2);

        dh_demo_rtk[n]=sqrt(SQR(e1[0])+SQR(e1[1]));
        dv_demo_rtk[n]=fabs(e1[2]);

        dh_chip_rtk[n]=sqrt(SQR(e2[0])+SQR(e2[1]));
        dv_chip_rtk[n]=fabs(e2[2]);
        n++;
    }
    if (n<=0) {
        perror("solution comparison fail\n");
        free(dv_chip_rtk);
        free(dh_chip_rtk);
        free(dv_demo_rtk);
        free(dh_demo_rtk);
        fclose(fp); return 0;
    }
    /* output compare solution results */
    qsort(dv_demo_rtk,n,sizeof(double),cmpres);
    qsort(dh_demo_rtk,n,sizeof(double),cmpres);
    qsort(dv_chip_rtk,n,sizeof(double),cmpres);
    qsort(dh_chip_rtk,n,sizeof(double),cmpres);

    fprintf(fp,"numbers of epochs/s: %d\n",n);
    fprintf(fp,"horizontal: Demo-RTK    Phone-RTK\n");
    fprintf(fp,"  68%%/m   : %6.3lf      %6.3lf\n",dh_demo_rtk[(int)(n*0.68)],dh_chip_rtk[(int)(n*0.68)]);
    fprintf(fp,"  95%%/m   : %6.3lf      %6.3lf\n",dh_demo_rtk[(int)(n*0.95)],dh_chip_rtk[(int)(n*0.95)]);
    fprintf(fp,"\n");
    fprintf(fp,"vertical  : Demo-RTK    Phone-RTK\n");
    fprintf(fp,"  68%%/m   : %6.3lf      %6.3lf\n",dv_demo_rtk[(int)(n*0.68)],dv_chip_rtk[(int)(n*0.68)]);
    fprintf(fp,"  95%%/m   : %6.3lf      %6.3lf\n",dv_demo_rtk[(int)(n*0.95)],dv_chip_rtk[(int)(n*0.95)]);
    fflush(fp);

    perror("solution comparison ok: stat.out|map.out\n");
    free(dv_chip_rtk);
    free(dh_chip_rtk);
    free(dv_demo_rtk);
    free(dh_demo_rtk);
    fclose(fp); return 1;
}
static int add_eph(nav_t *nav, const eph_t *eph)
{
    eph_t *nav_eph;

    if (nav->nmax<=nav->n) {
        nav->nmax+=1024;
        if (!(nav_eph=(eph_t *)realloc(nav->eph,sizeof(eph_t)*nav->nmax))) {
            free(nav->eph);
            nav->eph=NULL;
            nav->n=nav->nmax=0;
            return 0;
        }
        nav->eph=nav_eph;
    }
    nav->eph[nav->n++]=*eph;
    return 1;
}
static int add_geph(nav_t *nav, const geph_t *geph)
{
    geph_t *nav_geph;

    if (nav->ngmax<=nav->ng) {
        nav->ngmax+=1024;
        if (!(nav_geph=(geph_t *)realloc(nav->geph,sizeof(geph_t)*nav->ngmax))) {
            free(nav->geph);
            nav->geph=NULL;
            nav->ng=nav->ngmax=0;
            return 0;
        }
        nav->geph=nav_geph;
    }
    nav->geph[nav->ng++]=*geph;
    return 1;
}
static int decode_nav_gps(const char **val, int nval, eph_t *eph)
{
    eph->sat=satid2no(val[1]);
    eph->iode=atoi(val[2]);
    eph->iodc=atoi(val[3]);
    eph->sva =atoi(val[4]);
    eph->svh =atoi(val[5]);
    eph->week=atoi(val[6]);
    eph->toe =gpst2time(eph->week,atof(val[7]));
    eph->toc =gpst2time(eph->week,atof(val[8]));
    eph->ttr =gpst2time(eph->week,atof(val[9]));
    eph->A   =atof(val[10]);
    eph->e   =atof(val[11]);
    eph->i0  =atof(val[12]);
    eph->OMG0=atof(val[13]);
    eph->omg =atof(val[14]);
    eph->M0  =atof(val[15]);
    eph->deln=atof(val[16]);
    eph->OMGd=atof(val[17]);
    eph->idot=atof(val[18]);
    eph->crc =atof(val[19]);
    eph->crs =atof(val[20]);
    eph->cuc =atof(val[21]);
    eph->cus =atof(val[22]);
    eph->cic =atof(val[23]);
    eph->cis =atof(val[24]);
    eph->toes=atof(val[25]);
    eph->fit =atof(val[26]);
    eph->f0  =atof(val[27]);
    eph->f1  =atof(val[28]);
    eph->f2  =atof(val[29]);
    eph->tgd[0]=atof(val[30]);
    eph->code=atoi(val[31]);
    eph->flag=atoi(val[32]);
    return 2;
}
static int decode_nav_gal(const char **val, int nval, eph_t *eph)
{
    eph->sat=satid2no(val[1]);
    eph->iode=atoi(val[2]);
    eph->iodc=atoi(val[3]);
    eph->sva =atoi(val[4]);
    eph->svh =atoi(val[5]);
    eph->week=atoi(val[6]);
    eph->toe =gpst2time(eph->week,atof(val[7]));
    eph->toc =gpst2time(eph->week,atof(val[8]));
    eph->ttr =gpst2time(eph->week,atof(val[9]));
    eph->A   =atof(val[10]);
    eph->e   =atof(val[11]);
    eph->i0  =atof(val[12]);
    eph->OMG0=atof(val[13]);
    eph->omg =atof(val[14]);
    eph->M0  =atof(val[15]);
    eph->deln=atof(val[16]);
    eph->OMGd=atof(val[17]);
    eph->idot=atof(val[18]);
    eph->crc =atof(val[19]);
    eph->crs =atof(val[20]);
    eph->cuc =atof(val[21]);
    eph->cus =atof(val[22]);
    eph->cic =atof(val[23]);
    eph->cis =atof(val[24]);
    eph->toes=atof(val[25]);
    eph->fit =atof(val[26]);
    eph->f0  =atof(val[27]);
    eph->f1  =atof(val[28]);
    eph->f2  =atof(val[29]);
    eph->tgd[0]=atof(val[30]);
    eph->code=atoi(val[31]);
    eph->flag=atoi(val[32]);
    return 2;
}
static int decode_nav_bds(const char **val, int nval, eph_t *eph)
{
    eph->sat=satid2no(val[1]);
    eph->iode=atoi(val[2]);
    eph->iodc=atoi(val[3]);
    eph->sva =atoi(val[4]);
    eph->svh =atoi(val[5]);
    eph->week=atoi(val[6]);
    eph->toe =bdt2time(eph->week,atof(val[7]));
    eph->toc =bdt2time(eph->week,atof(val[8]));
    eph->ttr =bdt2time(eph->week,atof(val[9]));
    eph->A   =atof(val[10]);
    eph->e   =atof(val[11]);
    eph->i0  =atof(val[12]);
    eph->OMG0=atof(val[13]);
    eph->omg =atof(val[14]);
    eph->M0  =atof(val[15]);
    eph->deln=atof(val[16]);
    eph->OMGd=atof(val[17]);
    eph->idot=atof(val[18]);
    eph->crc =atof(val[19]);
    eph->crs =atof(val[20]);
    eph->cuc =atof(val[21]);
    eph->cus =atof(val[22]);
    eph->cic =atof(val[23]);
    eph->cis =atof(val[24]);
    eph->toes=atof(val[25]);
    eph->fit =atof(val[26]);
    eph->f0  =atof(val[27]);
    eph->f1  =atof(val[28]);
    eph->f2  =atof(val[29]);
    eph->tgd[0]=atof(val[30]);
    eph->code=atoi(val[31]);
    eph->flag=atoi(val[32]);
    return 2;
}
static int decode_nav_qzs(const char **val, int nval, eph_t *eph)
{
    eph->sat=satid2no(val[1]);
    eph->iode=atoi(val[2]);
    eph->iodc=atoi(val[3]);
    eph->sva =atoi(val[4]);
    eph->svh =atoi(val[5]);
    eph->week=atoi(val[6]);
    eph->toe =gpst2time(eph->week,atof(val[7]));
    eph->toc =gpst2time(eph->week,atof(val[8]));
    eph->ttr =gpst2time(eph->week,atof(val[9]));
    eph->A   =atof(val[10]);
    eph->e   =atof(val[11]);
    eph->i0  =atof(val[12]);
    eph->OMG0=atof(val[13]);
    eph->omg =atof(val[14]);
    eph->M0  =atof(val[15]);
    eph->deln=atof(val[16]);
    eph->OMGd=atof(val[17]);
    eph->idot=atof(val[18]);
    eph->crc =atof(val[19]);
    eph->crs =atof(val[20]);
    eph->cuc =atof(val[21]);
    eph->cus =atof(val[22]);
    eph->cic =atof(val[23]);
    eph->cis =atof(val[24]);
    eph->toes=atof(val[25]);
    eph->fit =atof(val[26]);
    eph->f0  =atof(val[27]);
    eph->f1  =atof(val[28]);
    eph->f2  =atof(val[29]);
    eph->tgd[0]=atof(val[30]);
    eph->code=atoi(val[31]);
    eph->flag=atoi(val[32]);
    return 2;
}
static int decode_nav_glo(const char **val, int nval, geph_t *geph, int week)
{
    geph->sat =satid2no(val[1]);
    geph->iode=atoi(val[2]);
    geph->frq =atoi(val[3]);
    geph->svh =atoi(val[4]);
    geph->sva =atoi(val[5]);
    geph->age =atoi(val[6]);
    geph->toe =gpst2time(week,atof(val[7]));
    geph->tof =gpst2time(week,atof(val[8]));
    geph->pos[0]=atof(val[9]);
    geph->pos[1]=atof(val[10]);
    geph->pos[2]=atof(val[11]);
    geph->vel[0]=atof(val[12]);
    geph->vel[1]=atof(val[13]);
    geph->vel[2]=atof(val[14]);
    geph->acc[0]=atof(val[15]);
    geph->acc[1]=atof(val[16]);
    geph->acc[2]=atof(val[17]);
    geph->taun  =atof(val[18]);
    geph->gamn  =atof(val[19]);
    geph->dtaun =atof(val[20]);
    return 2;
}
static int decode_nav(char *buff, araw_t *raw)
{
    char *val[36],*p,*q;
    int sat,sys,n=0;
    static int gps_week=0;

    for (p=buff;*p&&n<36;p=q+1) {
        if ((q=strchr(p,','))||(q=strchr(p,'\n'))) {val[n++]=p; *q='\0';}
        else break;
    }
    sat=satid2no(val[1]);
    sys=satsys(sat,NULL);

    if (sys==SYS_GPS) {
        gps_week=atoi(val[6]);
        raw->week=gps_week;
    }
    raw->ephsat=sat;
    switch (sys) {
        case SYS_GLO: return decode_nav_glo(val,n,raw->nav.geph,gps_week);
        case SYS_GPS: return decode_nav_gps(val,n,raw->nav.eph);
        case SYS_GAL: return decode_nav_gal(val,n,raw->nav.eph);
        case SYS_CMP: return decode_nav_bds(val,n,raw->nav.eph);
        case SYS_QZS: return decode_nav_qzs(val,n,raw->nav.eph);
        default:
            return -1;
    }
}
static int decode_raw(const char *buff, araw_t *raw)
{
    int i,ret;

    for (i=0;i<strlen(buff);i++) ret=input_android(raw,buff[i]);
    return ret;
}
static int decode_raw_data(char *buff, araw_t *raw)
{
    char *pbuf;
    int type=-1;

    if      ((pbuf=strstr(buff,"Eph"))) type=0;
    else if ((pbuf=strstr(buff,"Raw"))) type=1;

    switch (type) {
        case 0: return decode_nav(pbuf,raw);
        case 1: return decode_raw(pbuf,raw);
        default:
            return -1;
    }
}
static int update_data_nav(const araw_t *raw, nav_t *nav)
{
    int sys=satsys(raw->ephsat,NULL);
    if (sys==SYS_GLO) {
        if (!add_geph(nav,raw->nav.geph)) return 0;
    }
    else {
        if (!add_eph(nav,raw->nav.eph)) return 0;
    }
    uniqnav(nav);
    return 1;
}
static int update_data_obs(const araw_t *raw, aobs_t *obs)
{
    int i;

    for (obs->n=i=0;i<raw->obs.n;i++) {
        if (raw->obs.data[i].sat<=0||raw->obs.data[i].sat>MAXSAT) continue;
        obs->data[i]=raw->obs.data[i];
        obs->n++;
    }
    return 1;
}
static int update_data(int type, const araw_t *raw, nav_t *nav, aobs_t *obs)
{
    switch (type) {
        case 1: return update_data_obs(raw,obs);
        case 2: return update_data_nav(raw,nav);
        default:
            return 0;
    }
}
static int read_phone_sol(int format, const char *file, const char *ofile, solbuf_t *solbuf)
{
    const char *file_[8]; int ret;
    FILE *fp=NULL,*fp_out;
    araw_t raw={0};
    double rb[3]={-2164248.0388,4385009.9044,4081321.3324};

    /* reset solution buffer */
    freesolbuf(solbuf);

    /* open solution file */
    if (!(fp=fopen(file,"r"))) {
        perror("open solution file error\n");
        return 0;
    }
    if (!(fp_out=fopen(ofile,"w"))) {
        fp_out=NULL;
    }
    /* read solution from file */
    /* 1.gnsslogger format */
    if (format==FMT_POS) {
        initaraw(&raw);
        while (1) {
            ret=input_androidf(&raw,fp);
            if (ret==-2) break;
            if (ret==4) {
                if (raw.sol.stat==SOLQ_CHIP) {
                    if (fp_out) {
                        outsol(fp_out,&raw.sol,rb,&solopt_default);
                        fflush(fp_out);
                    }
                    addsol(solbuf,&raw.sol);
                }
            }
        }
    }
    /* 2.NMEA format */
    else if (format==FMT_NMA) {
        file_[0]=file;
        if (!readsol(file_,1,solbuf)) {
            perror("read solution fail\n");
            fclose(fp); return 0;
        }
    }
    /* 3.rtklib solution format */
    else if (format==FMT_SOL) {
        file_[0]=file;
        if (!readsol(file_,1,solbuf)) {
            perror("read solution fail\n");
            fclose(fp); return 0;
        }
    }
    sort_solbuf(solbuf);
    if (fp_out) fclose(fp_out);
    fclose(fp);
    freearaw(&raw); return solbuf->n;
}
extern sol_t* getphonesol(gtime_t time)
{
    sol_t *sol1=NULL;
    sol_t *sol2=NULL;
    int i,j;

    for (i=0;i<solbuf_phone.n;i++) {
        if (fabs(timediff(time,solbuf_phone.data[i].time))<0.5) {
            if (i-1>=0) {
                sol1=&solbuf_phone.data[i  ];
                sol2=&solbuf_phone.data[i-1];
                if (fabs(timediff(sol1->time,sol2->time))>5.0) {
                    return sol1;
                }
                for (j=0;j<3;j++) {
                    sol1->rr[j+3]=(sol2->rr[j]-sol1->rr[j])/timediff(sol2->time,sol1->time);
                }
            }
            return sol1;
        }
    }
    return NULL;
}
extern sol_t* getnlpsol(gtime_t time)
{
    double dt;
    int i;
    for (i=0;i<solbuf_nlp.n;i++) {
        dt=timediff(time,solbuf_nlp.data[i].time);
        if (dt<0.0) {
            return &solbuf_nlp.data[i];
        }
    }
    return NULL;
}
extern sol_t* getrtksol(gtime_t time)
{
    sol_t *sol1=NULL;
    sol_t *sol2=NULL;
    int i,j;

    for (i=0;i<solbuf_rtk.n;i++) {
        if (fabs(timediff(time,solbuf_rtk.data[i].time))<0.1) {
            if (i-1>=0) {
                sol1=&solbuf_rtk.data[i  ];
                sol2=&solbuf_rtk.data[i-1];
                if (fabs(timediff(sol1->time,sol2->time))>5.0) {
                    return sol1;
                }
                for (j=0;j<3;j++) {
                    sol1->rr[j+3]=(sol2->rr[j]-sol1->rr[j])/timediff(sol2->time,sol1->time);
                }
            }
            return sol1;
        }
    }
    return NULL;
}

static void read_refsol(const char *phone_file, const char *rtk_file, const char *nlp_file,
                        const char *phone_sol_file)
{
    char *solfile_[32];
    solfile_[0]=(char*)calloc(sizeof(char),1024);
#if READ_PHONE_SOLS
    /* read phone internal solution */
    read_phone_sol(FMT_POS,phone_file,phone_sol_file,&solbuf_phone);
    if (solbuf_phone.n<=0) {
        trace(2,"no phone internal solutions\n");
    }
#endif

#if READ_NLP_SOLS
    /* read nlp solution */
    strcpy(solfile_[0],nlp_file);
    readsol(solfile_,1,&solbuf_nlp);
    if (solbuf_nlp.n<=0) {
        trace(2,"no phone internal solutions\n");
    }
#endif

#if READ_RKT_SOLS
    /* read rtk solution */
    strcpy(solfile_[0],rtk_file);
    solbuf_rtk.time=solbuf_phone.time;
    readsol(solfile_,1,&solbuf_rtk);
    if (solbuf_rtk.n<=0) {
        trace(2,"no phone internal solutions\n");
    }
#endif
}

int main(int argc, char **argv)
{
    const char *logfile="";
    const char *phonefile="";
    const char *phonesolf=".\\phone.out";
    const char *rtkfile="";
    const char *navfile="";
    const char *nlpfile="";
    const char *solfile=".\\sol.out";
    FILE *fp_log=NULL;
    FILE *fp_sol=NULL;
    FILE *fp_google_sol=NULL;
    FILE *fp_rnx=NULL;
    char buff[1024],opt_file[1024];

    double ep[6]={2021,4,15,0,0,0},doy;
    doy=time2doy(epoch2time(ep));

    araw_t raw={0};
    aobs_t obs={0};
    nav_t nav={0};
    agpos_t agp={0};
    prcopt_t opt=prcopt_default;
    solopt_t sopt=solopt_default;
    solbuf_t solbuf={0};
    int i,trace_level=3,ret;

    opt.mode=PMODE_SINGLE;

    if (!(fp_log=fopen(logfile,"r"))) {
        fprintf(stderr,"open log file fail: %s\n",logfile);
        return 0;
    }
    if (!(fp_sol=fopen(solfile,"w"))) {
        fprintf(stderr,"open sol file fail: %s\n",solfile);
        return 0;
    }
    read_refsol(phonefile,rtkfile,nlpfile,phonesolf);

    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-o")&&i+1<argc) strcpy(opt_file,argv[++i]);
        else if (!strcmp(argv[i],"-t")&&i+1<argc) trace_level=atoi(argv[++i]);
    }
    if (trace>0) {
        traceopen(".\\trace.out");
        tracelevel(trace_level);
    }
    initaraw(&raw);
    agpinit(&agp,&opt);

    gtime_t t0={0};
    uniqnav(&nav);

    obs.data=(aobsd_t*)malloc(sizeof(aobsd_t)*MAXOBS);
    obs.nmax=128;

    while (1) {
        if (feof(fp_log)) break;
        if (!fgets(buff,sizeof(buff),fp_log)) break;

        ret=decode_raw_data(buff,&raw);
        update_data(ret,&raw,&nav,&obs);

        if (ret!=1) continue;

        for (i=0;i<NFREQ;i++) {
            fixobs(obs.data[0].time,&nav,obs.data,obs.n,NULL,i,&opt);
        }
        /* rtk positioning */
        agpos(&agp,obs.data,obs.n,&nav);

        if (agp.sol.stat) {
            outsol(fp_sol,&agp.sol,agp.rb,&sopt);
            addsol(&solbuf,&agp.sol);

            udepos(agp.sol.rr);
            for (i=0;i<obs.n;i++) {
                udsatrcvdis(obs.data[i].P[0],obs.data[i].sat);
            }
        }
    }
    cmpprc(&solbuf_phone,&solbuf,&solbuf_rtk);

    outputkml(&solbuf_phone,&solbuf,&solbuf_rtk);

    freearaw(&raw);
    agpfree(&agp);
    freenav(&nav,0xFF);
    fclose(fp_log);
    fclose(fp_sol);
    fclose(fp_google_sol);
    free(obs.data);
    freesolbuf(&solbuf);
    freesolbuf(&solbuf_phone);
    freesolbuf(&solbuf_nlp);
    freesolbuf(&solbuf_rtk);
    return 1;
}