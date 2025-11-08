/*------------------------------------------------------------------------------
 * agpos.c : Android GNSS positioning
 *-----------------------------------------------------------------------------*/
#include "agpos.h"

/* SPP positioning------------------------------------------------------------*/
extern int aspplsq(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav);
extern int asppflt(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav);
extern int asppfltnx(const prcopt_t *opt);
extern int asppflt_ISB(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav);

/* initialize agpos control ---------------------------------------------------
 * initialize agpos control struct
 * args   : agpos_t  *agp    IO  agpos control/result struct
 *          prcopt_t *opt    I   positioning options (see rtklib.h)
 * return : none
 *-----------------------------------------------------------------------------*/
extern void agpinit(agpos_t *agp, const prcopt_t *opt)
{
    assat_t ssat0={0};
    sol_t sol0={{0}};
    int i;

    for (i=0;i<6;i++) agp->rb[i]=0.0;
    agp->sol=sol0;
    agp->spp.flt.tt=0.0;

    if (opt->mode==PMODE_SINGLE&&opt->dynamics) {
        agp->spp.flt.nx=asppfltnx(opt);
    }
    agp->spp.flt.x =zeros(agp->spp.flt.nx,1);
    agp->spp.flt.xp=zeros(agp->spp.flt.nx,1);
    agp->spp.flt.P =zeros(agp->spp.flt.nx,agp->spp.flt.nx);
    agp->spp.flt.Pp=zeros(agp->spp.flt.nx,agp->spp.flt.nx);

    for (i=0;i<MAXSAT;i++) {
        agp->ssat[i]=ssat0;
    }
    agp->opt=*opt;
}
/* free agpos control ---------------------------------------------------------
 * free memory for rtk control struct
 * args   : agpos_t   *agp    IO  agpos control/result struct
 * return : none
 *-----------------------------------------------------------------------------*/
extern void agpfree(agpos_t *agp)
{
    agp->spp.flt.nx=0;

    if (agp->spp.flt.xp) free(agp->spp.flt.xp); agp->spp.flt.xp=NULL;
    if (agp->spp.flt.Pp) free(agp->spp.flt.Pp); agp->spp.flt.Pp=NULL;
    if (agp->spp.flt.x ) free(agp->spp.flt.x ); agp->spp.flt.x =NULL;
    if (agp->spp.flt.P ) free(agp->spp.flt.P ); agp->spp.flt.P =NULL;
}
/* precise positioning ---------------------------------------------------------
 * input observation data and navigation message, compute rover position by
 * precise positioning
 * args   : agpos_t *agp     IO  agpos control/result struct
 *          obsd_t *obs      I   observation data for an epoch
 *                               obs[i].rcv=1:rover,2:reference
 *                               sorted by receiver and satellte
 *          int    n         I   number of observation data
 *          nav_t  *nav      I   navigation messages
 * return : status (0:no solution,1:valid solution)
 * notes  : before calling function, base station position rtk->sol.rb[] should
 *          be properly set for relative mode except for moving-baseline
 *-----------------------------------------------------------------------------*/
extern int agpos(agpos_t *agp, const aobsd_t *obs, int n, const nav_t *nav)
{
    prcopt_t *opt=&agp->opt;
    int i,nu,nr;
    obsd_t obs_[MAXOBS];

    for (i=0;i<n;i++) aobs2obs(&obs[i],&obs_[i]);

    trace(3,"agpos: %s(%10.4lf) n=%d\n",time_str(obs[0].time,3),time2gpst(obs[0].time,NULL),n);
    trace(3,"obs=\n"); traceobs(3,obs_,n);

    /* update base station position */
    for (i=0;i<6;i++) {
        agp->rb[i]=i<3?opt->rb[i]:0.0;
    }
    /* count rover/base observations */
    for (nu=0;nu   <n&&obs[nu   ].rcv==1;nu++);
    for (nr=0;nu+nr<n&&obs[nu+nr].rcv==2;nr++);

    /* rover position by single point position */
    if (!aspplsq(agp,obs,nu,nav)) {
        trace(2,"LSQ single position error\n");
    }
    /* single point position filter */
    if (opt->dynamics) {
#if FLT_ISB_MDL
        /* single point using EKF filter in ISB model */
        asppflt_ISB(rtk,obs,nu,nav);
#else
        /* single point using EKF filter */
        asppflt(agp,obs,nu,nav);
#endif
    }
    /* single point positioning */
    if (opt->mode==PMODE_SINGLE) return 1;
    return 1;
}

