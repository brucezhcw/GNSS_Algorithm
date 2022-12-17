/*------------------------------------------------------------------------------
* pntpos.c : standard positioning
*
*          Copyright (C) 2007-2020 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2010/07/28 1.0  moved from rtkcmn.c
*                           changed api:
*                               pntpos()
*                           deleted api:
*                               pntvel()
*           2011/01/12 1.1  add option to include unhealthy satellite
*                           reject duplicated observation data
*                           changed api: ionocorr()
*           2011/11/08 1.2  enable snr mask for single-mode (rtklib_2.4.1_p3)
*           2012/12/25 1.3  add variable snr mask
*           2014/05/26 1.4  support galileo and beidou
*           2015/03/19 1.5  fix bug on ionosphere correction for GLO and BDS
*           2018/10/10 1.6  support api change of satexclude()
*           2020/11/30 1.7  support NavIC/IRNSS in pntpos()
*                           no support IONOOPT_LEX option in ioncorr()
*                           improve handling of TGD correction for each system
*                           use E1-E5b for Galileo dual-freq iono-correction
*                           use API sat2freq() to get carrier frequency
*                           add output of velocity estimation error in estvel()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define MAX(x,y)    ((x)>=(y)?(x):(y))

#if 0 /* enable GPS-QZS time offset estimation */
#define NX          (4+5)       /* # of estimated parameters */
#else
#define NX          (4+4)       /* # of estimated parameters */
#endif
#define MAXITR      10          /* max number of iteration for point pos */
#define ERR_ION     5.0         /* ionospheric delay Std (m) */
#define ERR_TROP    3.0         /* tropspheric delay Std (m) */
#define ERR_SAAS    0.3         /* Saastamoinen model error Std (m) */
#define ERR_BRDCI   0.5         /* broadcast ionosphere model error factor */
#define ERR_CBIAS   0.3         /* code bias error Std (m) */
#define REL_HUMI    0.7         /* relative humidity for Saastamoinen model */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */
# define MAX_GDOP   30          /* max gdop for valid solution  */

#define VAR_POS     SQR(30.0)	/* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0)	/* initial variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)	/* initial variance of receiver acc ((m/ss)^2) */
#define NX_F		(9 + 1 + 4)	/* 位置*3 速度*3 加速度*3 接收机钟差*1 GNSS系统间群延迟*4 */

/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, const obsd_t *obs, double el, int sys, int freq)
{
    double fact=1.5,varr,snr_rover,snrweight;

	if (freq == 1) fact *= 0.9;
	else if (freq == 2) fact *= 0.85;

    switch (sys) {
        case SYS_GPS: fact *= EFACT_GPS; break;
        case SYS_GLO: fact *= EFACT_GLO; break;
        case SYS_SBS: fact *= EFACT_SBS; break;
        case SYS_CMP: fact *= EFACT_CMP; break;
        case SYS_QZS: fact *= EFACT_QZS; break;
        case SYS_IRN: fact *= EFACT_IRN; break;
        default:      fact *= EFACT_GPS; break;
    }
    if (el<MIN_EL) el=MIN_EL;
    /* var = R^2*(a^2 + (b^2/sin(el) + c^2*(10^(0.1*(snr_max-snr_rover)))) + (d*rcv_std)^2) */
    varr=SQR(opt->err[1])+SQR(opt->err[2])/SQR(sin(el))/(sin(el));
	//snr_rover = SNR_UNIT*obs->SNR[freq];
	//pow(10, (snr_rover - 50) / 30)*(1 - (30/(pow(10,(50-10)/30)-1))*(snr_rover - 50) / (50-10));
	//snrweight = pow(10, (snr_rover - 50) / 30)*(1 - 1.460255716*(snr_rover - 50) / 40);
	//varr /= (snr_rover >= 50 ? 1 : snrweight);
	//varr=0.64 + 784 * pow(2.71828,-0.142 * obs->SNR[freq]);
    varr*=SQR(opt->eratio[0]);
    if (opt->err[7]>0.0) {
        varr+=SQR(opt->err[7]*0.01*(1<<(obs->Pstd[freq]+5)));  /* 0.01*2^(n+5) m */
    }
	if (opt->ionoopt == IONOOPT_IFLC) varr *= SQR(3.0); /* iono-free */
    return SQR(fact)*varr;
}
/* doppler measurement error variance ------------------------------------*/
static double varerr_dop(const prcopt_t *opt, const obsd_t *obs, double el, int sys, int freq)
{
    double fact=0.0003,varr,snr_rover;

	if (freq == 1) fact *= 1.0;
	else if (freq == 2) fact *= 0.9;

	switch (sys) {
        case SYS_GPS: fact *= EFACT_GPS; break;
        case SYS_GLO: fact *= EFACT_GLO; break;
        case SYS_SBS: fact *= EFACT_SBS; break;
        case SYS_CMP: fact *= EFACT_CMP; break;
        case SYS_QZS: fact *= EFACT_QZS; break;
        case SYS_IRN: fact *= EFACT_IRN; break;
        default:      fact *= EFACT_GPS; break;
    }

    //if (el<MIN_EL) el=MIN_EL;
    //varr=SQR(opt->err[1])+SQR(opt->err[2])/sin(el);
	snr_rover = SNR_UNIT*obs->SNR[freq];
	varr = fact*pow(10, 0.1*MAX(opt->err[5] - snr_rover, 0));
    
    return varr;
}
/* get group delay parameter (m) ---------------------------------------------*/
static double gettgd(int sat, const nav_t *nav, int type)
{
    int i,sys=satsys(sat,NULL);
    
    if (sys==SYS_GLO) {
        for (i=0;i<nav->ng;i++) {
            if (nav->geph[i].sat==sat) break;
        }
        return (i>=nav->ng)?0.0:-nav->geph[i].dtaun*CLIGHT;
    }
    else {
        for (i=0;i<nav->n;i++) {
            if (nav->eph[i].sat==sat) break;
        }
        return (i>=nav->n)?0.0:nav->eph[i].tgd[type]*CLIGHT;
    }
}
/* test SNR mask -------------------------------------------------------------*/
static int snrmask(const obsd_t *obs, const double *azel, const prcopt_t *opt)
{
    int f2;

    if (testsnr(0,0,azel[1],obs->SNR[0]*SNR_UNIT,&opt->snrmask)) {
        return 0;
    }
    if (opt->ionoopt==IONOOPT_IFLC) {
        f2=2;
        if (testsnr(0,f2,azel[1],obs->SNR[f2]*SNR_UNIT,&opt->snrmask)) return 0;
    }
    return 1;
}
/* iono-free or "pseudo iono-free" pseudorange with code bias correction -----*/
static double prange(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt,
                     double *var)
{
    double PC,P1,P2,gamma,b1,b2;
    int sat,sys,f2;
    
    sat=obs->sat;
    sys=satsys(sat,NULL);
    P1=obs->P[0];
    f2=2;
    P2=obs->P[f2];
    *var=0.0;
    
    if (P1==0.0 && P2==0.0) return 0.0;
    
    /* P1-C1,P2-C2 DCB correction */
    if (sys==SYS_GPS||sys==SYS_GLO) {
        //if (obs->code[0]==CODE_L1C) P1+=nav->cbias[sat-1][1]; /* C1->P1 */
        //if (obs->code[1]==CODE_L2C) P2+=nav->cbias[sat-1][2]; /* C2->P2 */
    }

    if (opt->ionoopt==IONOOPT_IFLC && P1!=0.0 && P2!=0.0) { /* dual-frequency */
		PC = P2;
        if(P2-P1>20 || P2-P1<1.5)
		{
			if (obs->SNR[2]>30) return P2;
			else if (obs->SNR[0]>35) return P1;
			else return 0;
		}

        if (sys==SYS_GPS||sys==SYS_QZS) { /* L1-L2 or L1-L5 */
			gamma = f2 == 1 ? SQR(FREQL1)/(SQR(FREQL1)-SQR(FREQL2)) : SQR(FREQL1)/(SQR(FREQL1)-SQR(FREQL5));
            return P2-gamma/(P2-P1);
        }
        else if (sys==SYS_GLO) { /* G1-G2 or G1-G3 */
			gamma = f2 == 1 ? SQR(FREQ1_GLO)/(SQR(FREQ1_GLO)-SQR(FREQ2_GLO)) : SQR(FREQ1_GLO)/(SQR(FREQ1_GLO)-SQR(FREQ3_GLO));
            return P2-gamma/(P2-P1);
        }
        else if (sys==SYS_GAL) { /* E1-E5b, E1-E5a */
			gamma = f2 == 1 ? SQR(FREQL1)/(SQR(FREQL1)-SQR(FREQE5b)) : SQR(FREQL1)/(SQR(FREQL1)-SQR(FREQL5));
            if (f2==1&&getseleph(SYS_GAL)) { /* F/NAV */
                P2-=gettgd(sat,nav,0)-gettgd(sat,nav,1); /* BGD_E5aE5b */
            }
            return P2-gamma/(P2-P1);
        }
        else if (sys==SYS_CMP) { /* B1-B2 or  B1-B2a */
			gamma = f2 == 1 ? SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQ2_CMP)) : SQR(FREQ1_CMP)/(SQR(FREQ1_CMP)-SQR(FREQL5));
            return P2-gamma/(P2-P1);
        }
        else if (sys==SYS_IRN) { /* L5-S */
			gamma = SQR(FREQs) / (SQR(FREQs) - SQR(FREQL5));
            return P2-gamma/(P2-P1);
        }
    }
    else { /* single-freq (L1/E1/B1) */
        *var=SQR(ERR_CBIAS);
        
		if (obs->SNR[2]>25) PC = P2;
		else if (obs->SNR[0]>30) PC = P1;
		else return 0;

   //     if (sys==SYS_GPS||sys==SYS_QZS) { /* L1 */
   //         b1=gettgd(sat,nav,0); /* TGD (m) */
			//return PC-b1;
   //     }
   //     else if (sys==SYS_GLO) { /* G1 */
   //         gamma=SQR(FREQ1_GLO/FREQ2_GLO);
   //         b1=gettgd(sat,nav,0); /* -dtaun (m) */
   //         return PC-b1/(gamma-1.0);
   //     }
   //     else if (sys==SYS_GAL) { /* E1 */
   //         if (getseleph(SYS_GAL)) b1=gettgd(sat,nav,0); /* BGD_E1E5a */
   //         else                    b1=gettgd(sat,nav,1); /* BGD_E1E5b */
   //         return PC-b1;
   //     }
   //     else if (sys==SYS_CMP) { /* B1I/B1Cp/B1Cd */
   //         if      (obs->code[0]==CODE_L2I) b1=gettgd(sat,nav,0); /* TGD_B1I */
   //         else if (obs->code[0]==CODE_L1P) b1=gettgd(sat,nav,2); /* TGD_B1Cp */
   //         else b1=gettgd(sat,nav,2)+gettgd(sat,nav,4); /* TGD_B1Cp+ISC_B1Cd */
   //         return PC-b1;
   //     }
   //     else if (sys==SYS_IRN) { /* L5 */
   //         gamma=SQR(FREQs/FREQL5);
   //         b1=gettgd(sat,nav,0); /* TGD (m) */
   //         return PC-gamma*b1;
   //     }
    }
	return PC;
}
/* psendorange with code bias correction -------------------------------------*/
static double prange_mulfreq(const obsd_t *obs, const nav_t *nav, const prcopt_t *opt, const int k, double *var)
{
	double P1, P2, gamma, b1, b2, P, tgd;
	int sat, sys;

	sat = obs->sat;
	sys = satsys(sat, NULL);
	P1 = obs->P[0];
	P2 = obs->P[2];
	*var = 0.0;

	if ((P1 == 0.0 && P2 == 0.0) || (opt->ionoopt == IONOOPT_IFLC && (P1 == 0.0 || P2 == 0.0)))
		return 0.0;

	/* P1-C1,P2-C2 DCB correction */
	if (sys == SYS_GPS || sys == SYS_GLO)
	{
		if (obs->code[0] == CODE_L1C)
			P1 += nav->cbias[sat - 1][1]; /* C1->P1 */
		if (obs->code[2] == CODE_L2C)
			P2 += nav->cbias[sat - 1][2]; /* C2->P2 */
	}
	if (opt->ionoopt == IONOOPT_IFLC)
	{ /* dual-frequency */

		if (sys == SYS_GPS || sys == SYS_QZS)
		{ /* L1-L2,G1-G2 */
			//gamma = SQR(FREQL1 / FREQL2);
			gamma = SQR(FREQL1 / FREQL5);
			return (P2 - gamma * P1) / (1.0 - gamma);
		}
		else if (sys == SYS_GLO)
		{ /* G1-G2 */
			gamma = SQR(FREQ1_GLO / FREQ2_GLO);
			return (P2 - gamma * P1) / (1.0 - gamma);
		}
		else if (sys == SYS_GAL)
		{ /* E1-E5b */
			//gamma = SQR(FREQL1 / FREQE5b);
			gamma = SQR(FREQL1 / FREQL5);
			//if (getseleph(SYS_GAL))
			{                                                    /* F/NAV */
				//P2 -= gettgd(sat, nav, 0) - gettgd(sat, nav, 1); /* BGD_E5aE5b */
			}
			return (P2 - gamma * P1) / (1.0 - gamma);
		}
		else if (sys == SYS_CMP)
		{ /* B1-B2 */
			//gamma = SQR(((obs->code[0] == CODE_L2I) ? FREQ1_CMP : FREQL1) / FREQ2_CMP);
			gamma = SQR(((obs->code[0] == CODE_L2I) ? FREQ1_CMP : FREQL1) / FREQL5);
			if (obs->code[0] == CODE_L2I)
				b1 = gettgd(sat, nav, 0); /* TGD_B1I */
			else if (obs->code[0] == CODE_L1P)
				b1 = gettgd(sat, nav, 2); /* TGD_B1Cp */
			else
				b1 = gettgd(sat, nav, 2) + gettgd(sat, nav, 4); /* TGD_B1Cp+ISC_B1Cd */
			//b2 = gettgd(sat, nav, 1);                           /* TGD_B2I/B2bI (m) */
			b2 = 0.0;
			return ((P2 - gamma * P1) - (b2 - gamma * b1)) / (1.0 - gamma);
		}
		else if (sys == SYS_IRN)
		{ /* L5-S */
			gamma = SQR(FREQL5 / FREQs);
			return (P2 - gamma * P1) / (1.0 - gamma);
		}
	}
	else
	{ /* single-freq */
		if (k < 0 || k >= (NFREQ + NEXOBS))
		{
			return 0.0;
		}
		P = obs->P[k];
		*var = SQR(ERR_CBIAS);
		tgd = 0.0;
		if (sys == SYS_GPS || sys == SYS_QZS)
		{
			gamma = SQR(FREQL1 / FREQL2);
			b1 = gettgd(sat, nav, 0); /* TGD (m) */
			switch (obs->code[k])
			{
				/*for L1C, we need to calibrate the DCB between P1 and C1. but here we ignore it*/
			case CODE_L1C:
			case CODE_L1P:
				tgd = b1;
				break;
				/*for L2P, calibrate tgd according the paper*/
			case CODE_L2P:
			case CODE_L2W:
				tgd = b1 * gamma;
				break;
				/*no L5 tgd info in broadcast nav, give the default value 0*/
			case CODE_L5Q:
				tgd = 0.0;
				break;
			}
		}
		else if (sys == SYS_GLO)
		{   /*for GLO, the logic may has error*/
			gamma = SQR(FREQ1_GLO / FREQ2_GLO);
			b1 = gettgd(sat, nav, 0); /* -dtaun (m) */
			b1 = b1 / (gamma - 1.0);
			switch (obs->code[k])
			{
			case CODE_L1C:
			case CODE_L1P:
				tgd = b1;
				break;
			case CODE_L2C:
			case CODE_L2P:
				tgd = b1 * gamma;
				break;
			}
		}

		else if (sys == SYS_GAL)
		{   /*for GAL sys, there is two ephemerises*/
			if (getseleph(SYS_GAL))
				b1 = gettgd(sat, nav, 0); /* BGD_E1E5a */
			else
				b1 = gettgd(sat, nav, 1); /* BGD_E1E5b */
			switch (obs->code[k])
			{
			case CODE_L1C:
				tgd = b1;
				break;
			case CODE_L5Q:
				gamma = SQR(FREQL1 / FREQL5);
				tgd = b1 * gamma;
				break;
			case CODE_L7Q:
				gamma = SQR(FREQL1 / FREQE5b);
				tgd = b1 * gamma;
				break;
			}
		}
		else if (sys == SYS_CMP)
		{
			switch (obs->code[k])
			{
			case CODE_L2I:
				tgd = gettgd(sat, nav, 0);
				break;
			case CODE_L7I:
				tgd = gettgd(sat, nav, 1);
				break;
				/*for BDS sys, the reference frequency is 6I. that's why the correction in 6I is zeror*/
			case CODE_L6I:
				tgd = 0.0;
				break;
			case CODE_L1P:
				tgd = gettgd(sat, nav, 2);
				break;
			case CODE_L5P:
				tgd = gettgd(sat, nav, 3);
				break;
			}
		}
		return P - tgd;
	}
	return P1;
}
/* ionospheric correction ------------------------------------------------------
* compute ionospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          int    sat       I   satellite number
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    ionoopt   I   ionospheric correction option (IONOOPT_???)
*          double *ion      O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric delay (L1) variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int ionocorr(gtime_t time, const nav_t *nav, int sat, const double *pos,
                    const double *azel, int ionoopt, double *ion, double *var)
{
    int err=0;

    trace(4,"ionocorr: time=%s opt=%d sat=%2d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),ionoopt,sat,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* SBAS ionosphere model */
    if (ionoopt==IONOOPT_SBAS) {
        if (sbsioncorr(time,nav,pos,azel,ion,var)) return 1;
        err=1;
    }
    /* IONEX TEC model */
    if (ionoopt==IONOOPT_TEC) {
        if (iontec(time,nav,pos,azel,1,ion,var)) return 1;
        err=1;
    }
    /* QZSS broadcast ionosphere model */
    if (ionoopt==IONOOPT_QZS&&norm(nav->ion_qzs,8)>0.0) {
        *ion=ionmodel(time,nav->ion_qzs,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    /* GPS broadcast ionosphere model */
    if (ionoopt==IONOOPT_BRDC||err==1) {
        *ion=ionmodel(time,nav->ion_gps,pos,azel);
        *var=SQR(*ion*ERR_BRDCI);
        return 1;
    }
    *ion=0.0;
    *var=ionoopt==IONOOPT_OFF?SQR(ERR_ION):0.0;
    return 1;
}
/* tropospheric correction -----------------------------------------------------
* compute tropospheric correction
* args   : gtime_t time     I   time
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad|m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    tropopt   I   tropospheric correction option (TROPOPT_???)
*          double *trp      O   tropospheric delay (m)
*          double *var      O   tropospheric delay variance (m^2)
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int tropcorr(gtime_t time, const nav_t *nav, const double *pos,
                    const double *azel, int tropopt, double *trp, double *var)
{
    trace(4,"tropcorr: time=%s opt=%d pos=%.3f %.3f azel=%.3f %.3f\n",
          time_str(time,3),tropopt,pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,
          azel[1]*R2D);
    
    /* Saastamoinen model */
    if (tropopt==TROPOPT_SAAS||tropopt==TROPOPT_EST||tropopt==TROPOPT_ESTG) {
        *trp=tropmodel(time,pos,azel,REL_HUMI);
        *var=SQR(ERR_SAAS/(sin(azel[1])+0.1));
        return 1;
    }
    /* SBAS (MOPS) troposphere model */
    if (tropopt==TROPOPT_SBAS) {
        *trp=sbstropcorr(time,pos,azel,var);
        return 1;
    }
    /* no correction */
    *trp=0.0;
    *var=tropopt==TROPOPT_OFF?SQR(ERR_TROP):0.0;
    return 1;
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode(int iter, const obsd_t *obs, int n, const double *rs,
                   const double *dts, const double *vare, const int *svh,
                   const nav_t *nav, const double *x, const prcopt_t *opt,
                   const ssat_t *ssat, double *v, double *H, double *var,
                   double *azel, int *vsat, double *resp, int *ns)
{
    gtime_t time;
    double r,freq,dion=0.0,dtrp=0.0,vmeas,vion=0.0,vtrp=0.0,rr[3],pos[3],dtr,e[3],P;
    int i,j,nv=0,sat,sys,mask[NX-3]={0};

    for (i=0;i<3;i++) rr[i]=x[i];
    dtr=x[3];
    
    ecef2pos(rr,pos);
    trace(3,"rescode: rr=%.3f %.3f %.3f\n",rr[0], rr[1], rr[2]);
    
    for (i=*ns=0;i<n&&i<MAXOBS;i++) {
        vsat[i]=0; azel[i*2]=azel[1+i*2]=resp[i]=0.0;
        time=obs[i].time;
        sat=obs[i].sat;
        if (!(sys=satsys(sat,NULL))) continue;
        
        /* reject duplicated observation data */
        if (i<n-1&&i<MAXOBS-1&&sat==obs[i+1].sat) {
            trace(2,"duplicated obs data %s sat=%d\n",time_str(time,3),sat);
            i++;
            continue;
        }
        /* excluded satellite? */
        if (satexclude(sat,vare[i],svh[i],opt)) continue;
        
        /* geometric distance and elevation mask*/
        if ((r=geodist(rs+i*6,rr,e))<=0.0) continue;
        if (satazel(pos,e,azel+i*2)<opt->elmin) continue;
        
        if (iter>0) {
            /* test SNR mask */
            if (!snrmask(obs+i,azel+i*2,opt)) continue;
        
            /* ionospheric correction */
            if (!ionocorr(time,nav,sat,pos,azel+i*2,opt->ionoopt,&dion,&vion)) {
                continue;
            }
            if ((freq=sat2freq(sat,obs[i].code[0],nav))==0.0) continue;
            dion*=SQR(FREQL1/freq);
            vion*=SQR(FREQL1/freq);
        
            /* tropospheric correction */
            if (!tropcorr(time,nav,pos,azel+i*2,opt->tropopt,&dtrp,&vtrp)) {
                continue;
            }
        }
        /* psendorange with code bias correction */
        if ((P=prange(obs+i,nav,opt,&vmeas))==0.0) continue;
        
        /* pseudorange residual */
        v[nv]=P-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
        trace(4,"sat=%d: v=%.3f P=%.3f r=%.3f dtr=%.6f dts=%.6f dion=%.3f dtrp=%.3f\n",
            sat,v[nv],P,r,dtr,dts[i*2],dion,dtrp);
        
        /* design matrix */
        for (j=0;j<NX;j++) {
            H[j+nv*NX]=j<3?-e[j]:(j==3?1.0:0.0);
        }
        /* time system offset and receiver bias correction */
        if      (sys==SYS_GLO) {v[nv]-=x[4]; H[4+nv*NX]=1.0; mask[1]=1;}
        else if (sys==SYS_GAL) {v[nv]-=x[5]; H[5+nv*NX]=1.0; mask[2]=1;}
        else if (sys==SYS_CMP) {v[nv]-=x[6]; H[6+nv*NX]=1.0; mask[3]=1;}
        else if (sys==SYS_IRN) {v[nv]-=x[7]; H[7+nv*NX]=1.0; mask[4]=1;}
#if 0 /* enable QZS-GPS time offset estimation */
        else if (sys==SYS_QZS) {v[nv]-=x[8]; H[8+nv*NX]=1.0; mask[5]=1;}
#endif
        else mask[0]=1;

        vsat[i]=1; resp[i]=v[nv]; (*ns)++;
        
        /* variance of pseudorange error */
        var[nv++]=varerr(opt,&obs[i],azel[1+i*2],sys,0)+vare[i]+vmeas+vion+vtrp;
        trace(4,"sat=%2d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n",obs[i].sat,
              azel[i*2]*R2D,azel[1+i*2]*R2D,resp[i],sqrt(var[nv-1]));
    }
    /* constraint to avoid rank-deficient */
    for (i=0;i<NX-3;i++) {
        if (mask[i]) continue;
        v[nv]=0.0;
        for (j=0;j<NX;j++) H[j+nv*NX]=j==i+3?1.0:0.0;
        var[nv++]=0.01;
    }
    return nv;
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode_mulfreq(int iter, const obsd_t *obs, int n, const double *rs,
	const double *dts, const double *vare, const int *svh,
	const nav_t *nav, const double *x, const prcopt_t *opt,
	double *v, double *H, double *var, double *azel, int *vsat,
	double *resp, int *ns)
{
	gtime_t time;
	double r, freq, dion_ref = 0.0, dion = 0.0, dtrp = 0.0, vmeas, vion = 0.0, vtrp = 0.0, rr[3], pos[3], dtr, e[3], P;
	int i, j, nv = 0, sat, sys, mask[NX - 3] = { 0 };
	int freq_idx, nf;

	trace(3, "rescode_mulfreq : n=%d\n", n);

	nf = (opt->spp_mode == SPP_MODE_LX) ? opt->nf : 1;
	for (i = 0; i < 3; i++)
		rr[i] = x[i];
	dtr = x[3];

	ecef2pos(rr, pos);
	for (i = *ns = 0; i < n && i < MAXOBS; i++)
	{
		vsat[i] = 0;
		azel[i * 2] = azel[1 + i * 2] = resp[i] = 0.0;
		time = obs[i].time;
		sat = obs[i].sat;
		if (!(sys = satsys(sat, NULL))) continue;
		/* reject duplicated observation data */
		if (i < n - 1 && i < MAXOBS - 1 && sat == obs[i + 1].sat)
		{
			trace(2, "duplicated obs data %s sat=%d\n", time_str(time, 3), sat);
			i++;
			continue;
		}
		/* excluded satellite? */
		if (satexclude(sat, vare[i], svh[i], opt)) continue;
		/* geometric distance */
		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0) continue;

		if (iter > 0)
		{
			/* test elevation mask */
			if (satazel(pos, e, azel + i * 2) < opt->elmin) continue;
			/* test SNR mask */
			if (!snrmask(obs + i, azel + i * 2, opt)) continue;
			/* ionospheric correction */
			if (!ionocorr(time, nav, sat, pos, azel + i * 2, opt->ionoopt, &dion_ref, &vion)) continue;
			/* tropospheric correction */
			if (!tropcorr(time, nav, pos, azel + i * 2, opt->tropopt, &dtrp, &vtrp)) continue;
		}

		for (freq_idx = 0; freq_idx < nf; freq_idx++)
		{
			dion = dion_ref;
			if (iter > 0)
			{
				if ((freq = sat2freq(sat, obs[i].code[freq_idx], nav)) == 0.0)
					continue;
				dion *= SQR(FREQL1 / freq);
				vion *= SQR(FREQL1 / freq);
			}

			if (opt->ionoopt == IONOOPT_IFLC){
				dion = 0.0;
				vion = 0.0;
			}

			/* psendorange with code bias correction */
			if ((P = prange_mulfreq(obs + i, nav, opt, freq_idx, &vmeas)) == 0.0) continue;
			/* pseudorange residual */
			v[nv] = P - (r + dtr - CLIGHT * dts[i * 2] + dion + dtrp);
			/* design matrix */
			for (j = 0; j < NX; j++)
			{
				H[j + nv * NX] = j < 3 ? -e[j] : (j == 3 ? 1.0 : 0.0);
			}
			/* time system offset and receiver bias correction */
			if (sys == SYS_GLO)
			{
				v[nv] -= x[4];
				H[4 + nv * NX] = 1.0;
				mask[1] = 1;
			}
			else if (sys == SYS_GAL)
			{
				v[nv] -= x[5];
				H[5 + nv * NX] = 1.0;
				mask[2] = 1;
			}
			else if (sys == SYS_CMP)
			{
				v[nv] -= x[6];
				H[6 + nv * NX] = 1.0;
				mask[3] = 1;
			}
			else if (sys == SYS_IRN)
			{
				v[nv] -= x[7];
				H[7 + nv * NX] = 1.0;
				mask[4] = 1;
			}
#if 0 /* enable QZS-GPS time offset estimation */
			else if (sys == SYS_QZS) { v[nv] -= x[8]; H[8 + nv*NX] = 1.0; mask[5] = 1; }
#endif
			else
				mask[0] = 1;

			vsat[i] = 1;
			resp[i] = v[nv];

			/* variance of pseudorange error */
			var[nv++] = varerr(opt, &obs[i], azel[1 + i * 2], sys, freq_idx) + vare[i] + vmeas + vion + vtrp;

			trace(4, "sat=%3d azel=%5.1f %4.1f res=%7.3f sig=%5.3f\n", obs[i].sat,
				azel[i * 2] * R2D, azel[1 + i * 2] * R2D, resp[i], sqrt(var[nv - 1]));
		}
		(*ns)++;
	}
	/* constraint to avoid rank-deficient */
	for (i = 0; i < NX - 3; i++)
	{
		if (mask[i])
			continue;
		v[nv] = 0.0;
		for (j = 0; j < NX; j++)
			H[j + nv * NX] = j == i + 3 ? 1.0 : 0.0;
		var[nv++] = 0.01;
	}
	return nv;
}
/* pseudorange residuals -----------------------------------------------------*/
static int rescode_mulfreq_ekf(int iter, const obsd_t *obs, int n, const double *rs,
	const double *dts, const double *vare, const int *svh,
	const nav_t *nav, const double *x, const prcopt_t *opt,
	double *v, double *H, double *var, double *azel, int *vsat,
	double *resp, int *ns, ssat_t *ssat)
{
	gtime_t time;
	double r, freq, dion_ref = 0.0, dion = 0.0, dtrp = 0.0, vmeas, vion = 0.0, vtrp = 0.0, rr[3], pos[3], dtr = 0.0, e[3], P;
	int i, j, nv = 0, sat, sys, mask[NX_F - 9] = { 0 };
	int freq_idx, nf;

	trace(3, "rescode_mulfreq_ekf : n=%d\n", n);

	nf = (opt->spp_mode == SPP_MODE_LX) ? opt->nf : 1;
	for (i = 0; i < 3; i++)
		rr[i] = x[i];
	dtr = x[9];

	ecef2pos(rr, pos);
	for (i = *ns = 0; i < n && i < MAXOBS; i++)
	{
		vsat[i] = 0;
		azel[i * 2] = azel[1 + i * 2] = resp[i] = 0.0;
		time = obs[i].time;
		sat = obs[i].sat;
		if (!(sys = satsys(sat, NULL))) continue;
		/* reject duplicated observation data */
		if (i < n - 1 && i < MAXOBS - 1 && sat == obs[i + 1].sat)
		{
			trace(2, "duplicated obs data %s sat=%d\n", time_str(time, 3), sat);
			i++;
			continue;
		}
		/* excluded satellite? */
		if (satexclude(sat, vare[i], svh[i], opt)) continue;
		/* geometric distance */
		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0) continue;

		if (iter > 0)
		{
			/* test elevation mask */
			if (satazel(pos, e, azel + i * 2) < opt->elmin) continue;
			/* test SNR mask */
			if (!snrmask(obs + i, azel + i * 2, opt)) continue;
			/* ionospheric correction */
			if (!ionocorr(time, nav, sat, pos, azel + i * 2, opt->ionoopt, &dion_ref, &vion)) continue;
			/* tropospheric correction */
			if (!tropcorr(time, nav, pos, azel + i * 2, opt->tropopt, &dtrp, &vtrp)) continue;
		}

		for (freq_idx = 0; freq_idx < nf; freq_idx++)
		{
			dion = dion_ref;
			if (iter > 0)
			{
				if ((freq = sat2freq(sat, obs[i].code[freq_idx], nav)) == 0.0)
					continue;
				dion *= SQR(FREQL1 / freq);
				vion *= SQR(FREQL1 / freq);
			}

			if (opt->ionoopt == IONOOPT_IFLC){
				dion = 0.0;
				vion = 0.0;
			}

			/* psendorange with code bias correction */
			if ((P = prange_mulfreq(obs + i, nav, opt, freq_idx, &vmeas)) == 0.0) continue;

			//if (ssat[obs[i].sat-1].lock_P[0][freq_idx]<10) continue;
			/* pseudorange residual */
			v[nv] = P - (r + dtr - CLIGHT * dts[i * 2] + dion + dtrp);

			/* design matrix */
			for (j = 0; j < NX_F; j++)
			{
				H[j + nv * NX_F] = j < 3 ? -e[j] : (j == 9 ? 1.0 : 0.0);
			}
			/* time system offset and receiver bias correction */
			if (sys == SYS_GLO)
			{
				v[nv] -= x[10];
				H[10 + nv * NX_F] = 1.0;
				mask[1] = 1;
			}
			else if (sys == SYS_GAL)
			{
				v[nv] -= x[11];
				H[11 + nv * NX_F] = 1.0;
				mask[2] = 1;
			}
			else if (sys == SYS_CMP)
			{
				v[nv] -= x[12];
				H[12 + nv * NX_F] = 1.0;
				mask[3] = 1;
			}
			else if (sys == SYS_IRN)
			{
				v[nv] -= x[13];
				H[13 + nv * NX_F] = 1.0;
				mask[4] = 1;
			}
			else
				mask[0] = 1;

			vsat[i] = 1;
			resp[i] = v[nv];

			/* variance of pseudorange error */
			var[nv] = varerr(opt, &obs[i], azel[1 + i * 2], sys, freq_idx);
			if (v[nv]<-20 || v[nv]>20) var[nv] *= SQR(30);
			else if (fabs(v[nv])>10)
			{
				var[nv] *= SQR(fabs(v[nv] / 10));
			}
			nv++;

			trace(3, "sat=%3d azel=%5.1f %4.1f f=%1d snr=%2d res=%7.3f sig=%5.3f\n", obs[i].sat,
				azel[i * 2] * R2D, azel[1 + i * 2] * R2D, freq_idx+1, obs[i].SNR[freq_idx], resp[i], sqrt(var[nv - 1]));
		}
		(*ns)++;
	}
	/* constraint to avoid rank-deficient */
	for (i = 0; i < NX_F - 9; i++)
	{
		if (mask[i]) continue;
		v[nv] = 0.0;
		for (j = 0; j < NX_F; j++) H[j + nv * NX_F] = j == i + 9 ? 1.0 : 0.0;
		var[nv++] = SQR(0.01);
	}
	return nv;
}
static int resdop_mulfreq_ekf(const obsd_t *obs, const int n, const double *rs, const double *dts,
                                 const nav_t *nav, const double *rr, const prcopt_t *opt, const double *x,
                                 const double *azel, const int *vsat, double *v, double *var, double *H)
{
    double freq, rate, pos[3], E[9], a[3], e[3], vs[3], cosel, sig, e_ref[3] = {0.0}, v_ref = 0.0, ele_ref=0.0;
    int i, j, nv = 0, sys, freq_idx, SNR_ref=0, ind_ref=-1, freq_ref=-1;;

    trace(3, "resdop_mulfreq_ekf  : n=%d\n", n);

    ecef2pos(rr, pos);
    xyz2enu(pos, E);

	for (i = 0; i < n && i < MAXOBS; i++)
    {

        if (!(sys = satsys(obs[i].sat, NULL)) || !vsat[i] || norm(rs + 3 + i * 6, 3) <= 0.0) continue;

        for (freq_idx = 0; freq_idx < opt->nf; freq_idx++)
        {
            freq = sat2freq(obs[i].sat, obs[i].code[freq_idx], nav);

            if (obs[i].D[freq_idx] == 0.0 || freq == 0.0) continue;

			if (ind_ref == -1 || (obs[i].SNR[freq_idx]>=SNR_ref && azel[1 + i * 2] * R2D>ele_ref) ||
								 (obs[i].SNR[freq_idx]>SNR_ref && azel[1 + i * 2] * R2D>ele_ref-5))
            {
				ele_ref = azel[1 + i * 2] * R2D;
				SNR_ref = obs[i].SNR[freq_idx];
				ind_ref = i;
				freq_ref = freq_idx;
            }
        }
    }

	if (ind_ref == -1) return 0;

	/* LOS (line-of-sight) vector in ECEF */
	cosel = cos(azel[1 + ind_ref * 2]);
	a[0] = sin(azel[ind_ref * 2]) * cosel;
	a[1] = cos(azel[ind_ref * 2]) * cosel;
	a[2] = sin(azel[1 + ind_ref * 2]);
	matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);

	/* satellite velocity relative to receiver in ECEF */
	for (j = 0; j < 3; j++)
	{
		vs[j] = rs[j + 3 + ind_ref * 6] - x[j];
	}

	/* range rate with earth rotation correction */
	rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + ind_ref * 6] * rr[0] + rs[1 + ind_ref * 6] * x[0] - rs[3 + ind_ref * 6] * rr[1] - rs[ind_ref * 6] * x[1]);

	freq = sat2freq(obs[ind_ref].sat, obs[ind_ref].code[freq_ref], nav);

	/* range rate residual (m/s). */
	v_ref = (-obs[ind_ref].D[freq_ref] * CLIGHT / freq - (rate - CLIGHT * dts[1 + ind_ref * 2]));
	for (j = 0; j < 3; j++)
	{
		e_ref[j] = -e[j];
	}
	trace(3, "reference sat=%3d azel=%5.1f %4.1f f=%1d SNR=%2d\n", obs[ind_ref].sat, azel[ind_ref * 2] * R2D, ele_ref, freq_ref + 1, SNR_ref);

    for (i = 0; i < n && i < MAXOBS; i++)
    {

        if (!(sys = satsys(obs[i].sat, NULL)) || !vsat[i] || norm(rs + 3 + i * 6, 3) <= 0.0) continue;
        /* LOS (line-of-sight) vector in ECEF */
        cosel = cos(azel[1 + i * 2]);
        a[0] = sin(azel[i * 2]) * cosel;
        a[1] = cos(azel[i * 2]) * cosel;
        a[2] = sin(azel[1 + i * 2]);
        matmul("TN", 3, 1, 3, 1.0, E, a, 0.0, e);
        /* satellite velocity relative to receiver in ECEF */
        for (j = 0; j < 3; j++)
        {
            vs[j] = rs[j + 3 + i * 6] - x[j];
        }

        /* range rate with earth rotation correction */
        rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + i * 6] * rr[0] + rs[1 + i * 6] * x[0] - rs[3 + i * 6] * rr[1] - rs[i * 6] * x[1]);

        for (freq_idx = 0; freq_idx < opt->nf; freq_idx++)
        {
            freq = sat2freq(obs[i].sat, obs[i].code[freq_idx], nav);

            if (obs[i].D[freq_idx] == 0.0 || freq == 0.0) continue;

			if (i == ind_ref && freq_idx == freq_ref) continue;

            /* range rate residual (m/s). */
            v[nv] = (-obs[i].D[freq_idx] * CLIGHT / freq - (rate - CLIGHT * dts[1 + i * 2])) - v_ref;

            /* variance of doppler error */
			var[nv] = varerr_dop(opt, &obs[i], azel[1 + i * 2], sys, freq_idx);
			if (v[nv]<-1.0 || v[nv]>1.0) var[nv] *= SQR(30);
			else if (fabs(v[nv])>0.5)
			{
				var[nv] *= SQR(fabs(v[nv]) / 0.5);
			}

            /* design matrix */
            for (j = 0; j < 3; j++)
            {
                H[j + 3 + nv * NX_F] = -e[j] - e_ref[j];
            }
            
            nv++;

			trace(3, "sat=%3d azel=%5.1f %4.1f f=%1d snr=%2d residual=%7.3f sig=%5.3f\n", obs[i].sat, azel[i * 2] * R2D, azel[1 + i * 2] * R2D,
				freq_idx + 1, obs[i].SNR[freq_idx], v[nv - 1], sqrt(var[nv - 1]));
        }
    }
    return nv;
}
/* validate solution ---------------------------------------------------------*/
static int valsol(const double *azel, const int *vsat, int n,
                  const prcopt_t *opt, const double *v, int nv, int nx,
                  char *msg)
{
    double azels[MAXOBS*2],dop[4],vv;
    int i,ns;
    
    trace(3,"valsol  : n=%d nv=%d\n",n,nv);
    
    /* Chi-square validation of residuals */
    vv=dot(v,v,nv);
    if (nv>nx&&vv>chisqr[nv-nx-1]) {
        sprintf(msg,"Warning: large chi-square error nv=%d vv=%.1f cs=%.1f",nv,vv,chisqr[nv-nx-1]);
        /* return 0; */ /* threshold too strict for all use cases, report error but continue on */
    }
    /* large GDOP check */
    for (i=ns=0;i<n;i++) {
        if (!vsat[i]) continue;
        azels[  ns*2]=azel[  i*2];
        azels[1+ns*2]=azel[1+i*2];
        ns++;
    }
    dops(ns,azels,opt->elmin,dop);
    if (dop[0]<=0.0||dop[0]>MAX_GDOP) {
        sprintf(msg,"gdop error nv=%d gdop=%.1f",nv,dop[0]);
        return 0;
    }
    return 1;
}
/* estimate receiver position ------------------------------------------------*/
static int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const double *vare, const int *svh, const nav_t *nav,
                  const prcopt_t *opt, const ssat_t *ssat, sol_t *sol, double *azel,
                  int *vsat, double *resp, char *msg)
{
    double x[NX]={0},dx[NX],Q[NX*NX],*v,*H,*var,sig;
    int i,j,k,info,stat,nv,ns;
    
    trace(3,"estpos  : n=%d\n",n);
    
    v=mat(n+4,1); H=mat(NX,n+4); var=mat(n+4,1);
    
    for (i=0;i<3;i++) x[i]=sol->rr[i];

    for (i=0;i<MAXITR;i++) {

        /* pseudorange residuals (m) */
        nv=rescode(i,obs,n,rs,dts,vare,svh,nav,x,opt,ssat,v,H,var,azel,vsat,resp,
                   &ns);
        
        if (nv<NX) {
            sprintf(msg,"lack of valid sats ns=%d",nv);
            break;
        }
        /* weight by variance (lsq uses sqrt of weight */
        for (j=0;j<nv;j++) {
            sig=sqrt(var[j]);
            v[j]/=sig;
            for (k=0;k<NX;k++) H[k+j*NX]/=sig;
        }
        /* least square estimation */
        if ((info=lsq(H,v,NX,nv,dx,Q))) {
            sprintf(msg,"lsq error info=%d",info);
            break;
        }
        for (j=0;j<NX;j++) {
            x[j]+=dx[j];
        }
        if (norm(dx,NX)<1E-4) {
            sol->type=0;
            sol->time=timeadd(obs[0].time,-x[3]/CLIGHT);
            sol->dtr[0]=x[3]/CLIGHT; /* receiver clock bias (s) */
            sol->dtr[1]=x[4]/CLIGHT; /* GLO-GPS time offset (s) */
            sol->dtr[2]=x[5]/CLIGHT; /* GAL-GPS time offset (s) */
            sol->dtr[3]=x[6]/CLIGHT; /* BDS-GPS time offset (s) */
            sol->dtr[4]=x[7]/CLIGHT; /* IRN-GPS time offset (s) */
            for (j=0;j<6;j++) sol->rr[j]=j<3?x[j]:0.0;
            for (j=0;j<3;j++) sol->qr[j]=(float)Q[j+j*NX];
            sol->qr[3]=(float)Q[1];    /* cov xy */
            sol->qr[4]=(float)Q[2+NX]; /* cov yz */
            sol->qr[5]=(float)Q[2];    /* cov zx */
            sol->ns=(uint8_t)ns;
            sol->age=sol->ratio=0.0;
            
            /* validate solution */
            if ((stat=valsol(azel,vsat,n,opt,v,nv,NX,msg))) {
                sol->stat=opt->sateph==EPHOPT_SBAS?SOLQ_SBAS:SOLQ_SINGLE;
            }
            free(v); free(H); free(var);
            return stat;
        }
    }
    if (i>=MAXITR) sprintf(msg,"iteration divergent i=%d",i);
    
    free(v); free(H); free(var);
    return 0;
}
static void init_spp(rtk_t*rtk){
	int i = 0, j = 0;
	if (!rtk->x_spp)
	{
		rtk->x_spp = zeros(NX_F, 1);
	}
	if (!rtk->P_spp)
	{
		rtk->P_spp = zeros(NX_F, NX_F);
	}
	for (i = 0; i < NX_F; i++)
	{
		rtk->x_spp[i] = 0.0;
		for (j = 0; j < NX_F; j++)
		{
			rtk->P_spp[i * NX_F + j] = 0.0;
		}
	}

	for (i = 0; i < MAXSAT; i++)
	{
		for (j = 0; j < rtk->opt.nf; j++)
		{
			rtk->ssat[i].lock_P[0][j] = 0;
			rtk->ssat[i].lock_P[1][j] = 0;
		}
	}
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
	int j;
	rtk->x_spp[i] = xi;
	for (j = 0; j < NX_F; j++)
	{
		rtk->P_spp[i + j * NX_F] = rtk->P_spp[j + i * NX_F] = i == j ? var : 0.0;
	}
}
/* temporal update of position/velocity/acceleration -------------------------*/
static void udpos_spp(rtk_t *rtk, double tt)
{
	double *F, *P, *FP, *x, *xp, pos[3], Q[9] = { 0 }, Qv[9], var = 0.0;
	int i, j, *ix, nx;

	trace(3, "udpos_spp   : tt=%.3f\n", tt);

	/* initialize position for first epoch */
	if (norm(rtk->x_spp, 3) <= 0.0)
	{
		for (i = 0; i < 3; i++)
			initx(rtk, rtk->sol.rr[i], VAR_POS, i);
		for (i = 3; i < 6; i++)
			initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
		for (i = 6; i < 9; i++)
			initx(rtk, 1E-6, VAR_ACC, i);
	}

	/* check variance of estimated postion */
	for (i = 0; i < 3; i++)
		var += rtk->P_spp[i + i * NX_F];
	var /= 3.0;
	trace(3, "EKF pos var: %.4f\n", var);
	if (var > VAR_POS)
	{
		/* reset position with large variance */
		for (i = 0; i < 3; i++)
			initx(rtk, rtk->sol.rr[i], VAR_POS, i);
		for (i = 3; i < 6; i++)
			initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
		for (i = 6; i < 9; i++)
			initx(rtk, 1E-6, VAR_ACC, i);
		trace(2, "reset EKF position due to large variance: var=%.3f\n", var);
		return;
	}
	/* generate valid state index */
	ix = imat(NX_F, 1);
	for (i = nx = 0; i < NX_F; i++)
	{
		if (i < 9) ix[nx++] = i;
	}
	if (nx < 9)
	{
		free(ix);
		return;
	}
	/* state transition of position/velocity/acceleration */
	F = eye(nx);
	P = mat(nx, nx);
	FP = mat(nx, nx);
	x = mat(nx, 1);
	xp = mat(nx, 1);

	for (i = 0; i < 6; i++)
	{
		F[i + (i + 3) * nx] = tt;
	}
	if (var < 0.9)
	{ /* include accel terms if EKF filter is converged */
		for (i = 0; i < 3; i++)
		{
			F[i + (i + 6) * nx] = SQR(tt) / 2.0;
		}
	}
	for (i = 0; i < nx; i++)
	{
		x[i] = rtk->x_spp[ix[i]];
		for (j = 0; j < nx; j++)
		{
			P[i + j * nx] = rtk->P_spp[ix[i] + ix[j] * NX_F];
		}
	}
	/* x=F*x, P=F*P*F+Q */
	matmul("NN", nx, 1, nx, 1.0, F, x, 0.0, xp);
	matmul("NN", nx, nx, nx, 1.0, F, P, 0.0, FP);
	matmul("NT", nx, nx, nx, 1.0, FP, F, 0.0, P);

	/*trace the mat*/
	trace(3, "xp=\n");
	tracemat(3, xp, 1, nx, 7, 3);
	trace(3, "Pp=\n");
	tracemat(3, P, nx, nx, 7, 3);

	for (i = 0; i < nx; i++)
	{
		rtk->x_spp[ix[i]] = xp[i];
		for (j = 0; j < nx; j++)
		{
			rtk->P_spp[ix[i] + ix[j] * NX_F] = P[i + j * nx];
		}
	}
	/* process noise added to only acceleration */
	Q[0] = Q[4] = SQR(1.5) * SQR(tt);
	Q[8] = SQR(2.0) * SQR(tt);
	ecef2pos(rtk->x_spp, pos);
	covecef(pos, Q, Qv);
	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		{
			rtk->P_spp[i + 6 + (j + 6) * NX_F] += Qv[i + j * 3];
		}
	free(ix);
	free(F);
	free(P);
	free(FP);
	free(x);
	free(xp);
}
static int udobs_rover(const obsd_t *obs, int n, rtk_t *rtk)
{
	int i,sat,freq;
	unsigned char track_index[MAXSAT];

	for (freq = 0; freq < rtk->opt.nf; freq++)
	{
		for (i = 0; i < MAXSAT; i++)
		{
			track_index[i] = 0;
		}
		for (i = 0; i < n && i < MAXOBS; i++)
		{
			sat = obs[i].sat;
			if (sat>0 && sat <= MAXSAT && obs[i].P[freq]!=0.0 && obs[i].SNR[freq]>0) track_index[sat - 1] = 1;
		}

		for (i = 0; i < MAXSAT; i++)
		{
			if (track_index[i] == 1)
			{
				rtk->ssat[i].lock_P[0][freq]++;
				trace(4, "udobs_rover: sat=%3d f=%d lock=%5d\n", i + 1, freq + 1, rtk->ssat[i].lock_P[0][freq]);
			}
			else rtk->ssat[i].lock_P[0][freq] = 0;
		}
	}
}
/* estimate receiver position ------------------------------------------------*/
static int estpos_ekf(const obsd_t *obs, int n, const double *rs, const double *dts,
	const double *vare, const int *svh, const nav_t *nav, double *azel, int *vsat,
	double *resp, rtk_t *rtk)
{
	double x[NX_F] = { 0 }, P[NX_F * NX_F], *v, *H, *var, *R, sig;
	int i, j, k, info, stat = 0, nv, ns;
	int row_num = n + 4 + n* rtk->opt.nf, col_num = NX_F;
	double dt = 0.0;
	char msg[128] = "";
	sol_t *sol = &rtk->sol;
	static int inittimes = 0;

	trace(3, "estpos_ekf: n=%d\n", n);

	dt = timediff(obs[0].time, sol->time);

	/* use LSQ spp result to init the filter state */
	if (norm(rtk->x_spp, 3) < 1000)
	{
		stat = pntpos(obs, n, nav, &rtk->opt, sol, azel, rtk->ssat, msg);

		if (stat == SOLQ_NONE || sol->ns < 6) return stat;

		if (inittimes++ > 30)
		{
			init_spp(rtk);
			for (i = 0; i < 3; i++)
				initx(rtk, sol->rr[i], VAR_POS, i);
			for (i = 3; i < 6; i++)
				initx(rtk, 0.1, VAR_VEL, i);
			for (i = 6; i < 9; i++)
				initx(rtk, 1E-2, VAR_ACC, i);

			for (i = 9; i < NX_F; i++)
			{
				initx(rtk, 0.1, SQR(30), i);
			}
		}
		return stat;
	}
	else
	{
		udobs_rover(obs, n, rtk);
		udpos_spp(rtk, dt);
		rtk->P_spp[9 + 9 * NX_F] += SQR(0.3);
		/* add the receiver clk process noise */
		for (i = 10; i < NX_F; i++){
			rtk->P_spp[i + i * NX_F] += SQR(0.1);
		}
	}

	if (rtk->opt.spp_mode == SPP_MODE_LX)
	{
		row_num = 2 * n * rtk->opt.nf + 4;
	}

	v = mat(row_num, 1);
	H = zeros(col_num, row_num);
	var = mat(row_num, 1);
	R = zeros(row_num, row_num);

	matcpy(x, rtk->x_spp, 1, NX_F);
	matcpy(P, rtk->P_spp, NX_F, NX_F);

	for (i = 0; i < 1; i++)
	{
		/* pseudorange residuals (m) */
		nv = rescode_mulfreq_ekf(1, obs, n, rs, dts, vare, svh, nav, x, &rtk->opt, v, H, var, azel, vsat, resp, &ns, rtk->ssat);

		if (ns < NX_F-6)
		{
			return -1;
		}

		nv += resdop_mulfreq_ekf(obs, n, rs, dts, nav, x, &rtk->opt, x + 3, azel, vsat, v + nv, var + nv, H + nv * NX_F);

		/* weighted by Std */
		for (j = 0; j < nv; j++)
		{
			R[j * nv + j] = var[j];
		}

		/*debug trace*/
		trace(3, "before_x=\n");
		tracemat(3, x, 1, col_num, 7, 3);
		trace(4, "before_P=\n");
		tracemat(4, P, col_num, col_num, 7, 3);
		trace(4, "H=\n");
		tracemat(4, H, col_num, nv, 7, 3);
		trace(4, "v=\n");
		tracemat(4, v, 1, nv, 7, 3);
		trace(4, "R=\n");
		tracemat(4, R, nv, nv, 7, 3);
		/* kalman estimation */
		if ((info = filter(x, P, H, v, R, NX_F, nv)))
		{
			return -2;
		}
		trace(3, "after_x=\n");
		tracemat(3, x, 1, col_num, 7, 3);
		trace(3, "after_P=\n");
		tracemat(3, P, col_num, col_num, 7, 3);
		stat = SOLQ_SINGLE;
	}

	matcpy(rtk->x_spp, x, 1, NX_F);
	matcpy(rtk->P_spp, P, NX_F, NX_F);

	sol->type = 0;
	sol->time = timeadd(obs[0].time, -x[9] / CLIGHT);
	sol->dtr[0] = x[9] / CLIGHT;  /* receiver clock bias (s) */
	sol->dtr[1] = x[10] / CLIGHT; /* GLO-GPS time offset (s) */
	sol->dtr[2] = x[11] / CLIGHT; /* GAL-GPS time offset (s) */
	sol->dtr[3] = x[12] / CLIGHT; /* BDS-GPS time offset (s) */
	sol->dtr[4] = x[13] / CLIGHT; /* IRN-GPS time offset (s) */
	for (j = 0; j < 6; j++)
	{
		sol->rr[j] = x[j];
	}
	for (j = 0; j < 3; j++)
	{
		sol->qr[j] = (float)P[j + j * NX_F];
	}
	sol->qr[3] = (float)P[1];        /* cov xy */
	sol->qr[4] = (float)P[2 + NX_F]; /* cov yz */
	sol->qr[5] = (float)P[2];        /* cov zx */
	for (j = 0; j < 3; j++)
	{
		sol->qv[j] = (float)P[(j + 3) + (j + 3) * NX_F];
	}
	sol->qv[3] = 0.0;        /* cov xy */
	sol->qv[4] = 0.0;        /* cov yz */
	sol->qv[5] = 0.0;        /* cov zx */
	sol->ns = (uint8_t)ns;
	sol->age = sol->ratio = 0.0;
	sol->stat = SOLQ_SINGLE;

	free(v);
	free(H);
	free(var);
	free(R);
	return stat;
}
/* RAIM FDE (failure detection and exclution) -------------------------------*/
static int raim_fde(const obsd_t *obs, int n, const double *rs,
                    const double *dts, const double *vare, const int *svh,
                    const nav_t *nav, const prcopt_t *opt, const ssat_t *ssat, 
                    sol_t *sol, double *azel, int *vsat, double *resp, char *msg)
{
    obsd_t *obs_e;
    sol_t sol_e={{0}};
    char tstr[32],name[16],msg_e[128];
    double *rs_e,*dts_e,*vare_e,*azel_e,*resp_e,rms_e,rms=100.0;
    int i,j,k,nvsat,stat=0,*svh_e,*vsat_e,sat=0;
    
    trace(3,"raim_fde: %s n=%2d\n",time_str(obs[0].time,0),n);
    
    if (!(obs_e=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    rs_e = mat(6,n); dts_e = mat(2,n); vare_e=mat(1,n); azel_e=zeros(2,n);
    svh_e=imat(1,n); vsat_e=imat(1,n); resp_e=mat(1,n); 
    
    for (i=0;i<n;i++) {
        
        /* satellite exclusion */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            obs_e[k]=obs[j];
            matcpy(rs_e +6*k,rs +6*j,6,1);
            matcpy(dts_e+2*k,dts+2*j,2,1);
            vare_e[k]=vare[j];
            svh_e[k++]=svh[j];
        }
        /* estimate receiver position without a satellite */
        if (!estpos(obs_e,n-1,rs_e,dts_e,vare_e,svh_e,nav,opt,ssat,&sol_e,azel_e,
                    vsat_e,resp_e,msg_e)) {
            trace(3,"raim_fde: exsat=%2d (%s)\n",obs[i].sat,msg);
            continue;
        }
        for (j=nvsat=0,rms_e=0.0;j<n-1;j++) {
            if (!vsat_e[j]) continue;
            rms_e+=SQR(resp_e[j]);
            nvsat++;
        }
        if (nvsat<5) {
            trace(3,"raim_fde: exsat=%2d lack of satellites nvsat=%2d\n",
                  obs[i].sat,nvsat);
            continue;
        }
        rms_e=sqrt(rms_e/nvsat);
        
        trace(3,"raim_fde: exsat=%2d rms=%8.3f\n",obs[i].sat,rms_e);
        
        if (rms_e>rms) continue;
        
        /* save result */
        for (j=k=0;j<n;j++) {
            if (j==i) continue;
            matcpy(azel+2*j,azel_e+2*k,2,1);
            vsat[j]=vsat_e[k];
            resp[j]=resp_e[k++];
        }
        stat=1;
        sol_e.eventime = sol->eventime;
        *sol=sol_e;
        sat=obs[i].sat;
        rms=rms_e;
        vsat[i]=0;
        strcpy(msg,msg_e);
    }
    if (stat) {
        time2str(obs[0].time,tstr,2); satno2id(sat,name);
        trace(2,"%s: %s excluded by raim\n",tstr+11,name);
    }
    free(obs_e);
    free(rs_e ); free(dts_e ); free(vare_e); free(azel_e);
    free(svh_e); free(vsat_e); free(resp_e);
    return stat;
}
/* range rate residuals ------------------------------------------------------*/
static int resdop(const obsd_t *obs, int n, const double *rs, const double *dts,
                  const nav_t *nav, const double *rr, const double *x,
                  const double *azel, const int *vsat, double err, double *v,
                  double *H)
{
    double freq,rate,pos[3],E[9],a[3],e[3],vs[3],cosel,sig;
    int i,j,nv=0;
    
    trace(3,"resdop  : n=%d\n",n);
    
    ecef2pos(rr,pos); xyz2enu(pos,E);
    
    for (i=0;i<n&&i<MAXOBS;i++) {
        
        freq=sat2freq(obs[i].sat,obs[i].code[0],nav);
        
        if (obs[i].D[0]==0.0||freq==0.0||!vsat[i]||norm(rs+3+i*6,3)<=0.0) {
            continue;
        }
        /* LOS (line-of-sight) vector in ECEF */
        cosel=cos(azel[1+i*2]);
        a[0]=sin(azel[i*2])*cosel;
        a[1]=cos(azel[i*2])*cosel;
        a[2]=sin(azel[1+i*2]);
        matmul("TN",3,1,3,1.0,E,a,0.0,e);
        
        /* satellite velocity relative to receiver in ECEF */
        for (j=0;j<3;j++) {
            vs[j]=rs[j+3+i*6]-x[j];
        }
        /* range rate with earth rotation correction */
        rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[0]-
                                      rs[3+i*6]*rr[1]-rs[  i*6]*x[1]);
        
        /* Std of range rate error (m/s) */
        sig=(err<=0.0)?1.0:err*CLIGHT/freq;
        
        /* range rate residual (m/s) */
        v[nv]=(-obs[i].D[0]*CLIGHT/freq-(rate+x[3]-CLIGHT*dts[1+i*2]))/sig;
        
        /* design matrix */
        for (j=0;j<4;j++) {
            H[j+nv*4]=((j<3)?-e[j]:1.0)/sig;
        }
        nv++;
    }
    return nv;
}
/* estimate receiver velocity ------------------------------------------------*/
static void estvel(const obsd_t *obs, int n, const double *rs, const double *dts,
                   const nav_t *nav, const prcopt_t *opt, sol_t *sol,
                   const double *azel, const int *vsat)
{
    double x[4]={0},dx[4],Q[16],*v,*H;
    double err=opt->err[4]; /* Doppler error (Hz) */
    int i,j,nv;
    
    v=mat(n,1); H=mat(4,n);
    
    for (i=0;i<MAXITR;i++) {
        
        /* range rate residuals (m/s) */
        if ((nv=resdop(obs,n,rs,dts,nav,sol->rr,x,azel,vsat,err,v,H))<4) {
            break;
        }
        /* least square estimation */
        if (lsq(H,v,4,nv,dx,Q)) break;
        
        for (j=0;j<4;j++) x[j]+=dx[j];
        
        if (norm(dx,4)<1E-6) {
            trace(3,"estvel : vx=%.3f vy=%.3f vz=%.3f, n=%d\n",x[0],x[1],x[2],n);
            matcpy(sol->rr+3,x,3,1);
            sol->qv[0]=(float)Q[0];  /* xx */
            sol->qv[1]=(float)Q[5];  /* yy */
            sol->qv[2]=(float)Q[10]; /* zz */
            sol->qv[3]=(float)Q[1];  /* xy */
            sol->qv[4]=(float)Q[6];  /* yz */
            sol->qv[5]=(float)Q[2];  /* zx */
            break;
        }
    }
    free(v); free(H);
}
/* single-point positioning ----------------------------------------------------
* compute receiver position, velocity, clock bias by single-point positioning
* with pseudorange and doppler observables
* args   : obsd_t *obs      I   observation data
*          int    n         I   number of observation data
*          nav_t  *nav      I   navigation data
*          prcopt_t *opt    I   processing options
*          sol_t  *sol      IO  solution
*          double *azel     IO  azimuth/elevation angle (rad) (NULL: no output)
*          ssat_t *ssat     IO  satellite status              (NULL: no output)
*          char   *msg      O   error message for error exit
* return : status(1:ok,0:error)
*-----------------------------------------------------------------------------*/
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel, ssat_t *ssat,
                  char *msg)
{
    prcopt_t opt_=*opt;
    double *rs,*dts,*var,*azel_,*resp;
    int i,stat,vsat[MAXOBS]={0},svh[MAXOBS];
    
    trace(3,"pntpos  : tobs=%s n=%d\n",time_str(obs[0].time,3),n);
    
    sol->stat=SOLQ_NONE;
    
    if (n<=0) {
        strcpy(msg,"no observation data");
        return 0;
    }
    sol->time=obs[0].time;
    msg[0]='\0';
    sol->eventime = obs[0].eventime;
    
    rs=mat(6,n); dts=mat(2,n); var=mat(1,n); azel_=zeros(2,n); resp=mat(1,n);
    
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].snr_rover[0]=0;
            ssat[i].snr_base[0]=0;
        }
        for (i=0;i<n;i++)
            ssat[obs[i].sat-1].snr_rover[0]=obs[i].SNR[0];
    }
    
    if (opt_.mode!=PMODE_SINGLE) { /* for precise positioning */
        opt_.ionoopt=IONOOPT_BRDC;
        opt_.tropopt=TROPOPT_SAAS;
    }
    /* satellite positons, velocities and clocks */
    satposs(sol->time,obs,n,nav,opt_.sateph,rs,dts,var,svh);
    
    /* estimate receiver position and time with pseudorange */
    stat=estpos(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    
    /* RAIM FDE */
    if (!stat&&n>=6&&opt->posopt[4]) {
        stat=raim_fde(obs,n,rs,dts,var,svh,nav,&opt_,ssat,sol,azel_,vsat,resp,msg);
    }
    /* estimate receiver velocity with Doppler */
    if (stat) {
        estvel(obs,n,rs,dts,nav,&opt_,sol,azel_,vsat);
    }
    if (azel) {
        for (i=0;i<n*2;i++) azel[i]=azel_[i];
    }
    if (ssat) {
        for (i=0;i<MAXSAT;i++) {
            ssat[i].vs=0;
            ssat[i].azel[0]=ssat[i].azel[1]=0.0;
            ssat[i].resp[0]=ssat[i].resc[0]=0.0;
        }
        for (i=0;i<n;i++) {
            ssat[obs[i].sat-1].azel[0]=azel_[  i*2];
            ssat[obs[i].sat-1].azel[1]=azel_[1+i*2];
            if (!vsat[i]) continue;
            ssat[obs[i].sat-1].vs=1;
            ssat[obs[i].sat-1].resp[0]=resp[i];
        }
    }
    free(rs); free(dts); free(var); free(azel_); free(resp);
    return stat;
}

extern int pntpos_ekf(const obsd_t *obs, const int n, const nav_t *nav, rtk_t *rtk)
{
	double *rs, *dts, *var, *azel_, *resp,pos[3];
	int i, stat, vsat[MAXOBS] = { 0 }, svh[MAXOBS];

	trace(3, "pntpos_ekf: tobs=%s nobs=%d\n", time_str(obs[0].time, 3), n);
	/* init spp state*/
	if (rtk->x_spp == NULL || rtk->P_spp == NULL){
		init_spp(rtk);
	}

	if (n <= 0)
	{
		return 0;
	}

	rs = mat(6, n);
	dts = mat(2, n);
	var = mat(1, n);
	azel_ = zeros(2, n);
	resp = mat(1, n);

	/* satellite positons, velocities and clocks */
	satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

	/* estimate receiver position with pseudorange */
	stat = estpos_ekf(obs, n, rs, dts, var, svh, nav, azel_, vsat, resp, rtk);
	ecef2pos(rtk->sol.rr, pos);
	trace(3, "estpos_ekf: stat=%d time=%s position=%.8lf %.8lf %.2lf\n", stat, time_str(obs[0].time, 3), pos[0] * R2D, pos[1] * R2D, pos[2]);
	printf("estpos_ekf: stat=%d time=%s nobs=%d position=%.8lf %.8lf %.2lf\n", stat, time_str(obs[0].time, 3), n, pos[0] * R2D, pos[1] * R2D, pos[2]);

	free(rs);
	free(dts);
	free(var);
	free(azel_);
	free(resp);
	return stat;
}

