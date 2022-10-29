/*------------------------------------------------------------------------------
*  GNSS_Algorithm.cpp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2022 by BruceZhcw, All rights reserved.
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include <iostream>
#include "rtklib.h"

#define PROGNAME    "GNSS_Algorithm"	/* program name */
#define MAXFILE     16                  /* max number of input files */

/* rnx2rtkp main -------------------------------------------------------------*/
int main(int argc, char **argv)
{
	prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = { "" };
	gtime_t ts = { 0 }, te = { 0 };
	double tint = 0.0, es[] = { 2022, 7, 21, 06, 44, 44 }, ee[] = { 2022, 7, 21, 06, 49, 01 }, pos[3] = {40.0680091,116.3355171,46};
	int i, j, n, ret;
	char *infile[MAXFILE] = {"F:\\GNSS_RTK\\GNSS_Algorithm\\GNSS_Algorithm\\data\\rover.obs",
							 "F:\\GNSS_RTK\\GNSS_Algorithm\\GNSS_Algorithm\\data\\rover.nav"};
	char *outfile =			{"F:\\GNSS_RTK\\GNSS_Algorithm\\GNSS_Algorithm\\data\\output.txt"};

	solopt.posf = SOLF_LLH;
	solopt.timef = 1;
	sprintf(solopt.prog, "%s ver.%s %s", PROGNAME, VER_RTKLIB, PATCH_LEVEL);
	sprintf(filopt.trace, "F:\\GNSS_RTK\\GNSS_Algorithm\\GNSS_Algorithm\\data\\%s.trace", PROGNAME);

	//ts = epoch2time(es);
	//te = epoch2time(ee);

	prcopt.snrmask.ena[0] = prcopt.snrmask.ena[1] = 1;
	for (i = 0; i < NFREQ; i++) for (j = 0; j < 9; j++)
		prcopt.snrmask.mask[i][j] = 20;

	prcopt.mode = PMODE_SINGLE;
	prcopt.spp_mode = SPP_MODE_LX;

	for (j = 0; j<2; j++) pos[j] *= D2R;
	pos2ecef(pos, prcopt.rb);
	matcpy(prcopt.ru, prcopt.rb, 3, 1);

	ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, 2, outfile, "", "");

	std::cout << "postpos over! return with: " << ret << std::endl;
	getchar();

	return ret;
}

