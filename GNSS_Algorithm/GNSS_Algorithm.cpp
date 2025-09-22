/*------------------------------------------------------------------------------
*  GNSS_Algorithm.cpp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2022 by BruceZhcw, All rights reserved.
*-----------------------------------------------------------------------------*/
#include <iostream>
#include <cstring>
#include "rtklib.h"

#define PROGNAME    "GNSS_Algorithm"	/* program name */
#define MAXFILE     16                  /* max number of input files */

int main(int argc, char **argv)
{
	if (argc < 2) {
		std::cerr << "Usage:   " << argv[0] << " <rinex directory>" << std::endl;
		std::cerr << "example: " << "./gnss_algorithm  /mnt/e/gnss/" << std::endl;
		return -1;
	}

	char directory[256];
	strcpy(directory, argv[1]);
	if (directory[strlen(directory)-1] != '/') {
		strcat(directory, "/");
	}

	prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = { "" };
	gtime_t ts = { 0 }, te = { 0 };
	double tint = 0.0, pos[3] = {40.0680091,116.3355171,46};
	int i, j, ret;

	char *infile[MAXFILE];
	char filepath1[512], filepath2[512], outfilepath[512], tracepath[512];

	sprintf(filepath1, "%srover.obs", directory);
	sprintf(filepath2, "%srover.nav", directory);
	sprintf(outfilepath, "%srover.pos", directory);
	sprintf(tracepath, "%s%s.trace", directory, PROGNAME);

	infile[0] = filepath1;
	infile[1] = filepath2;
	char *outfile = outfilepath;

	solopt.posf = SOLF_LLH;
	solopt.timef = 1;
	sprintf(solopt.prog, "%s ver.%s %s", PROGNAME, VER_RTKLIB, PATCH_LEVEL);
	strcpy(filopt.trace, tracepath);

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

	return ret;
}

