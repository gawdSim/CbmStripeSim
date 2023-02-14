/*
 * poissonregencells.h
 *
 *  Created on: Dec 13, 2012
 *      Author: consciousness
 */

#ifndef POISSONREGENCELLS_H_
#define POISSONREGENCELLS_H_

#include <iostream>
#include <fstream>
#include <algorithm> // for random_shuffle
#include <cstdlib> // for srand and rand, sorry Wen
#include <random>
#include <math.h>
#include <limits.h>

#include <cstdint>
#include "sfmt.h"
#include "mzone.h"
#include "connectivityparams.h"
#include "activityparams.h"

class PoissonRegenCells
{
public:
	PoissonRegenCells();
	PoissonRegenCells(int randSeed, std::fstream &fr_file_buf);
	~PoissonRegenCells();

	void init_fr_from_file(std::fstream &input_file_buf);
	void calcGRPoissActivity(uint32_t ts);
	void fill_rasters(uint32_t ts);
	void fill_psths(uint32_t ts);
	void save_rasters(std::string out_file);
	void save_psths(std::string out_file);
	const uint8_t *getGRAPs();
	const float **getGRFRs();
	uint32_t **getApBufGR();
	uint64_t **getApHistGR();

private:

	std::normal_distribution<float> *normDist;
	std::mt19937 *noiseRandGen;
	CRandomSFMT0 *randSeedGen;
	CRandomSFMT0 **randGens;

	unsigned int nThreads;

	float sPerTS;

	uint8_t **psths;
	uint8_t **rasters;
	float **gr_fr;
	uint8_t *aps;
	uint32_t **apBufs;
	uint64_t **apHists;
	int spikeTimer = 0;
};

#endif /* POISSONREGENCELLS_H_ */

