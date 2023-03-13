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
	PoissonRegenCells(int randSeed, std::fstream &psth_file_buf);
	~PoissonRegenCells();

	void calcGRPoissActivity(uint32_t ts);
	void calcGRPoissActivitySample(uint32_t ts);
	//void fill_rasters(uint32_t ts);
	void fill_psths(size_t ts);
	void fill_psths_sample(size_t ts);
	void save_sample_template_indices(std::string out_file);
	void save_template_indices(std::string out_file);
	//void save_rasters(std::string out_file);
	void save_psths(std::string out_file);
	void save_psths_sample(std::string out_file);
	const uint8_t *getGRAPs();
	const float **getGRFRs();
	uint32_t **getApBufGR();
	uint64_t **getApHistGR();

private:
	void init_fr_from_file(std::fstream &input_file_buf);
	void init_templates_from_psth_file(std::fstream &input_psth_file_buf);
	std::normal_distribution<float> *normDist;
	std::mt19937 *noiseRandGen;
	CRandomSFMT0 *randSeedGen;
	CRandomSFMT0 **randGens;

	unsigned int nThreads;

	float threshBase;
	float threshMax;
	float threshIncTau;
	float threshInc;
	float sPerTS;
	uint32_t expansion_factor;

	uint32_t psth_sample_size;
	// template indices indicate from which *input* cell template do we generate spikes
	// from. Will have dimensions num_gr / expansion_factor
	size_t *sample_template_indices;
	size_t *template_indices;

	// watch out: these are indices of expanded number of gr cells, so num_gr
	// (which is expanded by expansion_factor from inputs)
	size_t *sample_indices;
	
	uint8_t **psths;
	//uint8_t **rasters;
	float **gr_fr;
	float **gr_templates;
	// i made this var because is mildly faster :-) come back in 46 days pls and thank
	float **gr_templates_t;
	float *threshs;
	uint8_t *aps;
	uint32_t **apBufs;
	uint64_t **apHists;
	int spikeTimer = 0;
};

#endif /* POISSONREGENCELLS_H_ */

