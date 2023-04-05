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
#include "kernels.h"

class PoissonRegenCells
{
public:
	PoissonRegenCells();
	PoissonRegenCells(std::fstream &psth_file_buf,
					  cudaStream_t **streams);
	~PoissonRegenCells();

	void calcGRPoissActivity(size_t ts, cudaStream_t **streams, uint8_t streamN);
	//void fill_rasters(uint32_t ts);
	void fill_psths(size_t ts);
	//void save_rasters(std::string out_file);
	void save_psths(std::string out_file);
	const uint8_t *getGRAPs();
	const float **getGRFRs();
	uint32_t **get_ap_buf_gr_gpu();
	uint64_t **get_ap_hist_gr_gpu();

private:
	void init_fr_from_file(std::fstream &input_file_buf);
	void init_templates_from_psth_file(std::fstream &input_psth_file_buf);
	void initGRCUDA();
	void initCURAND(cudaStream_t **streams);

	std::normal_distribution<float> *normDist;
	std::mt19937 *noiseRandGen;
	CRandomSFMT0 *randSeedGen;
	CRandomSFMT0 **randGens;
	curandStateMRG32k3a **mrg32k3aRNGs; // device randGens

	uint32_t gpuIndStart = 0;
  uint64_t numGPUs = 2;

	float **grActRandNums;

	unsigned int nThreads;
	uint32_t num_template_ts;
	uint32_t num_trials;

  uint64_t numGROldPerGPU;
  uint64_t numGRPerGPU;
  uint32_t calcGRActNumBlocks;
  uint32_t calcGRActNumGRPerB;

  uint64_t numGRPerRandBatch;
  uint64_t updateGRRandNumBlocks;
  uint64_t updateGRRandNumGRPerB;

	float threshBase;
	float threshMax;
	float threshIncTau;
	float threshInc;
	float sPerTS;
	uint32_t expansion_factor;
  size_t num_gr_old;

	uint8_t **psths;
	//uint8_t **rasters;
	float **gr_fr;
	float **gr_templates_h;
	// i made this var because is mildly faster :-) come back in 46 days pls and thank
	float **gr_templates_t_h;
  float **gr_templates_t_d;
  size_t *gr_templates_t_pitch;

	float *threshs_h;
	uint8_t *aps_h;
	uint32_t *aps_buf_h;
	uint64_t *aps_hist_h;

	float **threshs_d;
	uint8_t **aps_d;
	uint32_t **aps_buf_d;

	uint64_t **aps_hist_d;
};

#endif /* POISSONREGENCELLS_H_ */

