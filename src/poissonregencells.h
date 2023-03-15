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
	PoissonRegenCells(int randSeed, std::fstream &psth_file_buf);
	~PoissonRegenCells();

  void calcGRPoissActivity(size_t ts);
	//void fill_rasters(uint32_t ts);
	void fill_psths(size_t ts);
	//void save_rasters(std::string out_file);
	void save_psths(std::string out_file);
	const uint8_t *getGRAPs();
	const float **getGRFRs();
	uint32_t **getApBufGR();
	uint64_t **getApHistGR();

private:
	void init_fr_from_file(std::fstream &input_file_buf);
	void init_templates_from_psth_file(std::fstream &input_psth_file_buf);
  void initGRCUDA();
  void initCUDAStreams();
  void initCURAND();

	std::normal_distribution<float> *normDist;
	std::mt19937 *noiseRandGen;
	CRandomSFMT0 *randSeedGen;
	CRandomSFMT0 **randGens;
	curandStateMRG32k3a **mrg32k3aRNGs; // device randGens

	uint32_t gpuIndStart = 0;
  uint64_t numGPUs = 1;
  cudaStream_t **streams;

	float **grActRandNums;

	unsigned int nThreads;

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
	float **threshs_d;
	uint8_t **aps_d;

	uint32_t **apBufs;
	uint64_t **apHists;
};

#endif /* POISSONREGENCELLS_H_ */

