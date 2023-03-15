/*
 * poissonregencells.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: consciousness
 *  Modified on: Jul 22, 2015
 *  	Author: evandelord
 */
#include <cstring>
#include <omp.h>

#include "logger.h"
#include "array_util.h"
#include "dynamic2darray.h"
#include "poissonregencells.h"

PoissonRegenCells::PoissonRegenCells(int randSeed, std::fstream &psth_file_buf)
{
	randSeedGen = new CRandomSFMT0(randSeed);

	nThreads=1;
	randGens=new CRandomSFMT0*[nThreads];

	for(unsigned int i=0; i<nThreads; i++)
	{
		randGens[i] = new CRandomSFMT0(randSeedGen->IRandom(0, INT_MAX));
	}

	// threshes need to be scaled
	threshBase = 0;
	threshMax  = 1;
	threshIncTau = 0.1; // hard code in the decay time constant
	threshInc = 1 - exp(-1.0 / threshIncTau); // for now, hard code in time-bin size
	sPerTS = msPerTimeStep / 1000;

	expansion_factor = 128; // by what factor are we expanding the granule cell number?
  num_gr_old = num_gr / expansion_factor;

	threshs_h = (float *)calloc(num_gr, sizeof(float));
	aps_h     = (uint8_t *)calloc(num_gr, sizeof(uint8_t));

  for (size_t i = 0; i < num_gr; i++) threshs_h[i]  = threshMax;

	//psths = allocate2DArray<uint8_t>(2000, psth_sample_size);
	//rasters = allocate2DArray<uint8_t>(2000, num_gr); // hardcoded, 2000 ts by 1000 trials by num_gr grs
	gr_templates_h = allocate2DArray<float>(num_gr / expansion_factor, 2000); // for now, only have firing rates from trials of 2000
	init_templates_from_psth_file(psth_file_buf);
  initGRCUDA();
}

PoissonRegenCells::~PoissonRegenCells()
{
	delete randSeedGen;
	for(uint32_t i=0; i<nThreads; i++)
	{
		delete randGens[i];
	}
	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaFree(mrg32k3aRNGs[i]);
		cudaFree(grActRandNums[i]);
		cudaDeviceSynchronize();
	}
	delete[] mrg32k3aRNGs;
	delete[] grActRandNums;

	delete[] randGens;

	//delete2DArray(psths);
	//delete2DArray(rasters);
  free(threshs_h);
  free(aps_h);
	delete2DArray(gr_templates_h);
	delete2DArray(gr_templates_t_h);

	for (uint32_t i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
    cudaFree(gr_templates_t_d[i]);
    cudaFree(threshs_d[i]);
    cudaFree(aps_d[i]);
		cudaDeviceSynchronize();
	}
  free(gr_templates_t_d);
  free(gr_templates_t_pitch);

  free(threshs_d);
  free(aps_d);
}

// Soon to be deprecated: was for loading in pre-computed smoothed-fr
//void PoissonRegenCells::init_fr_from_file(std::fstream &input_file_buf)
//{
//	input_file_buf.read((char *)gr_fr[0], num_gr * 2000 * sizeof(float));
//}

void PoissonRegenCells::init_templates_from_psth_file(std::fstream &input_psth_file_buf)
{
	/*
	 * Plan for the algorithm
	 *
	 * 1. read in psth data into a single array
	 * 2. determine categories:
	 *    a. which indices correspond to high firers during CS?
	 *    b. which indices correspond to cells that do not fire at all?
	 *    c. which indices correspond to cells that do fire, but below some cut-off avg fr?
	 * 3. reduce the number of templates given above categories:
	 *    a. full representation of gr templates from high CS firers
	 *    b. truncation of templates to 1 single "zero" template for low and non-firers,
	 *       or two separate: low fire and zero fire template
	 * 4. randomly assign each template to 100x its amount in the expanded sim
	 * 5. profit
	 */

	enum {ZERO_FIRERS, LOW_FIRERS, HIGH_FIRERS};
	uint32_t num_fire_categories[3] = {0};
	// expect input data comes from 1000 trials, 2000 ts a piece.
	const float THRESH_FR = 20.0; // in Hz
	const uint32_t THRESH_COUNT = (THRESH_FR / 1000.0) * 2000 * 1000;
	uint8_t **input_psths = allocate2DArray<uint8_t>(2000, num_gr_old);
	uint8_t *firing_categories = (uint8_t *)calloc(num_gr_old, sizeof(uint8_t));
	uint32_t *cell_spike_sums = (uint32_t *)calloc(num_gr_old, sizeof(uint32_t));
	LOG_INFO("loading cells from file.");
	input_psth_file_buf.read((char *)input_psths[0], 2000 * num_gr_old * sizeof(uint8_t));
	LOG_INFO("finished load.");
	LOG_INFO("transposing...");
	uint8_t **input_psths_t = transpose2DArray(input_psths,  2000, num_gr_old);
	LOG_INFO("finished transposing.");

	// determine firing rate categories
	LOG_INFO("determining firing categories...");
	for (size_t i = 0; i < num_gr_old; i++)
	{
		for (size_t j = 0; j < 2000; j++)
		{
			cell_spike_sums[i] += input_psths_t[i][j]; // add in the current time-steps accumulated spikes
		}
		if (cell_spike_sums[i] > THRESH_COUNT)
		{
			firing_categories[i] = HIGH_FIRERS;
			num_fire_categories[HIGH_FIRERS]++;
		}
		else if (cell_spike_sums[i] < THRESH_COUNT && cell_spike_sums[i] > 0)
		{
			firing_categories[i] = LOW_FIRERS;
			num_fire_categories[LOW_FIRERS]++;
		}
		else
		{
			firing_categories[i] = ZERO_FIRERS;
			num_fire_categories[ZERO_FIRERS]++;
		}
	}
	LOG_INFO("finished determining firing categories.");
	LOG_INFO("num high: %u", num_fire_categories[HIGH_FIRERS]);
	LOG_INFO("num low: %u", num_fire_categories[LOW_FIRERS]);
	LOG_INFO("num zero: %u", num_fire_categories[ZERO_FIRERS]);

	// create template pdfs
	LOG_INFO("creating templates...");
	for (size_t i = 0; i < num_gr_old; i++)
	{
		// NOTE: no need to include zero_firers as calloc sets their values to zero
		if (firing_categories[i] == LOW_FIRERS)
		{
			float mean_prob = (float)cell_spike_sums[i] / (2000 * 1000); // DEBUG
			for (size_t j = 0; j < 2000; j++)
			{
				gr_templates_h[i][j] = mean_prob; // just set every bin to same avg fr
			}
		}
		else // high firers
		{
			for (size_t j = 0; j < 2000; j++)
			{
				gr_templates_h[i][j] = (float)input_psths_t[i][j] / cell_spike_sums[i]; // the easiest transformation into a pdf
			}
		}
	}
	LOG_INFO("finished creating templates...");

	LOG_INFO("transposing templates...");
	gr_templates_t_h = transpose2DArray(gr_templates_h, num_gr_old, 2000);
	LOG_INFO("finished transposing templates.");

	delete2DArray(input_psths);
	delete2DArray(input_psths_t);
	free(firing_categories);
	free(cell_spike_sums);
}

void PoissonRegenCells::initGRCUDA()
{
	numGRPerGPU = num_gr / numGPUs;

  calcGRActNumGRPerB = 1024; // 1024 ie max num threads per block
  calcGRActNumBlocks = numGRPerGPU / calcGRActNumBlocks; // gives 2 ^ 17

  initCUDAStreams();
  initCURAND();

  gr_templates_t_d     = (float **)calloc(numGPUs, sizeof(float *));
  gr_templates_t_pitch = (size_t *)calloc(numGPUs, sizeof(size_t));

	threshs_d = (float **)calloc(numGPUs, sizeof(float *));
	aps_d     = (uint8_t **)calloc(numGPUs, sizeof(uint8_t *));

	LOG_INFO("Allocating device memory...");
	for (uint32_t i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaMallocPitch(&gr_templates_t_d[i], &gr_templates_t_pitch[i],
			num_gr_old * sizeof(float), 2000);
	  LOG_INFO("gr templates malloc pitch error: %s", cudaGetErrorString(cudaGetLastError()));
    cudaMalloc(&threshs_d[i], num_gr * sizeof(float));
	  LOG_INFO("threshs malloc error: %s", cudaGetErrorString(cudaGetLastError()));
    cudaMalloc(&aps_d[i], num_gr * sizeof(uint8_t));
	  LOG_INFO("aps malloc error: %s", cudaGetErrorString(cudaGetLastError()));
		cudaDeviceSynchronize();
	}
	LOG_INFO("alloc final error: %s", cudaGetErrorString(cudaGetLastError()));
	LOG_INFO("Finished allocating device memory.");


	LOG_INFO("Copying host to device memory...");
  for (uint32_t i = 0; i < numGPUs; i++)
  {
    cudaMemcpy2D(gr_templates_t_d[i], gr_templates_t_pitch[i], gr_templates_t_h,
      num_gr_old * sizeof(float), num_gr_old * sizeof(float), 2000, cudaMemcpyHostToDevice);
	  LOG_INFO("gr templates memcpy2D error: %s", cudaGetErrorString(cudaGetLastError()));
		cudaMemcpy(threshs_d[i], threshs_h, num_gr * sizeof(float), cudaMemcpyHostToDevice);
	  LOG_INFO("threshs memcpy error: %s", cudaGetErrorString(cudaGetLastError()));
		cudaMemcpy(aps_d[i], aps_h, num_gr * sizeof(uint8_t), cudaMemcpyHostToDevice);
	  LOG_INFO("aps memcpy error: %s", cudaGetErrorString(cudaGetLastError()));
  }
	LOG_INFO("memcpy final error: %s", cudaGetErrorString(cudaGetLastError()));
	LOG_INFO("Finished copying host to device memory.");
}

void PoissonRegenCells::initCUDAStreams()
{
	cudaError_t error;

	int maxNumGPUs;
	// TODO: use assert, try, and catch for these types of errors
	error = cudaGetDeviceCount(&maxNumGPUs);

	LOG_INFO("CUDA max num devices: %d", maxNumGPUs);
	LOG_INFO("%s", cudaGetErrorString(error));
	LOG_INFO("CUDA num devices: %d, starting at GPU %d", numGPUs, gpuIndStart);

	streams = new cudaStream_t*[numGPUs];

	for (int i = 0; i < numGPUs; i++)
	{
		error = cudaSetDevice(i + gpuIndStart);
		LOG_INFO("Selecting device %d", i);
		LOG_INFO("%s", cudaGetErrorString(error));
		streams[i] = new cudaStream_t[8];
		LOG_INFO("Resetting device %d", i);
		LOG_INFO("%s", cudaGetErrorString(error));
		cudaDeviceSynchronize();

		for (int j = 0; j < 8; j++)
		{
			error = cudaStreamCreate(&streams[i][j]);
			LOG_INFO("Initializing stream %d for device %d",j, i);
			LOG_INFO("%s", cudaGetErrorString(error));
		}
		cudaDeviceSynchronize();
		error = cudaGetLastError();
		LOG_INFO("Cuda device %d", i);
		LOG_INFO("%s", cudaGetErrorString(error));
	}
}

void PoissonRegenCells::initCURAND()
{
	// set up rng
	LOG_INFO("Initializing curand state...");

	CRandomSFMT cudaRNGSeedGen(time(0));

	int32_t curandInitSeed = cudaRNGSeedGen.IRandom(0, INT_MAX);

	mrg32k3aRNGs = new curandStateMRG32k3a*[numGPUs];
	grActRandNums = new float*[numGPUs];

  // we're just playing with some numbers rn
  numGRPerRandBatch = 32768; // 2 ^ 17
  updateGRRandNumGRPerB = 512;// gives 1024, max num thread per block
  updateGRRandNumBlocks = numGRPerRandBatch / updateGRRandNumGRPerB;
	dim3 updateGRRandGridDim(updateGRRandNumBlocks);
	dim3 updateGRRandBlockDim(updateGRRandNumGRPerB);

	for (uint8_t i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		// FIXME: instead of set seed, bring in a rand gen seed
		// FIXME: either need to access the threads for which
		//        mrg32k3aRNGs have been allocated by considering
		//        an offset in the plasticity kernel, or create
		//        enough rngs to cover that entire space
		cudaMalloc((void **)&mrg32k3aRNGs[i],
				updateGRRandNumGRPerB * updateGRRandNumBlocks * sizeof(curandStateMRG32k3a));
	  LOG_INFO("rng malloc error: %s", cudaGetErrorString(cudaGetLastError()));
		callCurandSetupKernel<curandStateMRG32k3a, dim3, dim3>
		(streams[i][1], mrg32k3aRNGs[i], curandInitSeed, updateGRRandGridDim, updateGRRandBlockDim);
	  LOG_INFO("curand setup error: %s", cudaGetErrorString(cudaGetLastError()));
		cudaMalloc((void **)&grActRandNums[i], numGRPerGPU * sizeof(float));
	  LOG_INFO("rand num malloc error: %s", cudaGetErrorString(cudaGetLastError()));
		cudaMemset(grActRandNums[i], 0, numGRPerGPU * sizeof(float));
	  LOG_INFO("rand num memset error: %s", cudaGetErrorString(cudaGetLastError()));
		cudaDeviceSynchronize();
	}
	LOG_INFO("Finished initializing curand state.");
	LOG_INFO("Last error: %s", cudaGetErrorString(cudaGetLastError()));
}

void PoissonRegenCells::calcGRPoissActivity(size_t ts)
{
  int gpuIndStart = 0;
  cudaError_t error;

	for(int i = 0; i < numGPUs; i++)
	{
    // generate the random numbers in batches
	  for (int j = 0; j < num_gr; j += numGRPerRandBatch)
	  {
	    error = cudaSetDevice(i + gpuIndStart);
	    callCurandGenerateUniformKernel<curandStateMRG32k3a>(streams[i][3],
	    	mrg32k3aRNGs[i], updateGRRandNumBlocks, updateGRRandNumGRPerB,
	    	grActRandNums[i], j);
    }
    callGRActKernel(streams[i][1], calcGRActNumBlocks, calcGRActNumGRPerB,
      threshs_d[i], aps_d[i], grActRandNums[i], gr_templates_t_d[i], gr_templates_t_pitch[i], num_gr_old,
      threshBase, threshMax, threshInc);
  }
}

//void PoissonRegenCells::fill_rasters(uint32_t ts)
//{
//	memcpy(rasters[ts], aps, num_gr * sizeof(uint8_t));
//}

//void PoissonRegenCells::fill_psths(size_t ts)
//{
//	for (size_t i = 0; i < num_gr; i++)
//	{
//		psths[ts][i] += aps[i];
//	}
//}

//void PoissonRegenCells::save_rasters(std::string out_file)
//{
//	write2DArray<uint8_t>(out_file, rasters, 2000, num_gr);
//}

void PoissonRegenCells::save_psths(std::string out_file)
{
	write2DArray<uint8_t>(out_file, psths, 2000, num_gr);
}

const uint8_t *PoissonRegenCells::getGRAPs() { return (const uint8_t *)aps_h; }

const float **PoissonRegenCells::getGRFRs() { return (const float **)gr_fr; }

uint32_t **PoissonRegenCells::getApBufGR() { return apBufs; }

uint64_t **PoissonRegenCells::getApHistGR() { return apHists; }

