/*
 * poissonregencells.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: consciousness
 *  Modified on: Jul 22, 2015
 *  	Author: evandelord
 */
#include <cstring>

#include "dynamic2darray.h"
#include "poissonregencells.h"

PoissonRegenCells::PoissonRegenCells(int randSeed, std::fstream &fr_file_buf)
{
	randSeedGen = new CRandomSFMT0(randSeed);

	nThreads=1;
	randGens=new CRandomSFMT0*[nThreads];

	for(unsigned int i=0; i<nThreads; i++)
	{
		randGens[i] = new CRandomSFMT0(randSeedGen->IRandom(0, INT_MAX));
	}

	sPerTS = msPerTimeStep / 1000;

	psths = allocate2DArray<uint8_t>(1300, num_gr);
	rasters = allocate2DArray<uint8_t>(1300, num_gr); // hardcoded, 1300 ts by 1000 trials by num_gr grs
	gr_fr = allocate2DArray<float>(num_gr, 1300); // for now, only have firing rates from trials of 1300
	aps  = (uint8_t *)calloc(num_gr, sizeof(uint8_t));
	init_fr_from_file(fr_file_buf);
}

PoissonRegenCells::~PoissonRegenCells()
{
	delete randSeedGen;
	for(uint32_t i=0; i<nThreads; i++)
	{
		delete randGens[i];
	}

	delete[] randGens;
	delete2DArray(psths);
	delete2DArray(rasters);
	delete2DArray(gr_fr);
	free(aps);
}

void PoissonRegenCells::init_fr_from_file(std::fstream &input_file_buf)
{
	input_file_buf.read((char *)gr_fr[0], num_gr * 1300 * sizeof(float));
}

void PoissonRegenCells::calcGRPoissActivity(uint32_t ts)
{
	for (uint32_t i = 0; i < num_gr; i++)
	{
		aps[i] = (randGens[0]->Random() < (gr_fr[i][ts]) * sPerTS);
	}
}

void PoissonRegenCells::fill_rasters(uint32_t ts)
{
	memcpy(rasters[ts], aps, num_gr * sizeof(uint8_t));
}

void PoissonRegenCells::fill_psths(uint32_t ts)
{
	for (uint32_t i = 0; i < num_gr; i++)
	{
		psths[ts][i] += aps[i];
	}
}

void PoissonRegenCells::save_rasters(std::string out_file)
{
	write2DArray<uint8_t>(out_file, rasters, 1300, num_gr);
}

void PoissonRegenCells::save_psths(std::string out_file)
{
	write2DArray<uint8_t>(out_file, psths, 1300, num_gr);
}

/*
 * Note: should be called something like "calcCollateralMFs." also, 
 * depending on the sets of mfs, why calculate this every time step. waste of time.
 */
//bool* PoissonRegenCells::calcTrueGRs(const float *frequencies)
//{
//	for (uint32_t i = 0; i < num_gr; i++)
//	{
//		if (frequencies[i] == -1) // indicates that this mf was not assigned a frequency,
//								  // so is a dcn collateral
//		{
//			isTrueGR[i] = false;
//		}
//	}
//	return isTrueGR;
//}

const uint8_t *PoissonRegenCells::getGRAPs() { return (const uint8_t *)aps; }

const float **PoissonRegenCells::getGRFRs() { return (const float **)gr_fr; }

uint32_t **PoissonRegenCells::getApBufGR() { return apBufs; }

uint64_t **PoissonRegenCells::getApHistGR() { return apHists; }

