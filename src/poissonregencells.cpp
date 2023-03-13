/*
 * poissonregencells.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: consciousness
 *  Modified on: Jul 22, 2015
 *  	Author: evandelord
 */
#include <cstring>

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
	threshDecTau = 0.1; // hard code in the decay time constant
	threshDec = 1 - exp(-1.0 / threshDecTau); // for now, hard code in time-bin size
	sPerTS = msPerTimeStep / 1000;

	psths = allocate2DArray<uint8_t>(2000, num_gr);
	rasters = allocate2DArray<uint8_t>(2000, num_gr); // hardcoded, 2000 ts by 1000 trials by num_gr grs
	gr_pdfs = allocate2DArray<float>(2000, num_gr); // for now, only have firing rates from trials of 2000
	threshs = (float *)calloc(num_gr, sizeof(float));
	aps  = (uint8_t *)calloc(num_gr, sizeof(uint8_t));
	init_templates_from_psth_file(psth_file_buf);
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
	delete2DArray(gr_pdfs);
	free(aps);
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

	//TODO: debug :-) 
	enum {ZERO_FIRERS, LOW_FIRERS, HIGH_FIRERS};
	uint32_t num_fire_categories[3];
	// expect input data comes from 1000 trials, 2000 ts a piece.
	const float THRESH_FR = 20; // in Hz
	const uint32_t THRESH_COUNT = (uint32_t)(THRESH_FR / 1000.0) * 2000 * 1000;
	uint32_t expansion_factor = 128; // by what factor are we expanding the granule cell number?
	uint8_t **input_psths = allocate2DArray<uint8_t>(2000, num_gr / expansion_factor);
	uint8_t **input_pdf_templates = allocate2DArray<uint8_t>(2000, num_gr / expansion_factor);
	uint8_t *firing_categories = (uint8_t *)calloc(num_gr / expansion_factor, sizeof(uint8_t));
	uint32_t *cell_spike_sums = (uint32_t *)calloc(num_gr / expansion_factor, sizeof(uint32_t));
	input_psth_file_buf.read((char *)input_psths[0], 2000 * num_gr / expansion_factor * sizeof(uint8_t));

	// determine firing rate categories
	for (uint32_t i = 0; i < num_gr / expansion_factor; i++)
	{
		for (uint32_t j = 0; j < 2000; j++)
		{
			cell_spike_sums[i] += input_psths[j][i]; // add in the current time-steps accumulated spikes
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

	// create template pdfs
	for (uint32_t i = 0; i < num_gr / expansion_factor; i++)
	{
		// NOTE: no need to include zero_firers as calloc sets their values to zero
		if (firing_categories[i] == LOW_FIRERS)
		{
			float mean_prob = (float)cell_spike_sums[i] / (2000 * 1000); // DEBUG
			for (uint32_t j = 0; j < 2000; j++)
			{
				input_pdf_templates[j][i] = mean_prob; // just set every bin to same avg fr
			}
		}
		else // high firers
		{
			for (uint32_t j = 0; j < 2000; j++)
			{
				input_pdf_templates[j][i] = (float)input_psths[j][i] / cell_spike_sums[i]; // the easiest transformation into a pdf
			}
		}
	}

	// assign template pdfs to expanded gr set deterministically
	uint32_t expansion_counter = 0;
	for (uint32_t i = 0; i < num_gr / expansion_factor; i++)
	{
		if (firing_categories[i] != ZERO_FIRERS) 
		{
			for (uint32_t j = 0; j < expansion_factor; j++)
			{
				for (uint32_t k = 0; k < 2000; k++)
				{
					gr_pdfs[k][expansion_counter + j] = input_pdf_templates[k][i];
				}
			}
		}
		expansion_counter += expansion_factor;
	}
	// finally, shuffle the results along column_axis
	shuffle_along_axis<float>(gr_pdfs, 2000, num_gr, 1);

	delete2DArray(input_psths);
	delete2DArray(input_pdf_templates);
	free(firing_categories);
	free(cell_spike_sums);
}

void PoissonRegenCells::calcGRPoissActivity(uint32_t ts)
{
	for (uint32_t i = 0; i < num_gr; i++)
	{
		threshs[i] += (threshMax - threshs[i]) * threshDec;
		aps[i] = randGens[0]->Random() < (gr_pdfs[ts][i] * threshs[i]);
		threshs[i] = aps[i] * threshBase + (aps[i] - 1) * threshs[i];
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
	write2DArray<uint8_t>(out_file, rasters, 2000, num_gr);
}

void PoissonRegenCells::save_psths(std::string out_file)
{
	write2DArray<uint8_t>(out_file, psths, 2000, num_gr);
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

