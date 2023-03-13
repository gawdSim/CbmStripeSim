/*
 * poissonregencells.cpp
 *
 *  Created on: Dec 13, 2012
 *      Author: consciousness
 *  Modified on: Jul 22, 2015
 *  	Author: evandelord
 */
#include <cstring>

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

	psth_sample_size = num_gr / expansion_factor;
	sample_template_indices = (size_t *)calloc(psth_sample_size, sizeof(size_t));
	template_indices = (size_t *)calloc(num_gr, sizeof(size_t));

	sample_indices = (size_t *)calloc(psth_sample_size, sizeof(size_t));

	psths = allocate2DArray<uint8_t>(2000, psth_sample_size);
	//rasters = allocate2DArray<uint8_t>(2000, num_gr); // hardcoded, 2000 ts by 1000 trials by num_gr grs
	gr_templates = allocate2DArray<float>(num_gr / expansion_factor, 2000); // for now, only have firing rates from trials of 2000
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
	free(template_indices);

	free(sample_indices);

	delete2DArray(psths);
	//delete2DArray(rasters);
	delete2DArray(gr_templates);
	delete2DArray(gr_templates_t);
	free(threshs);
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
	uint32_t num_fire_categories[3] = {0};
	// expect input data comes from 1000 trials, 2000 ts a piece.
	const float THRESH_FR = 20.0; // in Hz
	const uint32_t THRESH_COUNT = (THRESH_FR / 1000.0) * 2000 * 1000;
	uint8_t **input_psths = allocate2DArray<uint8_t>(2000, num_gr / expansion_factor);
	uint8_t *firing_categories = (uint8_t *)calloc(num_gr / expansion_factor, sizeof(uint8_t));
	uint32_t *cell_spike_sums = (uint32_t *)calloc(num_gr / expansion_factor, sizeof(uint32_t));
	LOG_INFO("loading cells from file.");
	input_psth_file_buf.read((char *)input_psths[0], 2000 * num_gr / expansion_factor * sizeof(uint8_t));
	LOG_INFO("finished load.");
	LOG_INFO("transposing...");
	uint8_t **input_psths_t = transpose2DArray(input_psths,  2000, num_gr / expansion_factor);
	LOG_INFO("finished transposing.");

	// determine firing rate categories
	LOG_INFO("determining firing categories...");
	for (size_t i = 0; i < num_gr / expansion_factor; i++)
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
	for (size_t i = 0; i < num_gr / expansion_factor; i++)
	{
		// NOTE: no need to include zero_firers as calloc sets their values to zero
		if (firing_categories[i] == LOW_FIRERS)
		{
			float mean_prob = (float)cell_spike_sums[i] / (2000 * 1000); // DEBUG
			for (size_t j = 0; j < 2000; j++)
			{
				gr_templates[i][j] = mean_prob; // just set every bin to same avg fr
			}
		}
		else // high firers
		{
			for (size_t j = 0; j < 2000; j++)
			{
				gr_templates[i][j] = (float)input_psths_t[i][j] / cell_spike_sums[i]; // the easiest transformation into a pdf
			}
		}
	}
	LOG_INFO("finished creating templates...");

	LOG_INFO("transposing templates...");
	gr_templates_t = transpose2DArray(gr_templates, num_gr / expansion_factor, 2000);
	LOG_INFO("finished transposing templates.");

	LOG_INFO("Generating template indices");
	size_t expansion_counter = 0;
	for (size_t i = 0; i < num_gr; i++)
	{
		template_indices[i] = expansion_counter;
		if (i % expansion_factor == 0 && i > 0)
		{
			expansion_counter++;
		}
	}
	LOG_INFO("finished generating template indices.");

	LOG_INFO("shuffling...");
	fisher_yates_shuffle<size_t>(template_indices, num_gr);
	LOG_INFO("finished shuffling...");

	LOG_INFO("generating random sample of non-unique template gr indices for psth saving...");
	// NOTE: the template indices are not guaranteed to be unique, as template_indices
	// contains expansion_factor multiples of each template id
	size_t curr_sample_index = 0;
	CRandomSFMT0 localRNG(time(0));
	uint8_t *chosen = (uint8_t *)calloc(num_gr, sizeof(uint8_t));
	while (curr_sample_index < psth_sample_size)
	{
		size_t test_index = localRNG.IRandom(0, num_gr-1);
		if (!chosen[test_index])
		{
			sample_template_indices[curr_sample_index] = template_indices[test_index];
			sample_indices[curr_sample_index] = test_index;
			curr_sample_index++;
			chosen[test_index] = 1;
		}
	}
	free(chosen);
	LOG_INFO("finished generation of random sample");

	delete2DArray(input_psths);
	delete2DArray(input_psths_t);
	free(firing_categories);
	free(cell_spike_sums);
}

void PoissonRegenCells::calcGRPoissActivity(uint32_t ts)
{
	for (uint32_t i = 0; i < num_gr; i++)
	{
		threshs[i] += (threshMax - threshs[i]) * threshInc;
		aps[i] = randGens[0]->Random() < (gr_templates_t[ts][template_indices[i]] * threshs[i]);
		threshs[i] = aps[i] * threshBase + (aps[i] - 1) * threshs[i];
	}
}

void PoissonRegenCells::calcGRPoissActivitySample(uint32_t ts)
{
	for (uint32_t i = 0; i < psth_sample_size; i++)
	{
		float temp_thresh = threshs[sample_indices[i]];
		uint8_t temp_ap = aps[sample_indices[i]];
		temp_thresh += (threshMax - temp_thresh) * threshInc;
		temp_ap = randGens[0]->Random() < (gr_templates_t[ts][sample_template_indices[i]] * temp_thresh);
		temp_thresh = temp_ap * threshBase + (temp_ap - 1) * temp_thresh;

		threshs[sample_indices[i]] = temp_thresh;
		aps[sample_indices[i]] = temp_ap;
	}
}

//void PoissonRegenCells::fill_rasters(uint32_t ts)
//{
//	memcpy(rasters[ts], aps, num_gr * sizeof(uint8_t));
//}

void PoissonRegenCells::fill_psths(size_t ts)
{
	for (size_t i = 0; i < num_gr; i++)
	{
		psths[ts][i] += aps[i];
	}
}

void PoissonRegenCells::fill_psths_sample(size_t ts)
{
	for (size_t i = 0; i < psth_sample_size; i++)
	{
		psths[ts][i] += aps[sample_indices[i]];
	}
}

void PoissonRegenCells::save_sample_template_indices(std::string out_file)
{
	std::fstream out_file_buf(out_file.c_str(), std::ios::out | std::ios::binary);

	if (!out_file_buf.is_open())
	{
		fprintf(stderr, "[ERROR]: Couldn't open '%s' for writing. Exiting...\n", out_file.c_str());
		exit(-1);
	}
	rawBytesRW((char *)sample_template_indices, psth_sample_size  * sizeof(size_t), false, out_file_buf);
	out_file_buf.close();
}

void PoissonRegenCells::save_template_indices(std::string out_file)
{
	std::fstream out_file_buf(out_file.c_str(), std::ios::out | std::ios::binary);

	if (!out_file_buf.is_open())
	{
		fprintf(stderr, "[ERROR]: Couldn't open '%s' for writing. Exiting...\n", out_file.c_str());
		exit(-1);
	}
	rawBytesRW((char *)template_indices, num_gr  * sizeof(size_t), false, out_file_buf);
	out_file_buf.close();
}

//void PoissonRegenCells::save_rasters(std::string out_file)
//{
//	write2DArray<uint8_t>(out_file, rasters, 2000, num_gr);
//}

void PoissonRegenCells::save_psths(std::string out_file)
{
	write2DArray<uint8_t>(out_file, psths, 2000, num_gr);
}

void PoissonRegenCells::save_psths_sample(std::string out_file)
{
	write2DArray<uint8_t>(out_file, psths, 2000, psth_sample_size);
}

const uint8_t *PoissonRegenCells::getGRAPs() { return (const uint8_t *)aps; }

const float **PoissonRegenCells::getGRFRs() { return (const float **)gr_fr; }

uint32_t **PoissonRegenCells::getApBufGR() { return apBufs; }

uint64_t **PoissonRegenCells::getApHistGR() { return apHists; }

