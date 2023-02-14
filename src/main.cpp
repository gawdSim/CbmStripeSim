#include <time.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <iomanip>

//#include "control.h"
#include "commandline.h"
#include "file_parse.h"
#include "logger.h"
#include "poissonregencells.h"

#define NUM_TRIALS 1000
#define NUM_TS 1300

int main(int argc, char **argv)
{
	logger_initConsoleLogger(stderr);
// for now, set the log level dependent on whether
// we are compiling for debug target or release target
#ifdef DEBUG
	logger_setLevel(LogLevel_DEBUG);
#else
	logger_setLevel(LogLevel_INFO);
#endif
	//parsed_commandline p_cl = {};
	//parse_and_validate_parsed_commandline(&argc, &argv, p_cl);

	if (argc < 2)
	{
		LOG_FATAL("Did not specify granule firing rate file. exiting..");
		exit(1);
	}
	std::string gr_fr_file = std::string(argv[1]);
	std::fstream gr_fr_file_buf(gr_fr_file.c_str(), std::ios::in | std::ios::binary);
	PoissonRegenCells gr_cells(0, gr_fr_file_buf);
	std::string out_rf_name, out_pf_name;
	double trial_start, trial_end;
	//omp_set_num_threads(1); /* for 4 gpus, 8 is the sweet spot. Unsure for 2. */
	for (uint32_t trial = 0; trial < NUM_TRIALS; trial++)
	{
		LOG_INFO("Trial number: %d", trial + 1);
		trial_start = omp_get_wtime();
		for (uint32_t ts = 0; ts < NUM_TS; ts++)
		{
			gr_cells.calcGRPoissActivity(ts);
			gr_cells.fill_rasters(ts);
			gr_cells.fill_psths(ts);
		}
		trial_end = omp_get_wtime();
		LOG_INFO("Trial %d took %0.2fs", trial + 1, trial_end - trial_start);
#ifdef DEBUG
		out_rf_name = "../../data/outputs/est_gr_fr_fourier_14_modes_GR_RASTER_02072023_trial_" + std::to_string(trial) + ".bin";
#else
		out_rf_name = "../data/outputs/est_gr_fr_fourier_14_modes_GR_RASTER_02072023_trial_" + std::to_string(trial) + ".bin";
#endif
		LOG_INFO("saving granule rasters to file..");
		gr_cells.save_rasters(out_rf_name);
	}
#ifdef DEBUG
	out_pf_name = "../../data/outputs/est_gr_fr_fourier_14_modes_GR_PSTH_02072023.bin";
#else
	out_pf_name = "../data/outputs/est_gr_fr_fourier_14_modes_GR_PSTH_02072023.bin";
#endif
	LOG_INFO("saving granule psths to file..");
	gr_cells.save_psths(out_pf_name);

	LOG_INFO("simulation finished. exiting..");
	return 0;
}

