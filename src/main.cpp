#include <time.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <iomanip>

#include "commandline.h"
#include "file_parse.h"
#include "logger.h"
#include "poissonregencells.h"
#include "simulation.h"

#define NUM_TRIALS 1000
#define NUM_TS 2800

#ifdef DEBUG
#define DATA_OUT_DIR "../../data/outputs/"
#else
#define DATA_OUT_DIR "../data/outputs/"
#endif

#define OUT_FILE_BASENAME "please_just_work"
#define BIN_EXT ".bin"

int main(int argc, char **argv)
{
	logger_initConsoleLogger(stderr);
// for now, set the log level dependent on whether
// we are compiling for debug target or release target
//#ifdef DEBUG
//	logger_setLevel(LogLevel_DEBUG);
//#else
//	logger_setLevel(LogLevel_INFO);
//#endif

	logger_setLevel(LogLevel_DEBUG);
	parsed_commandline p_cl = {};
	parse_and_validate_parsed_commandline(&argc, &argv, p_cl);
	Simulation stripe_sim(p_cl);

	if (p_cl.vis_mode == "TUI") {
		if (!p_cl.build_file.empty()) {
			stripe_sim.build_sim();
			stripe_sim.save_sim();
		} else if (!p_cl.session_file.empty()) {
			stripe_sim.run_session();
		}
	}
	else if (p_cl.vis_mode == "GUI") {
		LOG_DEBUG("GUI not implemented for stripe simulation. Exiting...");
	}
	return 0;
}

