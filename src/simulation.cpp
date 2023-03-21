#include <sys/stat.h> // mkdir (POSIX ONLY)
#include <omp.h> // omp_get_wtime

#include "logger.h"
#include "file_utility.h"
#include "connectivityparams.h"
#include "activityparams.h"
#include "simulation.h"

Simulation::Simulation() {}

Simulation::Simulation(parsed_commandline &p_cl) {
	if (!p_cl.build_file.empty()) {
		tokenized_file t_file;
		lexed_file l_file;
		parsed_build_file pb_file;
		tokenize_file(p_cl.build_file, t_file);
		lex_tokenized_file(t_file, l_file);
		parse_lexed_build_file(l_file, pb_file);
		if (!con_params_populated) populate_con_params(pb_file);
		data_out_path = OUTPUT_DATA_PATH + p_cl.output_basename;
		data_out_base_name = p_cl.output_basename;
		int status = mkdir(data_out_path.c_str(), 0775);
		if (status == -1) {
			LOG_FATAL("Could not create directory '%s'. Maybe it already exists. Exiting...", data_out_path.c_str());
			exit(10);
		}
		data_out_dir_created = true;
		create_out_sim_filename(); //default
	} else if (!p_cl.session_file.empty()) {
		init_sess(p_cl.session_file);
		set_plast_modes(p_cl.pfpc_plasticity);
		// assume that validated commandline opts includes 1) input file 2) session file 3) output directory name
		data_out_path = OUTPUT_DATA_PATH + p_cl.output_basename;
		data_out_base_name = p_cl.output_basename;
		// NOTE: make the output directory here, so in case of error, user not
		// run an entire simulation just to not have files save
		int status = mkdir(data_out_path.c_str(), 0775);
		if (status == -1) {
			LOG_FATAL("Could not create directory '%s'. Maybe it already exists. Exiting...", data_out_path.c_str());
			exit(10);
		}
		data_out_dir_created = true;
		create_out_sim_filename(); 
		init_sim(p_cl.input_psth_file, p_cl.input_sim_file);
	}
}

Simulation::~Simulation() {
	if (trials_data_initialized) delete_trials_data(td); // FIXME: this is silly and should be deprecated

	//if (grs) delete grs;
	if (sim_state) delete sim_state;
	if (sim_core) delete sim_core;
}

void Simulation::create_out_sim_filename() {
	if (data_out_dir_created) {
		out_sim_name = data_out_path + "/"
					 + data_out_base_name + "_"
					 + get_current_time_as_string("%m%d%Y")
					 + SIM_EXT;
		out_sim_filename_created = true;
	}
}


void Simulation::build_sim() {
	if (!sim_state) sim_state = new MZoneState(num_mzones);
}

void Simulation::init_sess(std::string sess_file) {
	LOG_DEBUG("Initializing session...");
	tokenized_file t_file;
	lexed_file l_file;
	tokenize_file(sess_file, t_file);
	lex_tokenized_file(t_file, l_file);
	parse_lexed_sess_file(l_file, s_file);
	translate_parsed_trials(s_file, td);

	trial_time = std::stoi(s_file.parsed_var_sections["trial_spec"].param_map["trialTime"].value);
	ms_pre_cs  = std::stoi(s_file.parsed_var_sections["trial_spec"].param_map["msPreCS"].value);
	ms_post_cs = std::stoi(s_file.parsed_var_sections["trial_spec"].param_map["msPostCS"].value);
	ms_measure = ms_pre_cs + td.cs_lens[0] + ms_post_cs; // assumes all trials have same cs len

	trials_data_initialized = true;
	LOG_DEBUG("Session initialized.");
}

void Simulation::set_plast_modes(std::string pfpc_plast) {
	if (pfpc_plast == "off") this->pf_pc_plast = OFF;
	else if (pfpc_plast == "graded") this->pf_pc_plast = GRADED;
	else if (pfpc_plast == "binary") this->pf_pc_plast = BINARY;
	else if (pfpc_plast == "abbott-cascade") this->pf_pc_plast = ABBOTT_CASCADE;
	else if (pfpc_plast == "mauk-cascade") this->pf_pc_plast = MAUK_CASCADE;
}

void Simulation::init_sim(std::string in_psth_filename, std::string in_sim_filename) {
	LOG_DEBUG("Initializing simulation...");
	//grs = new PoissonRegenCells(psth_file_buf);
	std::fstream sim_file_buf(in_sim_filename.c_str(), std::ios::in | std::ios::binary);
	read_con_params(sim_file_buf);
	populate_act_params(s_file);
	sim_state = new MZoneState(num_mzones, sim_file_buf);
	sim_file_buf.close();
	std::fstream psth_file_buf(in_psth_filename.c_str(), std::ios::in | std::ios::binary);
	sim_core  = new CBMSimCore(psth_file_buf, sim_state, gpu_index, num_gpu);
	psth_file_buf.close();
	sim_initialized = true;
	LOG_DEBUG("Simulation initialized.");
}

void Simulation::run_session() {

	double session_start, session_end, trial_start, trial_end;
	trial = 0;
	rast_ctr = 0;
	session_start = omp_get_wtime();
	while (trial < td.num_trials) {
		std::string trialName = td.trial_names[trial];

		uint32_t useCS        = td.use_css[trial];
		uint32_t onsetCS      = pre_collect_ts + td.cs_onsets[trial];
		uint32_t csLength     = td.cs_lens[trial];
		uint32_t useUS        = td.use_uss[trial];
		uint32_t onsetUS      = pre_collect_ts + td.us_onsets[trial];

		psth_ctr = 0;
		LOG_INFO("Trial number: %d", trial + 1);
		trial_start = omp_get_wtime();
		for (uint32_t ts = 0; ts < trial_time; ts++) {
			/* deliver the US */
			if (useUS == 1 && ts == onsetUS) {
				sim_core->updateErrDrive(0, 0.3);
			}
			sim_core->calcActivity(pf_pc_plast, ts);
		}
		trial_end = omp_get_wtime();
		LOG_INFO("'%s' took %0.2fs", trialName.c_str(), trial_end - trial_start);
		trial++;
	}
	session_end = omp_get_wtime();
	LOG_INFO("Session finished. took %0.2fs", session_end - session_start);
}

void Simulation::save_sim_to_file()
{
	if (out_sim_filename_created)
	{
		LOG_DEBUG("Saving simulation to file...");
		std::fstream outSimFileBuffer(out_sim_name.c_str(), std::ios::out | std::ios::binary);
		write_con_params(outSimFileBuffer);
		if (!sim_core) sim_state->writeState(outSimFileBuffer);
		else sim_core->writeState(outSimFileBuffer);
		outSimFileBuffer.close();
		LOG_DEBUG("simulation save to file %s.", out_sim_name.c_str());
	}
}

