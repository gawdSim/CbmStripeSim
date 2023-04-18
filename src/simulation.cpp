#include <sys/stat.h> // mkdir (POSIX ONLY)
#include <omp.h> // omp_get_wtime

#include "logger.h"
#include "file_utility.h"
#include "dynamic2darray.h"
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
		create_rast_or_psth_filenames(p_cl.raster_files, RASTER); //optional
		create_rast_or_psth_filenames(p_cl.psth_files, PSTH); //optional
		create_weights_filenames(p_cl.weights_files); //optional
		init_sim(p_cl.input_psth_file, p_cl.input_sim_file);
	}
}

Simulation::~Simulation() {
	if (trials_data_initialized) delete_trials_data(td); // FIXME: this is silly and should be deprecated

	if (sim_state) delete sim_state;
	if (sim_core) delete sim_core;

	if (raster_arrays_initialized) delete_rasts();
	if (psth_arrays_initialized)   delete_psths();
	if (spike_sums_initialized)    delete_spike_sums();
}

void Simulation::create_out_sim_filename() {
	if (data_out_dir_created) {
		out_sim_name = data_out_path + "/" + data_out_base_name + SIM_EXT;
		out_sim_filename_created = true;
	}
}

void Simulation::create_rast_or_psth_filenames(std::map<std::string, bool> &data_map, enum datatype data_type)
{
	if (data_out_dir_created)
	{
		switch (data_type) {
		   case RASTER:
				for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
				{
					if (data_map[CELL_IDS[i]])
					{
						rf_names[i] = data_out_path + "/" + data_out_base_name + RAST_EXT[i];
					}
				}
				raster_filenames_created = true;
				break;
			case PSTH:
				for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
				{
					if (data_map[CELL_IDS[i]])
					{
						pf_names[i] = data_out_path + "/" + data_out_base_name + PSTH_EXT[i];
					}
				}
				psth_filenames_created = true;
				break;
			default:
			  LOG_DEBUG("Something has gone terribly wrong.");
			  exit(1);
			  // something has gone terribly wrong
		}
	}
}

void Simulation::create_weights_filenames(std::map<std::string, bool> &weights_map) 
{
	if (data_out_dir_created)
	{
		if (weights_map["PFPC"])
		{
			pfpc_weights_file = data_out_path + "/" + data_out_base_name + ".pfpcw";
			pfpc_weights_filenames_created = true; // only useful so far for gui...
		}
	}
}

void Simulation::init_rast_cell_nums()
{
	rast_cell_nums[GR] = num_gr; 
	rast_cell_nums[BC] = num_bc;
	rast_cell_nums[SC] = num_sc;
	rast_cell_nums[PC] = num_pc;
	rast_cell_nums[IO] = num_io;
	rast_cell_nums[NC] = num_nc;
}

void Simulation::init_cell_spikes()
{
	for (uint32_t i; i < num_mzones; i++) {
		/* NOTE: incurs a call to cudaMemcpy from device to host,
		 * but initializing so is not repeatedly called */
		cell_spikes[GR] = sim_core->getGRRRRRRRs()->getGRAPs(); 
		cell_spikes[BC] = sim_core->getMZoneList()[i]->exportAPBC(); 
		cell_spikes[SC] = sim_core->getMZoneList()[i]->exportAPSC();
		cell_spikes[PC] = sim_core->getMZoneList()[i]->exportAPPC();
		cell_spikes[IO] = sim_core->getMZoneList()[i]->exportAPIO();
		cell_spikes[NC] = sim_core->getMZoneList()[i]->exportAPNC();
	}
}

void Simulation::init_spike_sums()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		spike_sums[i].non_cs_spike_sum = 0;
		spike_sums[i].cs_spike_sum = 0;
		spike_sums[i].non_cs_spike_counter = (uint32_t *)calloc(rast_cell_nums[i], sizeof(uint32_t));
		spike_sums[i].cs_spike_counter = (uint32_t *)calloc(rast_cell_nums[i], sizeof(uint32_t));
	}
	spike_sums_initialized = true;

}

void Simulation::init_rasts()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!rf_names[i].empty())
		{
			/* granules are saved every trial, so their raster size is msMeasure  x num_gr */
			uint32_t row_size = (CELL_IDS[i] == "GR") ? ms_measure : ms_measure * td.num_trials;
			rasters[i] = allocate2DArray<uint8_t>(row_size, rast_cell_nums[i]);
		}
	}
	raster_arrays_initialized = true;
}

void Simulation::init_psths()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!pf_names[i].empty())
			// TODO: make data type bigger for psth
			// FIXME: this guy will be too big for 2^24 gr cell!!!
			psths[i] = allocate2DArray<uint8_t>(ms_measure, rast_cell_nums[i]);
	}
	psth_arrays_initialized = true;
}

void Simulation::init_rast_save_funcs()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		rast_save_funcs[i] = [this, i]()
		{
			if (!rf_names[i].empty() && CELL_IDS[i] != "GR")
			{
				uint32_t row_size = (CELL_IDS[i] == "GR") ? this->ms_measure : this->ms_measure * this->td.num_trials;
				LOG_DEBUG("Saving %s raster to file...", CELL_IDS[i].c_str());
				write2DArray<uint8_t>(rf_names[i], this->rasters[i], row_size, this->rast_cell_nums[i]);
			}
		};
	}
}


void Simulation::init_psth_save_funcs()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		psth_save_funcs[i] = [this, i]()
		{
			if (!pf_names[i].empty())
			{
				LOG_DEBUG("Saving %s psth to file...", CELL_IDS[i].c_str());
				write2DArray<uint8_t>(pf_names[i], this->psths[i], this->ms_measure, this->rast_cell_nums[i]);
			}
		};
	}
}

void Simulation::update_spike_sums(uint32_t tts, float onset_cs, float offset_cs)
{
	// update cs spikes
	if (tts >= onset_cs && tts < offset_cs)
	{
		for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
		{
			for (uint32_t j = 0; j < rast_cell_nums[i]; j++)
			{
				spike_sums[i].cs_spike_sum += cell_spikes[i][j];
				spike_sums[i].cs_spike_counter[j] += cell_spikes[i][j];
			}
		}
	}
	// update non-cs spikes
	else if (tts < onset_cs)
	{
		for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
		{
			for (uint32_t j = 0; j < rast_cell_nums[i]; j++)
			{
				spike_sums[i].non_cs_spike_sum += cell_spikes[i][j];
				spike_sums[i].non_cs_spike_counter[j] += cell_spikes[i][j];
			}
		}
	}
}

void Simulation::calc_fire_rates(float onset_cs, float offset_cs)
{
	float non_cs_time_secs = (onset_cs - 1) / 1000.0; // why only pre-cs? (Ask Joe)
	float cs_time_secs = (offset_cs - onset_cs) / 1000.0;

	for (int i = 0; i < NUM_CELL_TYPES; i++)
	{
		// sort sums for medians 
		std::sort(spike_sums[i].cs_spike_counter,
			spike_sums[i].cs_spike_counter + rast_cell_nums[i]);
		std::sort(spike_sums[i].non_cs_spike_counter,
			spike_sums[i].non_cs_spike_counter + rast_cell_nums[i]);
		
		// calculate medians
		firing_rates[i].non_cs_median_fr =
			(spike_sums[i].non_cs_spike_counter[rast_cell_nums[i] / 2 - 1]
		   + spike_sums[i].non_cs_spike_counter[rast_cell_nums[i] / 2]) / (2.0 * non_cs_time_secs);
		firing_rates[i].cs_median_fr     =
			(spike_sums[i].cs_spike_counter[rast_cell_nums[i] / 2 - 1]
		   + spike_sums[i].cs_spike_counter[rast_cell_nums[i] / 2]) / (2.0 * cs_time_secs);
		
		// calculate means
		firing_rates[i].non_cs_mean_fr = spike_sums[i].non_cs_spike_sum / (non_cs_time_secs * rast_cell_nums[i]);
		firing_rates[i].cs_mean_fr     = spike_sums[i].cs_spike_sum / (cs_time_secs * rast_cell_nums[i]);
	}
}

void Simulation::fill_rasts(uint32_t rast_ctr, uint32_t psth_ctr)
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		uint32_t temp_counter = rast_ctr;
		if (!rf_names[i].empty())
		{
			/* GR spikes are only spikes not saved on host every time step:
			 * InNet::exportAPGR makes cudaMemcpy call before returning pointer to mem address */
			if (CELL_IDS[i] == "GR")
			{
				cell_spikes[i] = sim_core->getGRRRRRRRs()->getGRAPs(); 
				temp_counter = psth_ctr;
			}
			for (uint32_t j = 0; j < rast_cell_nums[i]; j++)
			{
				rasters[i][temp_counter][j] = cell_spikes[i][j];
			}
		}
	}
}

void Simulation::fill_psths(uint32_t psth_ctr)
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!pf_names[i].empty())
		{
			/* GR spikes are only spikes not saved on host every time step:
			 * InNet::exportAPGR makes cudaMemcpy call before returning pointer to mem address */
			if (CELL_IDS[i] == "GR")
			{
				cell_spikes[i] = sim_core->getGRRRRRRRs()->getGRAPs(); 
			}
			for (uint32_t j = 0; j < rast_cell_nums[i]; j++)
			{
				psths[i][psth_ctr][j] += cell_spikes[i][j];
			}
		}
	}
}

void Simulation::reset_spike_sums()
{
		for (int i = 0; i < NUM_CELL_TYPES; i++)
		{
			spike_sums[i].cs_spike_sum = 0;
			spike_sums[i].non_cs_spike_sum = 0;
			memset(spike_sums[i].cs_spike_counter, 0, rast_cell_nums[i] * sizeof(uint32_t));
			memset(spike_sums[i].non_cs_spike_counter, 0, rast_cell_nums[i] * sizeof(uint32_t));
		}
}

void Simulation::save_sim()
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

void Simulation::save_gr_rast()
{
	if (!rf_names[GR].empty())
	{
		std::string trial_raster_name = data_out_path + "/" + get_file_basename(rf_names[GR])
									  + "_trial_" + std::to_string(trial) + BIN_EXT;
		LOG_DEBUG("Saving granule raster to file...");
		write2DArray<uint8_t>(trial_raster_name, rasters[GR], num_gr, ms_measure);
	}
}

void Simulation::save_rasts()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!rf_names[i].empty())
			rast_save_funcs[i]();
	}
}

void Simulation::save_psths()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!pf_names[i].empty())
			psth_save_funcs[i]();
	}
}

void Simulation::save_pfpc_weights(int32_t trial)
{
	if (pfpc_weights_filenames_created)
	{
		if (trial != -1) // nonnegative indicates we want to append the trial to the file basename
		{
			pfpc_weights_file = data_out_path + "/" + get_file_basename(pfpc_weights_file)
									   + "_TRIAL_" + std::to_string(trial) + BIN_EXT;
		}
		LOG_DEBUG("Saving granule to purkinje weights to file...");
		if (!sim_core)
		{
			LOG_ERROR("Trying to write uninitialized weights to file.");
			LOG_ERROR("(Hint: Try initializing a sim or loading the weights first.)");
			return;
		}
		const float *pfpc_weights = sim_core->getMZoneList()[0]->exportPFPCWeights();
		std::fstream outPFPCFileBuffer(pfpc_weights_file.c_str(), std::ios::out | std::ios::binary);
		rawBytesRW((char *)pfpc_weights, num_gr * sizeof(float), false, outPFPCFileBuffer);
		outPFPCFileBuffer.close();
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

	trial_time = std::stoi(s_file.parsed_var_sections["trial_spec"].param_map["trialTime"]);
	ms_pre_cs  = std::stoi(s_file.parsed_var_sections["trial_spec"].param_map["msPreCS"]);
	ms_post_cs = std::stoi(s_file.parsed_var_sections["trial_spec"].param_map["msPostCS"]);
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
	std::fstream sim_file_buf(in_sim_filename.c_str(), std::ios::in | std::ios::binary);
	read_con_params(sim_file_buf);
	populate_act_params(s_file);
	sim_state = new MZoneState(num_mzones, sim_file_buf);
	sim_file_buf.close();
	std::fstream psth_file_buf(in_psth_filename.c_str(), std::ios::in | std::ios::binary);
	sim_core  = new CBMSimCore(psth_file_buf, sim_state, gpu_index, num_gpu);
	psth_file_buf.close();
	init_rast_cell_nums();
	init_cell_spikes();
	init_spike_sums();
	init_rasts(); 
	init_psths();
	init_rast_save_funcs();
	init_psth_save_funcs();
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
			//update_spike_sums(ts, onsetCS, onsetCS + csLength);
			if (ts >= onsetCS - ms_pre_cs && ts < onsetCS + csLength + ms_post_cs)
			{
				fill_rasts(rast_ctr, psth_ctr);
				fill_psths(psth_ctr);
				psth_ctr++;
				rast_ctr++;
			}
		}
		//calc_fire_rates(onsetCS, onsetCS + csLength);
		//LOG_INFO("PC Pop mean CS fr: %.2f", firing_rates[PC].cs_mean_fr);
		//LOG_INFO("BC Pop mean CS fr: %.2f", firing_rates[BC].cs_mean_fr);
		//LOG_INFO("SC Pop mean CS fr: %.2f", firing_rates[SC].cs_mean_fr);
		//LOG_INFO("NC Pop mean CS fr: %.2f", firing_rates[NC].cs_mean_fr);
		//LOG_INFO("IO Pop mean CS fr: %.2f", firing_rates[IO].cs_mean_fr);
		//reset_spike_sums();
		save_pfpc_weights(trial);
		trial_end = omp_get_wtime();
		LOG_INFO("'%s' took %0.2fs", trialName.c_str(), trial_end - trial_start);
		trial++;
	}
	save_sim();
	save_rasts();
	save_psths();
	session_end = omp_get_wtime();
	LOG_INFO("Session finished. took %0.2fs", session_end - session_start);
}

void Simulation::delete_spike_sums()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		free(spike_sums[i].non_cs_spike_counter);
		free(spike_sums[i].cs_spike_counter);
	}
}

void Simulation::delete_rasts()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!rf_names[i].empty()) delete2DArray<uint8_t>(rasters[i]);
	}
}

void Simulation::delete_psths()
{
	for (uint32_t i = 0; i < NUM_CELL_TYPES; i++)
	{
		if (!pf_names[i].empty()) delete2DArray<uint8_t>(psths[i]);
	}
}

