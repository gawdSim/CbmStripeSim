#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "file_parse.h"
#include "commandline.h"
#include "poissonregencells.h"
#include "mzonestate.h"
#include "cbmsimcore.h"

class Simulation
{
public:
	Simulation();
	Simulation(parsed_commandline &p_cl);
	~Simulation();

	// objects
	parsed_sess_file s_file;
	trials_data td;
	//PoissonRegenCells *grs = nullptr;
	MZoneState *sim_state  = nullptr;
	CBMSimCore *sim_core   = nullptr;

	bool data_out_dir_created     = false;
	bool out_sim_filename_created = false;
	bool trials_data_initialized  = false;
	bool sim_initialized          = false;

	std::string data_out_path;
	std::string data_out_base_name;
	std::string out_sim_name;

	uint8_t gpu_index = 0;
	uint8_t num_gpu   = 2;

	uint32_t trial;
	uint32_t rast_ctr;
	uint32_t psth_ctr;

	const uint32_t num_mzones = 1;
	const uint32_t pre_collect_ts = 2000;

	uint32_t trial_time;
	uint32_t ms_pre_cs;
	uint32_t ms_post_cs;
	uint32_t ms_measure;

	enum plasticity pf_pc_plast;

	void create_out_sim_filename();
	void build_sim();
	void init_sess(std::string sess_file);
	void set_plast_modes(std::string pfpc_plast);
	void init_sim(std::string in_psth_filename, std::string in_sim_filename);
	void run_session();

	void save_sim_to_file();

};

#endif /* SIMULATION_H_ */
