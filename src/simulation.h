#ifndef SIMULATION_H_
#define SIMULATION_H_
#include <functional>

#include "file_parse.h"
#include "commandline.h"
#include "poissonregencells.h"
#include "mzonestate.h"
#include "cbmsimcore.h"

#define NUM_CELL_TYPES 6
#define NUM_WEIGHTS_TYPES 2

enum cell_id {GR, BC, SC, PC, IO, NC};

// should put these in file_utility.h
const std::string RAST_EXT[NUM_CELL_TYPES] = {".grr", ".bcr", ".scr", ".pcr", ".ior", ".ncr"}; 
const std::string PSTH_EXT[NUM_CELL_TYPES] = {".grp", ".bcp", ".scp", ".pcp", ".iop", ".ncp"}; 
const std::string WEIGHTS_EXT[NUM_WEIGHTS_TYPES] = {".pfpcw", ".mfncw"}; 

// convenience array for getting string representations of the cell ids
const std::string CELL_IDS[NUM_CELL_TYPES] = {"GR", "BC", "SC", "PC", "IO", "NC"}; 

enum datatype {RASTER, PSTH};

struct cell_spike_sums
{
	uint32_t non_cs_spike_sum;
	uint32_t cs_spike_sum;
	uint32_t *non_cs_spike_counter;
	uint32_t *cs_spike_counter;
};

struct cell_firing_rates
{
	float non_cs_mean_fr;
	float non_cs_median_fr;
	float cs_mean_fr;
	float cs_median_fr;
};

class Simulation
{
public:
	Simulation();
	Simulation(parsed_commandline &p_cl);
	~Simulation();

	// objects
	parsed_sess_file s_file;
	trials_data td;
	MZoneState *sim_state  = nullptr;
	CBMSimCore *sim_core   = nullptr;

	bool data_out_dir_created     = false;
	bool out_sim_filename_created = false;
	bool raster_filenames_created  = false;
	bool psth_filenames_created    = false;

	bool pfpc_weights_filenames_created = false;

	bool trials_data_initialized  = false;
	bool sim_initialized          = false;

	bool raster_arrays_initialized = false;
	bool psth_arrays_initialized   = false;
	bool spike_sums_initialized    = false;

	std::string data_out_path;
	std::string data_out_base_name;
	std::string out_sim_name;
	std::string rf_names[NUM_CELL_TYPES];
	std::string pf_names[NUM_CELL_TYPES]; 

	std::string pfpc_weights_file = "";
	std::string mfnc_weights_file = "";

	uint8_t gpu_index = 0;
	uint8_t num_gpu   = 4;

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

	struct cell_spike_sums spike_sums[NUM_CELL_TYPES];
	struct cell_firing_rates firing_rates[NUM_CELL_TYPES];

	const uint8_t *cell_spikes[NUM_CELL_TYPES];
	uint32_t rast_cell_nums[NUM_CELL_TYPES];
	uint8_t **rasters[NUM_CELL_TYPES];
	uint8_t **psths[NUM_CELL_TYPES];
	uint32_t rast_sizes[NUM_CELL_TYPES]; 
	std::function<void()> psth_save_funcs[NUM_CELL_TYPES];
	std::function<void()> rast_save_funcs[NUM_CELL_TYPES];

	// data collection
	// filename creation
	void create_out_sim_filename();
	void create_rast_or_psth_filenames(std::map<std::string, bool> &data_map, enum datatype data_type);
	void create_weights_filenames(std::map<std::string, bool> &weights_map);

	// data initialization
	void init_rast_cell_nums();
	void init_cell_spikes();
	void init_spike_sums();
	void init_rasts(); 
	void init_psths();
	void init_rast_save_funcs();
	void init_psth_save_funcs();

	// filling data
	void update_spike_sums(uint32_t tts, float onset_cs, float offset_cs);
	void calc_fire_rates(float onset_cs, float offset_cs);
	void fill_rasts(uint32_t rast_ctr, uint32_t psth_ctr);
	void fill_psths(uint32_t psth_ctr);

	// resetting data
	void reset_spike_sums();

	// saving data
	void save_sim();
	void save_gr_rast();
	void save_rasts();
	void save_psths();
	void save_pfpc_weights(int32_t trial = -1);

	// initializing in build or run mode, and building and running
	void build_sim();
	void init_sess(std::string sess_file);
	void set_plast_modes(std::string pfpc_plast);
	void init_sim(std::string in_psth_filename, std::string in_sim_filename);
	void run_session();

	// destructor things
	void delete_spike_sums();
	void delete_rasts();
	void delete_psths();
};

#endif /* SIMULATION_H_ */
