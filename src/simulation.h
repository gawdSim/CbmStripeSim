#ifndef SIMULATION_H_
#define SIMULATION_H_
#include <functional>

#include "file_parse.h"
#include "commandline.h"
#include "poissonregencells.h"
#include "mzonestate.h"
#include "cbmsimcore.h"

#define NUM_CELL_TYPES 8
#define NUM_WEIGHTS_TYPES 2

enum cell_id {MF, GR, GO, BC, SC, PC, IO, NC};

// convenience array for getting string representations of the cell ids
const std::string CELL_IDS[NUM_CELL_TYPES] = {"MF", "GR", "GO", "BC", "SC", "PC", "IO", "NC"}; 

enum datatype {RASTER, PSTH};

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
	bool mfnc_weights_filenames_created = false;

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
	void initialize_rast_cell_nums();
	void initialize_cell_spikes();
	void initialize_spike_sums();
	void initialize_rasters(); 
	void initialize_psth_save_funcs();
	void initialize_raster_save_funcs();

	void initialize_psths();

	// filling data
	void fill_rasters(uint32_t raster_counter, uint32_t psth_counter);
	void fill_psths(uint32_t psth_counter);

	// saving data
	void save_sim();
	void save_weights();
	void save_gr_raster();
	void save_rasters();
	void save_psths();

	// initializing in build or run mode, and building and running
	void build_sim();
	void init_sess(std::string sess_file);
	void set_plast_modes(std::string pfpc_plast);
	void init_sim(std::string in_psth_filename, std::string in_sim_filename);
	void run_session();
};

#endif /* SIMULATION_H_ */
