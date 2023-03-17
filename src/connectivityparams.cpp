/*
 * connectivityparams.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: varicella
 */

#include "connectivityparams.h"

bool con_params_populated = false;

uint64_t gr_x                    = 0;//16 * 2048; 
uint64_t gr_y                    = 0;//8 * 512; 
uint64_t num_gr                  = 1 * 1048576;//128 * 1048576; 
int go_x                         = 0; 
int go_y                         = 0; 
int num_go                       = 0; 
int num_bc                       = 0; 
int num_sc                       = 0; 
int num_pc                       = 0; 
int num_nc                       = 0; 
int num_io                       = 0; 
int gr_pf_vel_in_gr_x_per_t_step = 0; 
int gr_af_delay_in_t_step        = 0; 
int num_p_bc_from_bc_to_pc       = 0; 
int num_p_pc_from_bc_to_pc       = 0; 
int num_p_bc_from_gr_to_bc       = 0; 
int num_p_bc_from_gr_to_bc_p2    = 0; 
int num_p_pc_from_pc_to_bc       = 0; 
int num_p_bc_from_pc_to_bc       = 0; 
int num_p_sc_from_sc_to_pc       = 0; 
int num_p_pc_from_sc_to_pc       = 0; 
int num_p_sc_from_gr_to_sc       = 0; 
int num_p_sc_from_gr_to_sc_p2    = 0; 
int num_p_pc_from_pc_to_nc       = 0; 
int num_p_nc_from_pc_to_nc       = 0; 
int num_p_pc_from_gr_to_pc       = 0; 
int num_p_pc_from_gr_to_pc_p2    = 0; 
int num_p_mf_from_mf_to_nc       = 0; 
int num_p_nc_from_mf_to_nc       = 0; 
int num_p_nc_from_nc_to_io       = 0; 
int num_p_io_from_nc_to_io       = 0; 
int num_p_io_from_io_to_pc       = 0; 
int num_p_io_in_io_to_io         = 0; 
int num_p_io_out_io_to_io        = 0; 

float msPerTimeStep            = 1.0;
float numPopHistBinsPC         = 0.0; 

//float eLeakGR      = 0.0; 
//float threshRestGR = 0.0;

float eLeakSC          = 0.0;
float threshRestSC     = 0.0;
float eLeakBC          = 0.0;
float threshRestBC     = 0.0;
float eLeakPC          = 0.0;
float threshRestPC     = 0.0;
float initSynWofGRtoPC = 0.0;
float eLeakIO          = 0.0;
float threshRestIO     = 0.0;
float eLeakNC          = 0.0;
float threshRestNC     = 0.0;

void populate_con_params(parsed_build_file &p_file)
{
	/* int con params */
	gr_x                         = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["gr_x"].value); 
	gr_y                         = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["gr_y"].value); 
	num_gr                       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_gr"].value); 
	num_bc                       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_bc"].value); 
	num_sc                       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_sc"].value);
	num_pc                       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_pc"].value); 
	num_nc                       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_nc"].value); 
	num_io                       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_io"].value); 
	gr_pf_vel_in_gr_x_per_t_step = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["gr_pf_vel_in_gr_x_per_t_step"].value);
	gr_af_delay_in_t_step        = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["gr_af_delay_in_t_step"].value);
	num_p_bc_from_bc_to_pc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_bc_from_bc_to_pc"].value);
	num_p_pc_from_bc_to_pc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_pc_from_bc_to_pc"].value);
	num_p_bc_from_gr_to_bc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_bc_from_gr_to_bc"].value);
	num_p_bc_from_gr_to_bc_p2    = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_bc_from_gr_to_bc_p2"].value);
	num_p_pc_from_pc_to_bc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_pc_from_pc_to_bc"].value);
	num_p_bc_from_pc_to_bc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_bc_from_pc_to_bc"].value);
	num_p_sc_from_sc_to_pc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_sc_from_sc_to_pc"].value);
	num_p_pc_from_sc_to_pc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_pc_from_sc_to_pc"].value);
	num_p_sc_from_gr_to_sc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_sc_from_gr_to_sc"].value);
	num_p_sc_from_gr_to_sc_p2    = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_sc_from_gr_to_sc_p2"].value);
	num_p_pc_from_pc_to_nc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_pc_from_pc_to_nc"].value);
	num_p_nc_from_pc_to_nc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_nc_from_pc_to_nc"].value);
	num_p_pc_from_gr_to_pc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_pc_from_gr_to_pc"].value);
	num_p_pc_from_gr_to_pc_p2    = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_pc_from_gr_to_pc_p2"].value);
	num_p_mf_from_mf_to_nc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_mf_from_mf_to_nc"].value);
	num_p_nc_from_mf_to_nc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_nc_from_mf_to_nc"].value);
	num_p_nc_from_nc_to_io       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_nc_from_nc_to_io"].value);
	num_p_io_from_nc_to_io       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_io_from_nc_to_io"].value);
	num_p_io_from_io_to_pc       = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_io_from_io_to_pc"].value);
	num_p_io_in_io_to_io         = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_io_in_io_to_io"].value);
	num_p_io_out_io_to_io        = std::stoi(p_file.parsed_var_sections["connectivity"].param_map["num_p_io_out_io_to_io"].value);

	/* float con params */
	msPerTimeStep            = std::stof(p_file.parsed_var_sections["connectivity"].param_map["msPerTimeStep"].value); 
	numPopHistBinsPC         = std::stof(p_file.parsed_var_sections["connectivity"].param_map["numPopHistBinsPC"].value); 

	/* act params */
	//eLeakGR      = std::stof(p_file.parsed_var_sections["activity"].param_map["eLeakGR"].value); 
	//threshRestGR = std::stof(p_file.parsed_var_sections["activity"].param_map["threshRestGR"].value); 
	
	eLeakSC          = std::stof(p_file.parsed_var_sections["activity"].param_map["eLeakSC"].value);
	threshRestSC     = std::stof(p_file.parsed_var_sections["activity"].param_map["threshRestSC"].value);
	eLeakBC          = std::stof(p_file.parsed_var_sections["activity"].param_map["eLeakBC"].value);
	threshRestBC     = std::stof(p_file.parsed_var_sections["activity"].param_map["threshRestBC"].value);
	eLeakPC          = std::stof(p_file.parsed_var_sections["activity"].param_map["eLeakPC"].value);
	threshRestPC     = std::stof(p_file.parsed_var_sections["activity"].param_map["threshRestPC"].value);
	initSynWofGRtoPC = std::stof(p_file.parsed_var_sections["activity"].param_map["initSynWofGRtoPC"].value);
	eLeakIO          = std::stof(p_file.parsed_var_sections["activity"].param_map["eLeakIO"].value);
	threshRestIO     = std::stof(p_file.parsed_var_sections["activity"].param_map["threshRestIO"].value);
	eLeakNC          = std::stof(p_file.parsed_var_sections["activity"].param_map["eLeakNC"].value);
	threshRestNC     = std::stof(p_file.parsed_var_sections["activity"].param_map["threshRestNC"].value);

	con_params_populated = true;
}

void read_con_params(std::fstream &in_param_buf)
{
	/* not checking whether these things are zeros or not... */
	in_param_buf.read((char *)&gr_x, sizeof(int));
	in_param_buf.read((char *)&gr_y, sizeof(int));
	in_param_buf.read((char *)&num_gr, sizeof(int));
	in_param_buf.read((char *)&num_bc, sizeof(int));
	in_param_buf.read((char *)&num_sc, sizeof(int));
	in_param_buf.read((char *)&num_pc, sizeof(int));
	in_param_buf.read((char *)&num_nc, sizeof(int));
	in_param_buf.read((char *)&num_io, sizeof(int));
	in_param_buf.read((char *)&gr_pf_vel_in_gr_x_per_t_step, sizeof(int));
	in_param_buf.read((char *)&gr_af_delay_in_t_step, sizeof(int));
	in_param_buf.read((char *)&num_p_bc_from_bc_to_pc, sizeof(int));
	in_param_buf.read((char *)&num_p_pc_from_bc_to_pc, sizeof(int));
	in_param_buf.read((char *)&num_p_bc_from_gr_to_bc, sizeof(int));
	in_param_buf.read((char *)&num_p_bc_from_gr_to_bc_p2, sizeof(int));
	in_param_buf.read((char *)&num_p_pc_from_pc_to_bc, sizeof(int));
	in_param_buf.read((char *)&num_p_bc_from_pc_to_bc, sizeof(int));
	in_param_buf.read((char *)&num_p_sc_from_sc_to_pc, sizeof(int));
	in_param_buf.read((char *)&num_p_pc_from_sc_to_pc, sizeof(int));
	in_param_buf.read((char *)&num_p_sc_from_gr_to_sc, sizeof(int));
	in_param_buf.read((char *)&num_p_sc_from_gr_to_sc_p2, sizeof(int));
	in_param_buf.read((char *)&num_p_pc_from_pc_to_nc, sizeof(int));
	in_param_buf.read((char *)&num_p_nc_from_pc_to_nc, sizeof(int));
	in_param_buf.read((char *)&num_p_pc_from_gr_to_pc, sizeof(int));
	in_param_buf.read((char *)&num_p_pc_from_gr_to_pc_p2, sizeof(int));
	in_param_buf.read((char *)&num_p_mf_from_mf_to_nc, sizeof(int));
	in_param_buf.read((char *)&num_p_nc_from_mf_to_nc, sizeof(int));
	in_param_buf.read((char *)&num_p_nc_from_nc_to_io, sizeof(int));
	in_param_buf.read((char *)&num_p_io_from_nc_to_io, sizeof(int));
	in_param_buf.read((char *)&num_p_io_from_io_to_pc, sizeof(int));
	in_param_buf.read((char *)&num_p_io_in_io_to_io, sizeof(int));
	in_param_buf.read((char *)&num_p_io_out_io_to_io, sizeof(int));

	in_param_buf.read((char *)&msPerTimeStep, sizeof(float));
	in_param_buf.read((char *)&numPopHistBinsPC, sizeof(float));

	/* act params */
	//in_param_buf.read((char *)&eLeakGR, sizeof(float)); 
	//in_param_buf.read((char *)&threshRestGR, sizeof(float));

	in_param_buf.read((char *)&eLeakSC, sizeof(float));
	in_param_buf.read((char *)&threshRestSC, sizeof(float));
	in_param_buf.read((char *)&eLeakBC, sizeof(float));
	in_param_buf.read((char *)&threshRestBC, sizeof(float));
	in_param_buf.read((char *)&eLeakPC, sizeof(float));
	in_param_buf.read((char *)&threshRestPC, sizeof(float));
	in_param_buf.read((char *)&initSynWofGRtoPC, sizeof(float));
	in_param_buf.read((char *)&eLeakIO, sizeof(float));
	in_param_buf.read((char *)&threshRestIO, sizeof(float));
	in_param_buf.read((char *)&eLeakNC, sizeof(float));
	in_param_buf.read((char *)&threshRestNC, sizeof(float));

	con_params_populated = true;
}

void write_con_params(std::fstream &out_param_buf)
{
	/* not checking whether these things are zeros or not... */
	out_param_buf.write((char *)&gr_x, sizeof(int));
	out_param_buf.write((char *)&gr_y, sizeof(int));
	out_param_buf.write((char *)&num_gr, sizeof(int));
	out_param_buf.write((char *)&num_bc, sizeof(int));
	out_param_buf.write((char *)&num_sc, sizeof(int));
	out_param_buf.write((char *)&num_pc, sizeof(int));
	out_param_buf.write((char *)&num_nc, sizeof(int));
	out_param_buf.write((char *)&num_io, sizeof(int));
	out_param_buf.write((char *)&gr_pf_vel_in_gr_x_per_t_step, sizeof(int));
	out_param_buf.write((char *)&gr_af_delay_in_t_step, sizeof(int));
	out_param_buf.write((char *)&num_p_bc_from_bc_to_pc, sizeof(int));
	out_param_buf.write((char *)&num_p_pc_from_bc_to_pc, sizeof(int));
	out_param_buf.write((char *)&num_p_bc_from_gr_to_bc, sizeof(int));
	out_param_buf.write((char *)&num_p_bc_from_gr_to_bc_p2, sizeof(int));
	out_param_buf.write((char *)&num_p_pc_from_pc_to_bc, sizeof(int));
	out_param_buf.write((char *)&num_p_bc_from_pc_to_bc, sizeof(int));
	out_param_buf.write((char *)&num_p_sc_from_sc_to_pc, sizeof(int));
	out_param_buf.write((char *)&num_p_pc_from_sc_to_pc, sizeof(int));
	out_param_buf.write((char *)&num_p_sc_from_gr_to_sc, sizeof(int));
	out_param_buf.write((char *)&num_p_sc_from_gr_to_sc_p2, sizeof(int));
	out_param_buf.write((char *)&num_p_pc_from_pc_to_nc, sizeof(int));
	out_param_buf.write((char *)&num_p_nc_from_pc_to_nc, sizeof(int));
	out_param_buf.write((char *)&num_p_pc_from_gr_to_pc, sizeof(int));
	out_param_buf.write((char *)&num_p_pc_from_gr_to_pc_p2, sizeof(int));
	out_param_buf.write((char *)&num_p_mf_from_mf_to_nc, sizeof(int));
	out_param_buf.write((char *)&num_p_nc_from_mf_to_nc, sizeof(int));
	out_param_buf.write((char *)&num_p_nc_from_nc_to_io, sizeof(int));
	out_param_buf.write((char *)&num_p_io_from_nc_to_io, sizeof(int));
	out_param_buf.write((char *)&num_p_io_from_io_to_pc, sizeof(int));
	out_param_buf.write((char *)&num_p_io_in_io_to_io, sizeof(int));
	out_param_buf.write((char *)&num_p_io_out_io_to_io, sizeof(int));

	out_param_buf.write((char *)&msPerTimeStep, sizeof(float));
	out_param_buf.write((char *)&numPopHistBinsPC, sizeof(float));

	//out_param_buf.write((char *)&eLeakGR, sizeof(float)); 
	//out_param_buf.write((char *)&threshRestGR, sizeof(float));

	out_param_buf.write((char *)&eLeakSC, sizeof(float));
	out_param_buf.write((char *)&threshRestSC, sizeof(float));
	out_param_buf.write((char *)&eLeakBC, sizeof(float));
	out_param_buf.write((char *)&threshRestBC, sizeof(float));
	out_param_buf.write((char *)&eLeakPC, sizeof(float));
	out_param_buf.write((char *)&threshRestPC, sizeof(float));
	out_param_buf.write((char *)&initSynWofGRtoPC, sizeof(float));
	out_param_buf.write((char *)&eLeakIO, sizeof(float));
	out_param_buf.write((char *)&threshRestIO, sizeof(float));
	out_param_buf.write((char *)&eLeakNC, sizeof(float));
	out_param_buf.write((char *)&threshRestNC, sizeof(float));
}

