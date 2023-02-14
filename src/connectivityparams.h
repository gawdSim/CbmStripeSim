/*
 * connectivityparams.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: varicella
 */
#ifndef CONNECTIVITYPARAMS_H_
#define CONNECTIVITYPARAMS_H_

#include <iostream>
#include <fstream>
#include <string>
#include "file_parse.h"
#include <cstdint>

extern bool con_params_populated;

extern int gr_x; 
extern int gr_y; 
extern int num_gr; 
extern int num_bc; 
extern int num_sc; 
extern int num_pc; 
extern int num_nc; 
extern int num_io; 
extern int gr_pf_vel_in_gr_x_per_t_step; 
extern int gr_af_delay_in_t_step; 
extern int num_p_bc_from_bc_to_pc; 
extern int num_p_pc_from_bc_to_pc; 
extern int num_p_bc_from_gr_to_bc; 
extern int num_p_bc_from_gr_to_bc_p2; 
extern int num_p_pc_from_pc_to_bc; 
extern int num_p_bc_from_pc_to_bc; 
extern int num_p_sc_from_sc_to_pc; 
extern int num_p_pc_from_sc_to_pc; 
extern int num_p_sc_from_gr_to_sc; 
extern int num_p_sc_from_gr_to_sc_p2; 
extern int num_p_pc_from_pc_to_nc; 
extern int num_p_nc_from_pc_to_nc; 
extern int num_p_pc_from_gr_to_pc; 
extern int num_p_pc_from_gr_to_pc_p2; 
extern int num_p_mf_from_mf_to_nc; 
extern int num_p_nc_from_mf_to_nc; 
extern int num_p_nc_from_nc_to_io; 
extern int num_p_io_from_nc_to_io; 
extern int num_p_io_from_io_to_pc; 
extern int num_p_io_in_io_to_io; 
extern int num_p_io_out_io_to_io; 

extern float msPerTimeStep;
extern float numPopHistBinsPC; 

//extern float eLeakGR; 
//extern float threshRestGR;

extern float eLeakSC;
extern float threshRestSC;
extern float eLeakBC;
extern float threshRestBC;
extern float eLeakPC;
extern float threshRestPC;
extern float initSynWofGRtoPC;
extern float eLeakIO;
extern float threshRestIO;
extern float eLeakNC;
extern float threshRestNC;

void populate_con_params(parsed_build_file &p_file);

void read_con_params(std::fstream &in_param_buf);

void write_con_params(std::fstream &out_param_buf);

#endif /* CONNECTIVITYPARAMS_H_ */


