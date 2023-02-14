/*
 * mzoneactivitystate.cpp
 *
 *  Created on: Nov 27, 2012
 *      Author: consciousness
 */

#include <iostream>
#include <algorithm> /* std::fill */

#include "logger.h"
#include "file_utility.h"
#include "sfmt.h"
#include "connectivityparams.h"
#include "activityparams.h"
#include "mzoneactivitystate.h"

MZoneActivityState::MZoneActivityState() {}

MZoneActivityState::MZoneActivityState(int randSeed)
{
	LOG_DEBUG("Allocating and initializing mzone activity state...");
	allocateMemory();
	initializeVals(randSeed);
	LOG_DEBUG("Finished allocating and initializing mzone activity state.");
}

MZoneActivityState::MZoneActivityState(std::fstream &infile)
{
	allocateMemory();
	stateRW(true, infile);
}

MZoneActivityState::~MZoneActivityState() {}

void MZoneActivityState::readState(std::fstream &infile)
{
	stateRW(true, infile);
}

void MZoneActivityState::writeState(std::fstream &outfile)
{
	stateRW(false, outfile);
}

void MZoneActivityState::allocateMemory()
{
	// stellate cells
	apSC           = std::make_unique<uint8_t[]>(num_sc);
	apBufSC        = std::make_unique<uint32_t[]>(num_sc);
	gPFSC          = std::make_unique<float[]>(num_sc);
	threshSC       = std::make_unique<float[]>(num_sc);
	vSC            = std::make_unique<float[]>(num_sc);

	// basket cells
	apBC      = std::make_unique<uint8_t[]>(num_bc);
	apBufBC   = std::make_unique<uint32_t[]>(num_bc);
	inputPCBC = std::make_unique<uint32_t[]>(num_bc);
	gPFBC     = std::make_unique<float[]>(num_bc);
	gPCBC     = std::make_unique<float[]>(num_bc);
	vBC       = std::make_unique<float[]>(num_bc);
	threshBC  = std::make_unique<float[]>(num_bc);

	// purkinje cells
	apPC          = std::make_unique<uint8_t[]>(num_pc);
	apBufPC       = std::make_unique<uint32_t[]>(num_pc);
	inputBCPC     = std::make_unique<uint32_t[]>(num_pc);
	inputSCPC     = std::make_unique<uint32_t[]>(num_pc);
	pfSynWeightPC = std::make_unique<float[]>(num_pc * num_p_pc_from_gr_to_pc);
	inputSumPFPC  = std::make_unique<float[]>(num_pc);
	gPFPC         = std::make_unique<float[]>(num_pc);
	gBCPC         = std::make_unique<float[]>(num_pc);
	gSCPC         = std::make_unique<float[]>(num_pc);
	vPC           = std::make_unique<float[]>(num_pc);
	threshPC      = std::make_unique<float[]>(num_pc);
	histPCPopAct = std::make_unique<uint32_t[]>(numPopHistBinsPC);

	// inferior olivary cells
	apIO      = std::make_unique<uint8_t[]>(num_io);
	apBufIO   = std::make_unique<uint8_t[]>(num_io);
	inputNCIO = std::make_unique<uint8_t[]>(num_io * num_p_io_from_nc_to_io);
	gNCIO     = std::make_unique<float[]>(num_io * num_p_io_from_nc_to_io);
	threshIO  = std::make_unique<float[]>(num_io);
	vIO       = std::make_unique<float[]>(num_io);
	vCoupleIO = std::make_unique<float[]>(num_io);
	pfPCPlastTimerIO = std::make_unique<int32_t[]>(num_io);

	// nucleus cells
	apNC          = std::make_unique<uint8_t[]>(num_nc);
	apBufNC       = std::make_unique<uint32_t[]>(num_nc);
	inputPCNC     = std::make_unique<uint8_t[]>(num_nc * num_p_nc_from_pc_to_nc);
	gPCNC         = std::make_unique<float[]>(num_nc * num_p_nc_from_pc_to_nc);
	threshNC      = std::make_unique<float[]>(num_nc);
	vNC           = std::make_unique<float[]>(num_nc);
	synIOPReleaseNC = std::make_unique<float[]>(num_nc);
}

void MZoneActivityState::initializeVals(int randSeed)
{
	// only actively initializing those arrays whose initial values we want
	// differ from the default initilizer value

	// sc
	std::fill(vSC.get(), vSC.get() + num_sc, eLeakSC);
	std::fill(threshSC.get(), threshSC.get() + num_sc, threshRestSC);

	// bc
	std::fill(vBC.get(), vBC.get() + num_bc, eLeakBC);
	std::fill(threshBC.get(), threshBC.get() + num_bc, threshRestBC);

	// pc
	std::fill(vPC.get(), vPC.get() + num_pc, eLeakPC);
	std::fill(threshPC.get(), threshPC.get() + num_pc, threshRestPC);

	std::fill(pfSynWeightPC.get(), pfSynWeightPC.get()
		+ num_pc * num_p_pc_from_gr_to_pc, initSynWofGRtoPC);

	std::fill(histPCPopAct.get(), histPCPopAct.get() + (int)numPopHistBinsPC, 0);

	histPCPopActSum     = 0;
	histPCPopActCurBinN = 0;
	pcPopAct            = 0;

	// IO
	std::fill(vIO.get(), vIO.get() + num_io, eLeakIO);
	std::fill(threshIO.get(), threshIO.get() + num_io, threshRestIO);

	errDrive = 0;
	
	// NC
	std::fill(vNC.get(), vNC.get() + num_nc, eLeakNC);
	std::fill(threshNC.get(), threshNC.get() + num_nc, threshRestNC);

}

void MZoneActivityState::stateRW(bool read, std::fstream &file)
{
	// stellate cells
	rawBytesRW((char *)apSC.get(), num_sc * sizeof(uint8_t), read, file);
	rawBytesRW((char *)apBufSC.get(), num_sc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)gPFSC.get(), num_sc * sizeof(float), read, file);
	rawBytesRW((char *)threshSC.get(), num_sc * sizeof(float), read, file);
	rawBytesRW((char *)vSC.get(), num_sc * sizeof(float), read, file);

	// basket cells
	rawBytesRW((char *)apBC.get(), num_bc * sizeof(uint8_t), read, file);
	rawBytesRW((char *)apBufBC.get(), num_bc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)inputPCBC.get(), num_bc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)gPFBC.get(), num_bc * sizeof(float), read, file);
	rawBytesRW((char *)gPCBC.get(), num_bc * sizeof(float), read, file);
	rawBytesRW((char *)vBC.get(), num_bc * sizeof(float), read, file);
	rawBytesRW((char *)threshBC.get(), num_bc * sizeof(float), read, file);

	// purkinje cells
	rawBytesRW((char *)apPC.get(), num_pc * sizeof(uint8_t), read, file);
	rawBytesRW((char *)apBufPC.get(), num_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)inputBCPC.get(), num_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)inputSCPC.get(), num_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pfSynWeightPC.get(),
		num_pc * num_p_pc_from_gr_to_pc * sizeof(float), read, file);
	rawBytesRW((char *)inputSumPFPC.get(), num_pc * sizeof(float), read, file);
	rawBytesRW((char *)gPFPC.get(), num_pc * sizeof(float), read, file);
	rawBytesRW((char *)gBCPC.get(), num_pc * sizeof(float), read, file);
	rawBytesRW((char *)gSCPC.get(), num_pc * sizeof(float), read, file);
	rawBytesRW((char *)vPC.get(), num_pc * sizeof(float), read, file);
	rawBytesRW((char *)threshPC.get(), num_pc * sizeof(float), read, file);
	rawBytesRW((char *)histPCPopAct.get(), numPopHistBinsPC * sizeof(uint32_t), read, file);

	rawBytesRW((char *)&histPCPopActSum, sizeof(uint32_t), read, file);
	rawBytesRW((char *)&histPCPopActCurBinN, sizeof(uint32_t), read, file);
	rawBytesRW((char *)&pcPopAct, sizeof(uint32_t), read, file);
	
	// inferior olivary cells
	rawBytesRW((char *)apIO.get(), num_io * sizeof(uint8_t), read, file);
	rawBytesRW((char *)apBufIO.get(), num_io * sizeof(uint8_t), read, file);
	rawBytesRW((char *)inputNCIO.get(),
		num_io * num_p_io_from_nc_to_io * sizeof(uint8_t), read, file);
	rawBytesRW((char *)gNCIO.get(),
		num_io * num_p_io_from_nc_to_io * sizeof(float), read, file);
	rawBytesRW((char *)threshIO.get(), num_io * sizeof(float), read, file);
	rawBytesRW((char *)vIO.get(), num_io * sizeof(float), read, file);
	rawBytesRW((char *)vCoupleIO.get(), num_io * sizeof(float), read, file);
	rawBytesRW((char *)pfPCPlastTimerIO.get(), num_io * sizeof(int32_t), read, file);

	rawBytesRW((char *)&errDrive, sizeof(float), read, file);

	// nucleus cells
	rawBytesRW((char *)apNC.get(), num_nc * sizeof(uint8_t), read, file);
	rawBytesRW((char *)apBufNC.get(), num_nc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)inputPCNC.get(),
		num_nc * num_p_nc_from_pc_to_nc * sizeof(uint8_t), read, file);
	rawBytesRW((char *)gPCNC.get(),
		num_nc * num_p_nc_from_pc_to_nc * sizeof(float), read, file);
	rawBytesRW((char *)threshNC.get(), num_nc * sizeof(float), read, file);
	rawBytesRW((char *)vNC.get(), num_nc * sizeof(float), read, file);
	rawBytesRW((char *)synIOPReleaseNC.get(), num_nc * sizeof(float), read, file);
}

