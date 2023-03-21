/*
 * cbmsimcore.h
 *
 *  Created on: Dec 14, 2011
 *      Author: consciousness
 */

#ifndef CBMSIMCORE_H_
#define CBMSIMCORE_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <limits.h>
#include <time.h>
#include <cstdint>

#include "poissonregencells.h"
#include "mzonestate.h"
#include "sfmt.h"
#include "mzone.h"

/* TODO: consider altering this code so that CBMSimCore does not keep local copies
 *       of the state classes. Consider whether transferring data between classes
 *       by using classes as arguments would be just as fast as we have things now.
 *       The advantage would be that we would use less memory and it would simplify the code.
 */

class CBMSimCore
{
public:
	CBMSimCore();
	CBMSimCore(std::fstream &psth_file_buf, MZoneState *state, int gpuIndStart = -1, int numGPUP2 = -1);
	~CBMSimCore();

	void calcActivity(enum plasticity pf_pc_plast);
	//void updateMFInput(const uint8_t *mfIn);
	//void updateTrueMFs(bool *isTrueMF);
	//void updateGRStim(int startGRStim, int numGRStim);
	void updateErrDrive(unsigned int zoneN, float errDriveRelative);

	void writeToState();
	void writeState(std::fstream& outfile);

	MZone** getMZoneList();

protected:
	void initCUDAStreams();
	void initAuxVars();

	void syncCUDA(std::string title);

	PoissonRegenCells *grs = nullptr;
	MZoneState *simState = nullptr;

	uint32_t numZones;

	MZone **zones;

	cudaStream_t **streams;
	int gpuIndStart;
	int numGPUs;

private:
	bool isGRStim    =  false;
	int numGRStim    =  0;
	int startGRStim  =  0;

	uint32_t curTime;

	void construct(std::fstream &psth_file_buf, MZoneState *state, int *mzoneRSeed,
		int gpuIndStart, int numGPUP2);
};

#endif /* CBMSIMCORE_H_ */

