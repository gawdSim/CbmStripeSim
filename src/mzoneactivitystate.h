/*
 * mzoneactivitystate.h
 *
 *  Created on: Nov 26, 2012
 *      Author: consciousness
 */

#ifndef MZONEACTIVITYSTATE_H_
#define MZONEACTIVITYSTATE_H_

#include <fstream>
#include <memory> /* unique_ptr, make_unique */
#include <cstdint>

class MZoneActivityState
{
public:
	MZoneActivityState();
	MZoneActivityState(int randSeed);
	MZoneActivityState(std::fstream &infile);

	~MZoneActivityState();
	
	void readState(std::fstream &infile);
	void writeState(std::fstream &outfile);

	//stellate cells
	std::unique_ptr<uint8_t[]> apSC{nullptr};
	std::unique_ptr<uint32_t[]> apBufSC{nullptr};
	std::unique_ptr<float[]> gPFSC{nullptr};
	std::unique_ptr<float[]> threshSC{nullptr};
	std::unique_ptr<float[]> vSC{nullptr};

	//basket cells
	std::unique_ptr<uint8_t[]> apBC{nullptr};
	std::unique_ptr<uint32_t[]> apBufBC{nullptr};
	std::unique_ptr<uint32_t[]> inputPCBC{nullptr};
	std::unique_ptr<float[]> gPFBC{nullptr};
	std::unique_ptr<float[]> gPCBC{nullptr};
	std::unique_ptr<float[]> vBC{nullptr};
	std::unique_ptr<float[]> threshBC{nullptr};

	//purkinje cells
	std::unique_ptr<uint8_t[]> apPC{nullptr};
	std::unique_ptr<uint32_t[]> apBufPC{nullptr};
	std::unique_ptr<uint32_t[]> inputBCPC{nullptr};
	std::unique_ptr<uint32_t[]> inputSCPC{nullptr};
	std::unique_ptr<float[]> pfSynWeightPC{nullptr};
	std::unique_ptr<float[]> inputSumPFPC{nullptr};
	std::unique_ptr<float[]> gPFPC{nullptr};
	std::unique_ptr<float[]> gBCPC{nullptr};
	std::unique_ptr<float[]> gSCPC{nullptr};
	std::unique_ptr<float[]> vPC{nullptr};
	std::unique_ptr<float[]> threshPC{nullptr};
	std::unique_ptr<uint32_t[]> histPCPopAct{nullptr};

	uint32_t histPCPopActSum;
	uint32_t histPCPopActCurBinN;
	uint32_t pcPopAct;

	//inferior olivary cells
	std::unique_ptr<uint8_t[]> apIO{nullptr};
	std::unique_ptr<uint8_t[]> apBufIO{nullptr}; /* should this be a byte array? */
	std::unique_ptr<uint8_t[]> inputNCIO{nullptr};
	std::unique_ptr<float[]> gNCIO{nullptr};
	std::unique_ptr<float[]> threshIO{nullptr};
	std::unique_ptr<float[]> vIO{nullptr};
	std::unique_ptr<float[]> vCoupleIO{nullptr};
	std::unique_ptr<int32_t[]> pfPCPlastTimerIO{nullptr};

	float errDrive;

	//nucleus cells
	std::unique_ptr<uint8_t[]> apNC{nullptr};
	std::unique_ptr<uint32_t[]> apBufNC{nullptr};
	std::unique_ptr<uint8_t[]> inputPCNC{nullptr};
	std::unique_ptr<float[]> gPCNC{nullptr};
	std::unique_ptr<float[]> threshNC{nullptr};
	std::unique_ptr<float[]> vNC{nullptr};
	std::unique_ptr<float[]> synIOPReleaseNC{nullptr};

private:
	void allocateMemory();
	void initializeVals(int randSeed);
	void stateRW(bool read, std::fstream &file);
};

#endif /* MZONEACTIVITYSTATE_H_ */

