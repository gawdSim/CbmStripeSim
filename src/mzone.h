/*
 * mzone.h
 *
 *  Created on: Jun 13, 2011
 *      Author: consciousness
 */

#ifndef MZONE_H_
#define MZONE_H_

#include <cstdint>
#include "mzoneconnectivitystate.h"
#include "mzoneactivitystate.h"
#include "kernels.h"

class MZone
{
public:
	MZone();
	MZone(MZoneConnectivityState *cs, MZoneActivityState *as, uint32_t **apBufGR,
			uint64_t **apHistGR, int randSeed, int gpuIndStart, int numGPUs);
	~MZone();

	void writeToState();
	void cpyPFPCSynWCUDA();

	void setErrDrive(float errDriveRelative);
	//void updateMFActivities(const uint8_t *actMF);
	//void updateTrueMFs(bool *trueMF);

	void calcPCActivities();
	void runSCActivitiesCUDA(cudaStream_t **sts, int streamN);
	void runBCActivitiesCUDA(cudaStream_t **sts, int streamN);
	void calcIOActivities();
	void calcNCActivities();

	void updatePCOut();
	void updateBCPCOut();
	void updateSCPCOut();
	void updateIOOut();
	void updateNCOut();
	//void updateMFNCOut();
	//void updateMFNCSyn(const uint8_t *histMF, uint32_t t);

	void runUpdatePFPCOutCUDA(cudaStream_t **sts, int streamN);
	void runSumPFPCCUDA(cudaStream_t **sts, int streamN);
	void cpyPFPCSumCUDA(cudaStream_t **sts, int streamN);
	void runPFPCPlastCUDA(cudaStream_t **sts, int streamN, uint32_t t);

	void runSumPFSCCUDA(cudaStream_t **sts, int streamN);

	void runUpdatePFBCOutCUDA(cudaStream_t **sts, int streamN);
	void runUpdatePFSCOutCUDA(cudaStream_t **sts, int streamN);

	void runSumPFBCCUDA(cudaStream_t **sts, int streamN);
	void cpyPFBCSumGPUtoHostCUDA(cudaStream_t **sts, int streamN);

	void setGRPCPlastSteps(float ltdStep, float ltpStep);
	void resetGRPCPlastSteps();

	const uint8_t* exportAPNC();
	const uint8_t* exportAPSC();
	const uint8_t* exportAPBC();
	const uint8_t* exportAPPC();
	const uint8_t* exportAPIO();

	const float* exportVmBC();
	const float* exportVmPC();
	const float* exportVmNC();
	const float* exportVmIO();
	const float* exportgBCPC();
	const float* exportgPFPC();
	const float* exportPFPCWeights();
	//const float* exportMFDCNWeights();

	void load_pfpc_weights_from_file(std::fstream &in_file_buf);
	//void load_mfdcn_weights_from_file(std::fstream &in_file_buf);

	const uint32_t* exportAPBufBC();
	const uint32_t* exportAPBufPC();
	const uint8_t* exportAPBufIO();
	const uint32_t* exportAPBufNC();

private:
	MZoneConnectivityState *cs;
	MZoneActivityState *as;

	CRandomSFMT0 *randGen;

	int gpuIndStart;
	int numGPUs;
	int numGRPerGPU;
	int numSCPerGPU;
	int numBCPerGPU;

	unsigned int updatePFPCNumGRPerB;
	unsigned int updatePFPCNumBlocks;

	unsigned int updatePFPCSynWNumGRPerB;
	unsigned int updatePFPCSynWNumBlocks;

	unsigned int updatePFBCSCNumGRPerB;
	unsigned int updatePFBCSCNumBlocks;

	unsigned int calcSCActNumSCPerB;
	unsigned int calcSCActNumBlocks;

	unsigned int calcBCActNumBCPerB;
	unsigned int calcBCActNumBlocks;

	/* ======== not used ====== */
	unsigned int updateGRBCOutNumGRPerR;
	unsigned int updateGRBCOutNumGRRows;

	unsigned int sumGRBCOutNumBCPerB;
	unsigned int sumGRBCOutNumBlocks;
	/* ======== not used ====== */

	//mossy fiber variables
	//const uint8_t *apMFInput;
	//const uint8_t *histMFInput;
	//bool *isTrueMF;

	//stellate cell variables
	//host variables
	uint8_t **apSCGPU;
	uint32_t **apBufSCGPU;
	float **gPFSCGPU;
	float **threshSCGPU;
	float **vSCGPU;
	uint32_t **inputPFSCGPU;
	size_t *inputPFSCGPUP;
	uint32_t **inputSumPFSCGPU;
	//end gpu related variables
	//end stellate cell variables

	//basket cell variables
	//host variables
	uint32_t *inputSumPFBCH;

	uint8_t **apBCGPU;
	uint32_t **apBufBCGPU;
	float **gPFBCGPU;
	float **gPCBCGPU;
	float **threshBCGPU;
	float **vBCGPU;

	uint32_t **inputPFBCGPU;
	size_t *inputPFBCGPUP;
	uint32_t **inputSumPFBCGPU;

	uint32_t **inputPCBCGPU;
	size_t *inputPCBCGPUP;
	uint32_t **inputSumPCBCGPU;
	//end gpu related variables
	//end basket cell variables

	//purkinje cell variables
	float **pfSynWeightPCGPU;
	float *pfSynWeightPCLinear;
	float **inputPFPCGPU;
	size_t *inputPFPCGPUPitch;
	float **inputSumPFPCMZGPU;
	float *inputSumPFPCMZH;

	uint32_t **apBufGRGPU;
	uint64_t **histGRGPU;
	uint32_t **delayMaskGRGPU;

	//IO cell variables
	float *pfPCPlastStepIO;
	float tempGRPCLTDStep;
	float tempGRPCLTPStep;

	void initCUDA();
	void initBCCUDA();
	void initSCCUDA();
	void testReduction();
};

#endif /* MZONE_H_ */

