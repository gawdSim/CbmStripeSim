/*
 * kernels.h
 *
 *  Created on: Jun 6, 2011
 *      Author: consciousness
 */

#ifndef KERNELS_H_
#define KERNELS_H_

#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

#include <cstdint>

void callTestKernel(cudaStream_t &st, float *a, float *b, float *c);

void callGRActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		float *vGPU, float *gKCaGPU, float *gLeakGRPGU, float *gNMDAGRGPU, float*gNMDAIncGRGPU,
		float *threshGPU, uint32_t *apBufGPU, uint8_t *apOutGRGPU, uint32_t *apGRGPU,
		int *apMFtoGRGPU, float *gESumGPU, float *gISumGPU, float eLeak, float eGOIn,
		float gAMPAInc, float threshBase, float threshMax, float threshDecay);

void callPCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numPCPerBlock,
	float *vPC, float *gPFPC, float *gSCPC, float *gBCPC, float *threshPC, uint8_t *apPC, uint32_t *apBufPC,
	float *gInputSumPFPC, float *gInputSumSCPC, float *gInputSumBCPC, float eLeakPC, float gLeakPC, float gIncPFPC,
	float gIncSCPC, float gIncBCPC, float eSCtoPC, float eBCtoPC, float gDecPFPC, float gDecSCPC, float gDecBCPC,
	float threshRestPC, float threshMaxPC, float threshDecayPC);

void callSCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numSCPerBlock, 
	float *vSC, float *gPFSC, float *threshSC, uint8_t *apSC, uint32_t *apBufSC, float *gInputSumPFSC,
	float eLeakSC, float gLeakSC, float gIncPFSC, float gDecPFSC, float threshRestSC, float threshMaxSC,
	float threshDecaySC);

void callBCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numBCPerBlock, 
	float *vBC, float *gPFBC, float *gPCBC, float *threshBC, uint8_t *apBC, uint32_t *apBufBC,
	float *gInputSumPFBC, float *gInputSumPCBC, float eLeakBC, float gLeakBC, float gIncPFBC,
	float gIncPCBC, float ePCtoBC, float gDecPFBC, float gDecPCBC, float threshRestBC,
	float threshMaxBC, float threshDecayBC);

void callNCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numNCPerBlock,
	float *vNC, float *threshNC, uint8_t *apNC, uint32_t *apBufNC, float *gInputSumPCNC,
	float *gInputSumMFNCNMDA, float *gInputSumMFNCAMPA, float eLeakNC, float gLeakNC,
	float ePCtoNC, float threshRestNC, float threshMaxNC, float threshDecayNC);

template<typename Type, bool inMultiP, bool outMultiP>
void callSumKernel(cudaStream_t &st, Type *inGPU, size_t inGPUP, Type *outSumGPU, size_t outSumGPUP,
		unsigned int nOutCells, unsigned int nOutCols, unsigned int rowLength);

template<typename Type>
void callBroadcastKernel(cudaStream_t &st, Type *broadCastVal, Type *outArray,
		unsigned int nBlocks, unsigned int rowLength);

void callUpdateGROutBCKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock, unsigned int numBC,
		uint32_t *apBufGPU, uint32_t *grInBCGPU, uint32_t grInBCGPUPitch,
		uint32_t *delayMasksGPU, uint32_t delayMasksGPUPitch,
		uint32_t *conGRtoBCGPU, size_t conGRtoBCGPUPitch,
		int32_t *numBCPerGRGPU);

void callSumGRBCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGOPerBlock,
		unsigned int numGROutRows, uint32_t *grInBCGPU,  size_t grInBCGPUPitch, uint32_t *grInBCSGPU);

void callUpdatePFBCSCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint32_t *delayMaskGPU,
		uint32_t *inPFBCGPU, size_t inPFBCGPUPitch, unsigned int numPFInPerBCP2,
		uint32_t *inPFSCGPU, size_t inPFSCGPUPitch, unsigned int numPFInPerSCP2);

void callUpdatePFPCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint32_t *delayMaskGPU,
		float *pfPCSynWGPU, float *inPFPCGPU, size_t inPFPCGPUPitch, unsigned int numPFInPerPCP2);

void callUpdateGRHistKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint64_t *historyGPU, uint32_t apBufGRHistMask);

void callUpdatePFPCPlasticityIOKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		float *synWeightGPU, uint64_t *historyGPU, unsigned int pastBinNToCheck,
		int offSet, float pfPCPlastStep);

#endif /* KERNELS_H_ */

