/*
 * mzone.cpp
 *
 *   Created on: Jun 14, 2011
 *   Author: consciousness
 *
 */
#include <math.h>
#include <iostream>
#include <fstream>

#include "logger.h"
#include "connectivityparams.h"
#include "activityparams.h"
#include "dynamic2darray.h"
#include "sfmt.h"
#include "file_utility.h"
#include "mzone.h"

MZone::MZone() {}

MZone::MZone(MZoneConnectivityState *cs, MZoneActivityState *as, uint32_t **apBufGR,
			uint64_t **apHistGR, int randSeed, int gpuIndStart, int numGPUs)
{
	randGen = new CRandomSFMT0(randSeed);

	// shallow copies. caller owns the data.
	this->cs = cs; 
	this->as = as; 

	// NOTE if we turn these guys into unique ptrs, we'll have to refactor
	// consider ownership: who should own these guys? maybe they should be global to both
	// innet and mzone (so within cbmsimcore) and fed in as const args to the respective
	// functions that call update kernels (06/16/2022)

	this->apBufGRGPU = apBufGR;
	this->histGRGPU  = apHistGR;

	delayMaskGRGPU = new uint32_t*[numGPUs];

	pfSynWeightPCLinear = new float[num_gr];
	pfPCPlastStepIO     = new float[num_io];

	tempGRPCLTDStep = synLTDStepSizeGRtoPC;
	tempGRPCLTPStep = synLTPStepSizeGRtoPC;

	this->numGPUs     = numGPUs;
	this->gpuIndStart = gpuIndStart;

	LOG_DEBUG("Initializing CUDA...");
	initCUDA();
}

MZone::~MZone()
{
	LOG_DEBUG("Deleting mzone gpu arrays...");

	delete randGen;

	delete[] pfSynWeightPCLinear;
	delete[] pfPCPlastStepIO;

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaFreeHost(inputSumPFPCMZH[i]);

		//free cuda device memory
		cudaFree(delayMaskGRGPU[i]);
		cudaFree(gr_pc_con_in_d[i]);
		cudaFree(pfSynWeightPCGPU[i]);
		cudaFree(inputPFPCGPU[i]);
		cudaFree(inputSumPFPCMZGPU[i]);
		cudaDeviceSynchronize();
	}

	delete[] delayMaskGRGPU;
	delete[] inputSumPFPCMZH;
	delete[] inputPFPCSumH;

	delete[] gr_pc_con_in_d;
	delete[] pfSynWeightPCGPU;
	delete[] inputPFPCGPU;
	delete[] inputPFPCGPUPitch;
	delete[] inputSumPFPCMZGPU;

	// sc
	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaFreeHost(inputSumPFSCMZH[i]);

		cudaSetDevice(i+gpuIndStart);
		//cudaFree(apSCGPU[i]);
		//cudaFree(apBufSCGPU[i]);
		//cudaFree(gPFSCGPU[i]);
		//cudaFree(threshSCGPU[i]);
		//cudaFree(vSCGPU[i]);
		cudaFree(gr_sc_con_in_d[i]);
		cudaFree(inputPFSCGPU[i]);
		cudaFree(inputSumPFSCGPU[i]);
		cudaDeviceSynchronize();
	}

	delete[] inputSumPFSCMZH;
	delete[] inputSumPFSCH;

	//delete[] apSCGPU;
	//delete[] apBufSCGPU;
	//delete[] gPFSCGPU;
	//delete[] threshSCGPU;
	//delete[] vSCGPU;
	delete[] gr_sc_con_in_d;
	delete[] inputPFSCGPU;
	delete[] inputPFSCGPUPitch;
	delete[] inputSumPFSCGPU;

	//bc
	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaFreeHost(inputSumPFBCMZH[i]);

		//cudaFree(apBCGPU[i]);
		//cudaFree(apBufBCGPU[i]);
		//cudaFree(gPFBCGPU[i]);
		//cudaFree(gPCBCGPU[i]);
		//cudaFree(threshBCGPU[i]);
		//cudaFree(vBCGPU[i]);
		cudaFree(gr_bc_con_in_d[i]);
		cudaFree(inputPFBCGPU[i]);
		cudaFree(inputSumPFBCGPU[i]);
		cudaDeviceSynchronize();
	}

	delete[] inputSumPFBCMZH;
	delete[] inputSumPFBCH;

	//delete[] apBCGPU;
	//delete[] apBufBCGPU;
	//delete[] gPFBCGPU;
	//delete[] gPCBCGPU;
	//delete[] threshBCGPU;
	//delete[] vBCGPU;
	delete[] gr_bc_con_in_d;
	delete[] inputPFBCGPU;
	delete[] inputPFBCGPUPitch;
	delete[] inputSumPFBCGPU;

	LOG_DEBUG("Finished deleting mzone gpu arrays.");
}

void MZone::initCUDA()
{
	int maxNumGPUs;
	cudaGetDeviceCount(&maxNumGPUs);

	numGRPerGPU = num_gr / numGPUs;
	numSCPerGPU = num_sc / numGPUs;
	numBCPerGPU = num_bc / numGPUs;

	updatePFPCNumGRPerB = 512 * (num_pc > 512) + num_pc * (num_pc <= 512);
	updatePFPCNumBlocks = numGRPerGPU / updatePFPCNumGRPerB;

	sumPFPCOutNumPCPerB = 1024 * (num_pc > 1024) + num_pc * (num_pc <= 1024);
	sumPFPCOutNumBlocks = num_pc / sumPFPCOutNumPCPerB;

	updatePFBCNumGRPerB = 512 * (num_bc > 512) + num_bc * (num_bc <= 512);
	updatePFBCNumBlocks = numGRPerGPU / updatePFBCNumGRPerB;

	sumPFBCOutNumBCPerB = 1024 * (num_bc > 1024) + num_bc * (num_bc <= 1024);
	sumPFBCOutNumBlocks = num_bc / sumPFBCOutNumBCPerB;

	updatePFSCNumGRPerB = 512 * (num_sc > 512) + num_sc * (num_sc <= 512);
	updatePFSCNumBlocks = numGRPerGPU / updatePFSCNumGRPerB;

	sumPFSCOutNumSCPerB = 1024 * (num_sc > 1024) + num_sc * (num_sc <= 1024);
	sumPFSCOutNumBlocks = num_sc / sumPFSCOutNumSCPerB;

	updatePFPCSynWNumGRPerB = 512 * (num_p_pc_from_gr_to_pc > 512) +
			num_p_pc_from_gr_to_pc * (num_p_pc_from_gr_to_pc <= 512);
	updatePFPCSynWNumBlocks = num_p_pc_from_gr_to_pc / updatePFPCSynWNumGRPerB;

	updatePFBCSCNumGRPerB = 512;
	updatePFBCSCNumBlocks = numGRPerGPU / updatePFBCSCNumGRPerB;

	calcSCActNumSCPerB = 512;
	calcSCActNumBlocks = numSCPerGPU / calcSCActNumSCPerB;

	calcBCActNumBCPerB = 512;
	calcBCActNumBlocks = numBCPerGPU / calcBCActNumBCPerB;

	/* ======== not used ====== */
	updateGRBCOutNumGRPerR=512*(num_bc>512)+num_bc*(num_bc<=512);
	updateGRBCOutNumGRRows=numGRPerGPU/updateGRBCOutNumGRPerR;

	sumGRBCOutNumBCPerB=1024*(num_bc>1024)+num_bc*(num_bc<=1024);
	sumGRBCOutNumBlocks=num_bc/sumGRBCOutNumBCPerB;
	/* ======== not used ====== */

	for (int i = 0; i < num_pc; i++)
	{
		for (int j = 0; j < num_p_pc_from_gr_to_pc; j++)
		{
			// TODO: get rid of pfSynWeightLinear and use our linearized version directly
			pfSynWeightPCLinear[i * num_p_pc_from_gr_to_pc + j] = as->pfSynWeightPC[i * num_p_pc_from_gr_to_pc + j];
		}
	}
	inputSumPFPCMZH = new float*[numGPUs];
	inputPFPCSumH   = new float[num_pc];

	gr_pc_con_in_d = new uint32_t*[numGPUs];
	pfSynWeightPCGPU = new float*[numGPUs];
	inputPFPCGPU = new float*[numGPUs];
	inputPFPCGPUPitch = new size_t[numGPUs];
	inputSumPFPCMZGPU = new float*[numGPUs];

	for (int i = 0; i < numGPUs; i++)
	{
		int cpyStartInd = i * numGRPerGPU;
		int cpySize     = numGRPerGPU;

		cudaSetDevice(i + gpuIndStart);

		cudaMallocHost((void **)&inputSumPFPCMZH[i], num_pc * sizeof(float));
		cudaMemset(inputSumPFPCMZH[i], 0, num_pc * sizeof(float));

		// conduction delay variables
		cudaMalloc((void **)&delayMaskGRGPU[i], numGRPerGPU * sizeof(uint32_t));
		cudaMemcpy(delayMaskGRGPU[i], &(cs->pGRDelayMaskfromGRtoBSP[cpyStartInd]),
			cpySize * sizeof(uint32_t), cudaMemcpyHostToDevice);

		cudaMalloc((void **)&gr_pc_con_in_d[i], numGRPerGPU * sizeof(uint32_t));
		cudaMemcpy(gr_pc_con_in_d[i], &(cs->pGRfromGRtoPC[cpyStartInd]), cpySize * sizeof(uint32_t), cudaMemcpyHostToDevice);

		//allocate device cuda memory
		cudaMalloc((void **)&pfSynWeightPCGPU[i], numGRPerGPU * sizeof(float));
		cudaMallocPitch((void **)&inputPFPCGPU[i], (size_t *)&inputPFPCGPUPitch[i],
				num_pc * sizeof(float), updatePFPCNumBlocks);
		cudaMalloc((void **)&inputSumPFPCMZGPU[i], num_pc * sizeof(float));

		cudaDeviceSynchronize();
		//initialize device cuda memory
		cudaMemcpy(pfSynWeightPCGPU[i], &pfSynWeightPCLinear[cpyStartInd],
				numGRPerGPU*sizeof(float), cudaMemcpyHostToDevice);

		for (int j = 0; j < updatePFPCNumBlocks; j++)
		{
			cudaMemset(((char *)inputPFPCGPU[i] + j * inputPFPCGPUPitch[i]),
					0, num_pc * sizeof(float));
		}
		cudaMemset(inputSumPFPCMZGPU[i], 0, num_pc * sizeof(float));

		cudaDeviceSynchronize();
	}
	initBCCUDA();
	LOG_DEBUG("Initialized BC CUDA");
	LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));

	initSCCUDA();
	LOG_DEBUG("Initialized SC CUDA");
	LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));
	
	//testReduction();
	//LOG_DEBUG("Finished Test.");
}

void MZone::initBCCUDA()
{
	//apBCGPU = new uint8_t*[numGPUs];
	//apBufBCGPU = new uint32_t*[numGPUs];
	//gPFBCGPU = new float*[numGPUs];
	//gPCBCGPU = new float*[numGPUs];
	//threshBCGPU = new float*[numGPUs];
	//vBCGPU = new float*[numGPUs];

	gr_bc_con_in_d = new uint32_t*[numGPUs];
	
	inputPFBCGPU    = new uint32_t*[numGPUs];
	inputPFBCGPUPitch   = new size_t[numGPUs];
	inputSumPFBCGPU = new uint32_t*[numGPUs];

	//allocate host memory
	LOG_DEBUG("Allocating BC cuda variables...");

	inputSumPFBCMZH = new uint32_t*[numGPUs];
	inputSumPFBCH   = new uint32_t[num_bc];

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		int cpyStartInd = i * numGRPerGPU;
		int cpySize     = numGRPerGPU;

		cudaMallocHost((void **)&inputSumPFBCMZH[i], num_bc * sizeof(uint32_t));
		cudaMemset(inputSumPFBCMZH[i], 0, num_bc * sizeof(uint32_t));

		//cudaMalloc((void **)&apBCGPU[i], numBCPerGPU * sizeof(uint8_t));
		//cudaMalloc((void **)&apBufBCGPU[i], numBCPerGPU * sizeof(uint32_t));
		//cudaMalloc((void **)&gPFBCGPU[i], numBCPerGPU * sizeof(float));
		//cudaMalloc((void **)&gPCBCGPU[i], numBCPerGPU * sizeof(float));
		//cudaMalloc((void **)&threshBCGPU[i], numBCPerGPU * sizeof(float));
		//cudaMalloc((void **)&vBCGPU[i], numBCPerGPU * sizeof(float));

		cudaMalloc((void **)&gr_bc_con_in_d[i], numGRPerGPU * sizeof(uint32_t));
		cudaMemcpy(gr_bc_con_in_d[i], &(cs->pGRfromGRtoBC[cpyStartInd]), cpySize * sizeof(uint32_t), cudaMemcpyHostToDevice);

		cudaMallocPitch((void **)&inputPFBCGPU[i], (size_t *)&inputPFBCGPUPitch[i],
			num_bc * sizeof(uint32_t), updatePFBCNumBlocks);
		cudaMalloc((void **)&inputSumPFBCGPU[i], num_bc * sizeof(uint32_t));
		cudaDeviceSynchronize();
	}		
	LOG_DEBUG("Finished BC variable cuda allocation");
	LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));

	// initialize BC vars
	LOG_DEBUG("Initializing BC cuda variables...");
	for (int i = 0; i < numGPUs; i++)
	{
		int cpyStartInd = numBCPerGPU * i;
		int cpySize     = numBCPerGPU;
		cudaSetDevice(i + gpuIndStart);
		//cudaMemset(apBCGPU[i], 0, numBCPerGPU * sizeof(uint8_t));
		//cudaMemset(apBufBCGPU[i], 0, numBCPerGPU * sizeof(uint32_t));
		//cudaMemset(gPFBCGPU[i], 0, numBCPerGPU * sizeof(float));
		//cudaMemset(gPCBCGPU[i], 0, numBCPerGPU * sizeof(float));
		//cudaMemcpy(threshBCGPU[i], &(as->threshBC[cpyStartInd]), cpySize * sizeof(float), cudaMemcpyHostToDevice);
		//cudaMemcpy(vBCGPU[i], &(as->vBC[cpyStartInd]), cpySize * sizeof(float), cudaMemcpyHostToDevice);

		for (int j = 0; j < updatePFBCNumBlocks; j++)
		{
			cudaMemset(((char *)inputPFBCGPU[i] + j * inputPFBCGPUPitch[i]),
				0, num_bc * sizeof(uint32_t));
		}
		cudaMemset(inputSumPFBCGPU[i], 0, num_bc * sizeof(uint32_t));
		cudaDeviceSynchronize();
	}
	LOG_DEBUG("Finished initializing BC cuda variables...");
}

void MZone::initSCCUDA()
{
	//apSCGPU = new uint8_t*[numGPUs];
	//apBufSCGPU = new uint32_t*[numGPUs];
	//gPFSCGPU = new float*[numGPUs];
	//threshSCGPU = new float*[numGPUs];
	//vSCGPU = new float*[numGPUs];

	gr_sc_con_in_d = new uint32_t*[numGPUs];

	inputPFSCGPU    = new uint32_t*[numGPUs];
	inputPFSCGPUPitch   = new size_t[numGPUs];
	inputSumPFSCGPU = new uint32_t*[numGPUs];

	//allocate host memory
	LOG_DEBUG("Allocating SC cuda variables...");

	inputSumPFSCMZH = new uint32_t*[numGPUs];
	inputSumPFSCH   = new uint32_t[num_sc];

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		int cpyStartInd = i * numGRPerGPU;
		int cpySize     = numGRPerGPU;

		cudaMallocHost((void **)&inputSumPFSCMZH[i], num_sc * sizeof(uint32_t));
		cudaMemset(inputSumPFSCMZH[i], 0, num_sc * sizeof(uint32_t));

		//cudaMalloc((void **)&apSCGPU[i], numSCPerGPU * sizeof(uint8_t));
		//cudaMalloc((void **)&apBufSCGPU[i], numSCPerGPU * sizeof(uint32_t));
		//cudaMalloc((void **)&gPFSCGPU[i], numSCPerGPU * sizeof(float));
		//cudaMalloc((void **)&threshSCGPU[i], numSCPerGPU * sizeof(float));
		//cudaMalloc((void **)&vSCGPU[i], numSCPerGPU * sizeof(float));

		cudaMalloc((void **)&gr_sc_con_in_d[i], numGRPerGPU * sizeof(uint32_t));
		cudaMemcpy(gr_sc_con_in_d[i], &(cs->pGRfromGRtoSC[cpyStartInd]), cpySize * sizeof(uint32_t), cudaMemcpyHostToDevice);

		cudaMallocPitch((void **)&inputPFSCGPU[i], (size_t *)&inputPFSCGPUPitch[i],
				num_sc * sizeof(uint32_t), updatePFSCNumBlocks);
		cudaMalloc((void **)&inputSumPFSCGPU[i], num_sc * sizeof(uint32_t));
		cudaDeviceSynchronize();
	}
	LOG_DEBUG("Finished SC variable cuda allocation");
	LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));

	// initialize SC vars
	LOG_DEBUG("Initializing SC cuda variables...");
	for (int i = 0; i < numGPUs; i++)
	{
		int cpyStartInd = numSCPerGPU * i;
		int cpySize     = numSCPerGPU;
		cudaSetDevice(i + gpuIndStart);
		//cudaMemset(apSCGPU[i], 0, numSCPerGPU * sizeof(uint8_t));
		//cudaMemset(apBufSCGPU[i], 0, numSCPerGPU * sizeof(uint32_t));
		//cudaMemset(gPFSCGPU[i], 0, numSCPerGPU * sizeof(float));
		//cudaMemcpy(threshSCGPU[i], &(as->threshSC[cpyStartInd]), cpySize * sizeof(float), cudaMemcpyHostToDevice);
		//cudaMemcpy(vSCGPU[i], &(as->vSC[cpyStartInd]), cpySize * sizeof(float), cudaMemcpyHostToDevice);
		for(int j =0; j < updatePFSCNumBlocks; j++)
		{
			cudaMemset(((char *)inputPFSCGPU[i] + j * inputPFSCGPUPitch[i]),
				0, num_sc * sizeof(uint32_t));
		}
		cudaMemset(inputSumPFSCGPU[i], 0, num_sc * sizeof(uint32_t));
		cudaDeviceSynchronize();
	}
	LOG_DEBUG("Finished initializing SC cuda variables...");
}

void MZone::writeToState()
{
	// TODO: write everything to state...only doing weights and pfpc input sums :/
	//cpyPFPCSynWCUDA();

	//for (int i = 0; i < num_pc; i++)
	//{
	//	as->inputSumPFPC[i] = inputSumPFPCMZH[i];
	//}
}

void MZone::cpyPFPCSynWCUDA()
{
	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaMemcpy((void *)&pfSynWeightPCLinear[i*numGRPerGPU], pfSynWeightPCGPU[i],
			numGRPerGPU * sizeof(float), cudaMemcpyDeviceToHost);
	}

	for (int i = 0; i < num_pc; i++)
	{
		for (int j = 0; j < num_p_pc_from_gr_to_pc; j++)
		{
			as->pfSynWeightPC[i * num_p_pc_from_gr_to_pc + j] = pfSynWeightPCLinear[i * num_p_pc_from_gr_to_pc + j];
		}
	}
}

void MZone::setErrDrive(float errDriveRelative)
{
	as->errDrive = errDriveRelative * maxExtIncVIO;
}

void MZone::calcPCActivities()
{
#pragma omp parallel \
default (none), shared (num_pc, gIncGRtoPC, gDecGRtoPC, gIncBCtoPC, gDecBCtoPC, \
	  gIncSCtoPC, gDecSCtoPC, gLeakPC, eLeakPC, eBCtoPC, eSCtoPC, \
	  threshDecPC, threshRestPC, threshMaxPC)
	{
#pragma omp for
		for (int i = 0; i < num_pc; i++)
		{
			inputPFPCSumH[i] = 0;
			for (int j = 0; j < numGPUs; j++)
			{
				inputPFPCSumH[i] += inputSumPFPCMZH[j][i];
			}

			as->gPFPC[i] += inputPFPCSumH[i] * gIncGRtoPC;
			as->gPFPC[i] *= gDecGRtoPC;
			as->gBCPC[i] += as->inputBCPC[i] * gIncBCtoPC;
			as->gBCPC[i] *= gDecBCtoPC;
			as->gSCPC[i] += as->inputSCPC[i] * gIncSCtoPC;
			as->gSCPC[i] *= gDecSCtoPC;

			as->vPC[i] += (gLeakPC * (eLeakPC - as->vPC[i]))
						- (as->gPFPC[i] * as->vPC[i])
						+ (as->gBCPC[i] * (eBCtoPC - as->vPC[i]))
						+ (as->gSCPC[i] * (eSCtoPC - as->vPC[i]));

			as->threshPC[i] += threshDecPC * (threshRestPC - as->threshPC[i]);

			as->apPC[i]    = as->vPC[i] > as->threshPC[i];
			as->apBufPC[i] = (as->apBufPC[i] << 1) | (as->apPC[i] * 0x00000001);

			as->threshPC[i] = as->apPC[i] * threshMaxPC + (!as->apPC[i]) * as->threshPC[i];
			//as->pcPopAct   += as->apPC[i];
		}
	}
}

//void MZone::runSCActivitiesCUDA(cudaStream_t **sts, int streamN)
//{
//	for (int i = 0; i < numGPUs; i++)
//	{
//		cudaSetDevice(i + gpuIndStart);
//		callSCActKernel(sts[i][streamN], calcSCActNumBlocks, calcSCActNumSCPerB, 
//		vSCGPU[i], gPFSCGPU[i], threshSCGPU[i], apSCGPU[i], apBufSCGPU[i], inputSumPFSCGPU[i],
//		eLeakSC, gLeakSC, gIncGRtoSC, gDecGRtoSC, threshRestSC, threshMaxSC, threshDecSC);
//	}
//}

//void MZone::runBCActivitiesCUDA(cudaStream_t **sts, int streamN)
//{
//	for (int i = 0; i < numGPUs; i++)
//	{
//		cudaSetDevice(i + gpuIndStart);
//		callBCActKernel(sts[i][streamN], calcBCActNumBlocks, calcBCActNumBCPerB, 
//			vBCGPU[i], gPFBCGPU[i], gPCBCGPU[i], threshBCGPU[i], apBCGPU[i], apBufBCGPU[i],
//			inputSumPFBCGPU[i], inputSumPCBCGPU[i], eLeakBC, gLeakBC, gIncGRtoBC,
//			gIncPCtoBC, ePCtoBC, gDecGRtoBC, gDecPCtoBC, threshRestBC,
//			threshMaxBC, threshDecBC);
//	}
//}

void MZone::calcBCActivities()
{
#pragma omp parallel \
default (none), shared (num_bc, gIncGRtoBC, gDecGRtoBC, gIncPCtoBC, gDecPCtoBC, \
	  gLeakBC, eLeakBC, ePCtoBC, threshDecBC, threshRestBC, threshMaxBC)
	{
#pragma omp for
		for (int i = 0; i < num_bc; i++)
		{
			inputSumPFBCH[i] = 0;
			for (int j = 0; j < numGPUs; j++)
			{
				inputSumPFBCH[i] += inputSumPFBCMZH[j][i];
			}

			as->gPFBC[i] += inputSumPFBCH[i] * gIncGRtoBC;
			as->gPFBC[i] *= gDecGRtoBC;
			as->gPCBC[i] += as->inputPCBC[i] * gIncPCtoBC;
			as->gPCBC[i] *= gDecPCtoBC;

			as->vBC[i] += (gLeakBC * (eLeakBC - as->vBC[i]))
						- (as->gPFBC[i] * as->vBC[i])
						+ (as->gPCBC[i] * (ePCtoBC - as->vBC[i]));

			as->threshBC[i] += threshDecBC * (threshRestBC - as->threshBC[i]);
			as->apBC[i] = as->vBC[i] > as->threshBC[i];
			as->apBufBC[i] = (as->apBufBC[i] << 1) | (as->apBC[i] * 0x00000001);

			as->threshBC[i] = as->apBC[i] * threshMaxBC + (1-as->apBC[i]) * (as->threshBC[i]);
		}
	}
}

void MZone::calcSCActivities()
{
#pragma omp parallel \
default (none), shared (num_sc, gIncGRtoSC, gDecGRtoSC, gLeakSC, eLeakSC, threshDecSC, threshRestSC, threshMaxSC)
	{
#pragma omp for
		for (int i = 0; i < num_sc; i++)
		{
			inputSumPFSCH[i] = 0;
			for (int j = 0; j < numGPUs; j++)
			{
				inputSumPFSCH[i] += inputSumPFSCMZH[j][i];
			}

			as->gPFSC[i] += inputSumPFSCH[i] * gIncGRtoSC;
			as->gPFSC[i] *= gDecGRtoSC;

			as->vSC[i] += gLeakSC * (eLeakSC - as->vSC[i]) - as->gPFSC[i] * as->vSC[i];
			as->threshSC[i] += threshDecSC * (threshRestSC - as->threshSC[i]);

			as->apSC[i] = as->vSC[i] > as->threshSC[i];
			as->apBufSC[i] = (as->apBufSC[i] << 1) | (as->apSC[i] * 0x00000001);

			as->threshSC[i] = as->apSC[i] * threshMaxSC + (1-as->apSC[i]) * as->threshSC[i];
		}
	}
}

void MZone::calcIOActivities()
{
	clock_t t;
	t = clock();
	srand(t);
	
	float r = static_cast<float> (rand()) / static_cast<float> (RAND_MAX);
	float gNoise = (r - 0.5) * 2.0;

	for (int i = 0; i < num_io; i++)
	{
		float gNCSum;
		gNCSum = 0;

		for (int j = 0; j < num_p_io_from_nc_to_io; j++)
		{
			as->gNCIO[i * num_p_io_from_nc_to_io + j] *= exp(-msPerTimeStep /
				(-gDecTSofNCtoIO * exp(-as->gNCIO[i * num_p_io_from_nc_to_io + j] / gDecTTofNCtoIO)
				 + gDecT0ofNCtoIO));
			as->gNCIO[i * num_p_io_from_nc_to_io + j] += as->inputNCIO[i * num_p_io_from_nc_to_io + j]
			   * gIncNCtoIO * exp(-as->gNCIO[i * num_p_io_from_nc_to_io + j] / gIncTauNCtoIO);
			gNCSum += as->gNCIO[i * num_p_io_from_nc_to_io + j];

			as->inputNCIO[i * num_p_io_from_nc_to_io + j] = 0;
		}

		gNCSum = 1.5 * gNCSum / 3.1;

		as->vIO[i] += gLeakIO * (eLeakIO - as->vIO[i])
					+ gNCSum * (eNCtoIO - as->vIO[i]) + as->vCoupleIO[i]
					+ as->errDrive + gNoise;

		as->apIO[i] = as->vIO[i] > as->threshIO[i];
		as->apBufIO[i] = (as->apBufIO[i] << 1) | (as->apIO[i] * 0x00000001);

		as->threshIO[i] = threshMaxIO * as->apIO[i] +
				(1.0-as->apIO[i]) * (as->threshIO[i] + threshDecIO * (threshRestIO - as->threshIO[i]));
	}
	as->errDrive = 0;
}

//FIXME: getting no spikes in NC. maybe has something to do with this: gPCNCSum *= (msPerTimeStep / (float)num_p_nc_from_pc_to_nc);
void MZone::calcNCActivities()
{
	for (int i = 0; i < num_nc; i++)
	{
		float gPCNCSum = 0;
		int inputPCNCSum = 0;

		for (int j = 0; j < num_p_nc_from_pc_to_nc; j++)
		{
			inputPCNCSum += as->inputPCNC[i * num_p_nc_from_pc_to_nc + j];

			as->gPCNC[i * num_p_nc_from_pc_to_nc + j] = as->gPCNC[i * num_p_nc_from_pc_to_nc + j] * gDecPCtoNC + 
				as->inputPCNC[i * num_p_nc_from_pc_to_nc + j] * gIncAvgPCtoNC
				* (1 - as->gPCNC[i * num_p_nc_from_pc_to_nc + j]);
			gPCNCSum += as->gPCNC[i * num_p_nc_from_pc_to_nc + j];
		}

		gPCNCSum *= (msPerTimeStep / (float)num_p_nc_from_pc_to_nc);
		as->vNC[i] += gLeakNC * (eLeakNC - as->vNC[i]) + gPCNCSum * (ePCtoNC - as->vNC[i]);
		
		as->threshNC[i] += threshDecNC * (threshRestNC - as->threshNC[i]);
		as->apNC[i] = as->vNC[i] > as->threshNC[i];
		as->apBufNC[i] = (as->apBufNC[i] << 1) | (as->apNC[i] * 0x00000001);

		as->threshNC[i] = as->apNC[i] * threshMaxNC + (1.0-as->apNC[i]) * as->threshNC[i];
	}
}

void MZone::updatePCOut()
{
	for (int i = 0; i < num_bc; i++)
	{
		as->inputPCBC[i] = 0;
	}

	for (int i = 0; i < num_pc; i++)
	{
		for (int j = 0; j < num_p_pc_from_pc_to_bc; j++)
		{
			as->inputPCBC[cs->pPCfromPCtoBC[i][j]] += as->apPC[i];
		}
	}

	for (int i = 0; i < num_nc; i++)
	{
		for (int j = 0; j < num_p_nc_from_pc_to_nc; j++)
		{
			as->inputPCNC[i * num_p_nc_from_pc_to_nc + j] = as->apPC[cs->pNCfromPCtoNC[i][j]];
		}
	}
}

void MZone::updateBCPCOut()
{
	for (int i = 0; i < num_pc; i++) as->inputBCPC[i] = 0;
	
	for (int i = 0; i < num_bc; i++)
	{
		if (as->apBC[i]) 
		{
			for (int j = 0; j < num_p_bc_from_bc_to_pc; j++)
			{
				/* if there is a bc spike, obtain the pcs that this connects with
				 * and increment the input that the pc gets 
				 */
				as->inputBCPC[cs->pBCfromBCtoPC[i][j]]++;
			}
		}
	}
}

void MZone::updateSCPCOut()
{
	for (int i = 0; i < num_pc; i++) as->inputSCPC[i] = 0;

	for (int i = 0; i < num_sc; i++)
	{
		if (as->apSC[i])
		{
			for (int j = 0; j < num_p_sc_from_sc_to_pc; j++)
			{
				as->inputSCPC[cs->pSCfromSCtoPC[i][j]]++;
			}
		}
	}
}

void MZone::updateIOOut()
{
	for (int i = 0; i < num_io; i++)
	{
		as->pfPCPlastTimerIO[i] = (1 - as->apIO[i]) * (as->pfPCPlastTimerIO[i] + 1) + as->apIO[i] * tsLTPEndAPIO;
		as->vCoupleIO[i] = 0;
		for (int j = 0; j < num_p_io_in_io_to_io; j++)
		{
			as->vCoupleIO[i] += coupleRiRjRatioIO * (as->vIO[cs->pIOInIOIO[i][j]] - as->vIO[i]);
		}
	}
}

void MZone::updateNCOut()
{
	for (int i = 0; i < num_nc; i++)
	{
		as->synIOPReleaseNC[i] *= exp(-msPerTimeStep / 
				(relPDecTSofNCtoIO * exp(-as->synIOPReleaseNC[i] / relPDecTTofNCtoIO) +
				 relPDecT0ofNCtoIO));
		as->synIOPReleaseNC[i] += as->apNC[i] * relPIncNCtoIO *
				exp(-as->synIOPReleaseNC[i] / relPIncTauNCtoIO);
	}

	for (int i = 0; i < num_io; i++)
	{
		for (int j = 0; j < num_p_io_from_nc_to_io; j++)
		{
			as->inputNCIO[i * num_p_io_from_nc_to_io + j] = (randGen->Random() < as->synIOPReleaseNC[cs->pIOfromNCtoIO[i][j]]);
		}
	}
}

//void MZone::updateMFNCOut()
//{
//	for (int i = 0; i < num_nc; i++)
//	{
//		for (int j = 0; j < num_p_nc_from_mf_to_nc; j++)
//		{
//			as->inputMFNC[i * num_p_nc_from_mf_to_nc + j] = apMFInput[cs->pNCfromMFtoNC[i][j]];
//		}
//	}
//}

//void MZone::updateMFNCSyn(const uint8_t *histMF, uint32_t t)
//{
//
//#ifdef DEBUGOUT
//	float sumSynW = 0;
//#endif
//	if (t % (uint32_t)tsPerPopHistBinPC == 0) return;
//
//	as->histPCPopActSum = (as->histPCPopActSum) - (as->histPCPopAct[as->histPCPopActCurBinN]) + (as->pcPopAct);
//	as->histPCPopAct[as->histPCPopActCurBinN] = as->pcPopAct;
//	as->pcPopAct = 0;
//	as->histPCPopActCurBinN++;
//	as->histPCPopActCurBinN %= (uint32_t)numPopHistBinsPC;
//
//	float avgAllAPPC = ((float)as->histPCPopActSum) / numPopHistBinsPC;
//
//#ifdef DEBUGOUT
//	std::cout << "avgAllAPPC: " << avgAllAPPC << std::endl;
//#endif
//
//	bool doLTD = false;
//	bool doLTP = false;
//	if (avgAllAPPC >= synLTDPCPopActThreshMFtoNC && !as->noLTDMFNC)
//	{
//		doLTD = true;
//		as->noLTDMFNC = true;
//	}
//	else if (avgAllAPPC < synLTDPCPopActThreshMFtoNC)
//	{
//		as->noLTDMFNC = false;
//	}
//
//	if (avgAllAPPC <= synLTPPCPopActThreshMFtoNC && !as->noLTPMFNC)
//	{
//		doLTP = true;
//		as->noLTPMFNC = true;
//	}
//	else if (avgAllAPPC > synLTPPCPopActThreshMFtoNC)
//	{
//		as->noLTPMFNC = false;
//	}
//
//	for (int i = 0; i < num_nc; i++)
//	{
//		for(int j = 0; j < num_p_nc_from_mf_to_nc; j++)
//		{
//			float synWDelta = histMF[cs->pNCfromMFtoNC[i][j]] * (doLTD * synLTDStepSizeMFtoNC +
//					doLTP * synLTPStepSizeMFtoNC);
//			as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] += synWDelta;
//			as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] *= as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] > 0;
//			as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] *= as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] <= 1; 
//			as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] += as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] > 1;
//			
//			//Now uses isTrueMF to take collaterals into account
//			as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j] *= isTrueMF[cs->pNCfromMFtoNC[i][j]];
//#ifdef DEBUGOUT
//			sumSynW += as->mfSynWeightNC[i * num_p_nc_from_mf_to_nc + j];
//#endif
//		}
//	}
//#ifdef DEBUGOUT
//	std::cout << sumSynW / NUM_MF << std::endl;
//#endif
//}

void MZone::runUpdatePFPCOutCUDA(cudaStream_t **sts, int streamN)
{
	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		callUpdatePFPCOutKernel(sts[i][streamN], updatePFPCNumBlocks, updatePFPCNumGRPerB,
				num_pc, apBufGRGPU[i], gr_pc_con_in_d[i], delayMaskGRGPU[i], pfSynWeightPCGPU[i], inputPFPCGPU[i],
				inputPFPCGPUPitch[i]);
	}
}

void MZone::runSumPFPCCUDA(cudaStream_t **sts, int streamN)
{
	cudaError_t error;
	for (uint32_t i = 0; i < numGPUs; i++)
	{
		error=cudaSetDevice(i+gpuIndStart);
		callSumPFPCOutKernel(sts[i][streamN], sumPFPCOutNumBlocks, sumPFPCOutNumPCPerB,
			  updatePFPCNumBlocks, inputPFPCGPU[i], inputPFPCGPUPitch[i], inputSumPFPCMZGPU[i]);
	}
}

void MZone::cpyPFPCSumCUDA(cudaStream_t **sts, int streamN)
{
	for (uint32_t i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i+gpuIndStart);
		cudaMemcpyAsync(inputSumPFPCMZH[i], inputSumPFPCMZGPU[i], num_pc * sizeof(float),
			cudaMemcpyDeviceToHost, sts[i][streamN]);
	}
}

void MZone::runPFPCPlastCUDA(cudaStream_t **sts, int streamN, uint32_t t)
{
	cudaError_t error;
	if (t % (uint32_t)tsPerHistBinGR == 0)
	{
		int curGROffset;
		int curGPUInd;
		int curIOInd;

		int numGRPerIO;

		curGROffset = 0;
		curGPUInd   = 0;
		curIOInd    = 0;

		numGRPerIO = num_gr / num_io;

		for (int i = 0; i < num_io; i++)
		{
			if (as->pfPCPlastTimerIO[i] < (tsLTDstartAPIO + (int)tsLTDDurationIO) &&
					as->pfPCPlastTimerIO[i] >= tsLTDstartAPIO)
			{
				pfPCPlastStepIO[i] = tempGRPCLTDStep;
			}
			else if (as->pfPCPlastTimerIO[i] >= tsLTPstartAPIO ||
					as->pfPCPlastTimerIO[i] < tsLTPEndAPIO)
			{
				pfPCPlastStepIO[i] = tempGRPCLTPStep;
			}
			else
			{
				pfPCPlastStepIO[i] = 0;
			}
		}

		error = cudaSetDevice(curGPUInd + gpuIndStart);
		for (int i = 0; i < num_gr; i += num_p_pc_from_gr_to_pc)
		{
			if (i >= (curGPUInd + 1) * numGRPerGPU)
			{
				curGPUInd++;
				curGROffset = 0;
				error = cudaSetDevice(curGPUInd+gpuIndStart);
			}
			if (i >= (curIOInd + 1) * numGRPerIO)
			{
				curIOInd++;
			}
			callUpdatePFPCPlasticityIOKernel(sts[curGPUInd][(streamN + curIOInd) % 8],
					updatePFPCSynWNumBlocks, updatePFPCSynWNumGRPerB, pfSynWeightPCGPU[curGPUInd],
					histGRGPU[curGPUInd], grPCHistCheckBinIO, curGROffset, pfPCPlastStepIO[curIOInd]);

			curGROffset += num_p_pc_from_gr_to_pc;
		}
	}
}

void MZone::runSumPFSCCUDA(cudaStream_t **sts, int streamN)
{
	for(int i=0; i<numGPUs; i++)
	{
		cudaSetDevice(i+gpuIndStart);
		callSumPFBCSCOutKernel(sts[i][streamN], sumPFSCOutNumBlocks, sumPFSCOutNumSCPerB,
			  updatePFSCNumBlocks, inputPFSCGPU[i], inputPFSCGPUPitch[i], inputSumPFSCGPU[i]);
	}
}

void MZone::cpyPFSCSumCUDA(cudaStream_t **sts, int streamN)
{
	cudaError_t error;
	for(int i=0; i<numGPUs; i++)
	{
		cudaSetDevice(i+gpuIndStart);
		cudaMemcpyAsync(inputSumPFSCMZH[i], inputSumPFSCGPU[i], num_sc * sizeof(uint32_t),
			cudaMemcpyDeviceToHost, sts[i][streamN]);
	}
}

void MZone::runUpdatePFBCOutCUDA(cudaStream_t **sts, int streamN)
{
	cudaError_t error;
	for(int i=0; i<numGPUs; i++)
	{
		error=cudaSetDevice(i+gpuIndStart);
		callUpdatePFBCOutKernel(sts[i][streamN], updatePFBCNumBlocks, updatePFBCNumGRPerB,
			num_bc, apBufGRGPU[i], gr_bc_con_in_d[i], delayMaskGRGPU[i], inputPFBCGPU[i],
			inputPFBCGPUPitch[i]);
	}
}

void MZone::runUpdatePFSCOutCUDA(cudaStream_t **sts, int streamN)
{
	cudaError_t error;
	for(int i=0; i<numGPUs; i++)
	{
		error=cudaSetDevice(i+gpuIndStart);

		callUpdatePFSCOutKernel(sts[i][streamN], updatePFSCNumBlocks, updatePFSCNumGRPerB,
			num_sc, apBufGRGPU[i], gr_sc_con_in_d[i], delayMaskGRGPU[i], inputPFSCGPU[i],
			inputPFSCGPUPitch[i]);
	}
}

void MZone::runSumPFBCCUDA(cudaStream_t **sts, int streamN)
{
	cudaError_t error;
	for(int i=0; i<numGPUs; i++)
	{
		error=cudaSetDevice(i+gpuIndStart);

		callSumPFBCSCOutKernel(sts[i][streamN], sumPFBCOutNumBlocks, sumPFBCOutNumBCPerB,
			  updatePFBCNumBlocks, inputPFBCGPU[i], inputPFBCGPUPitch[i], inputSumPFBCGPU[i]);
	}
}

void MZone::cpyPFBCSumCUDA(cudaStream_t **sts, int streamN)
{
	cudaError_t error;
	for(int i=0; i<numGPUs; i++)
	{
		cudaSetDevice(i+gpuIndStart);
		cudaMemcpyAsync(inputSumPFBCMZH[i], inputSumPFBCGPU[i], num_bc * sizeof(uint32_t),
			cudaMemcpyDeviceToHost, sts[i][streamN]);
	}
}

void MZone::setGRPCPlastSteps(float ltdStep, float ltpStep)
{
	tempGRPCLTDStep = ltdStep;
	tempGRPCLTPStep = ltpStep;
}

void MZone::resetGRPCPlastSteps()
{
	tempGRPCLTDStep = synLTDStepSizeGRtoPC;
	tempGRPCLTPStep = synLTPStepSizeGRtoPC;
}

const float* MZone::exportPFPCWeights()
{
	cpyPFPCSynWCUDA();
	return (const float *)pfSynWeightPCLinear; 
}

//const float* MZone::exportMFDCNWeights()
//{
//	return (const float *)as->mfSynWeightNC.get(); 
//}

void MZone::load_pfpc_weights_from_file(std::fstream &in_file_buf)
{
	rawBytesRW((char *)pfSynWeightPCLinear, num_gr * sizeof(float), true, in_file_buf);
	
	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		int cpyStartInd = i * numGRPerGPU;
		cudaMemcpy(pfSynWeightPCGPU[i], &pfSynWeightPCLinear[cpyStartInd],
				numGRPerGPU * sizeof(float), cudaMemcpyHostToDevice);
	}
}

//void MZone::load_mfdcn_weights_from_file(std::fstream &in_file_buf)
//{
//	rawBytesRW((char *)as->mfSynWeightNC.get(), num_nc * num_p_nc_from_mf_to_nc * sizeof(float), true, in_file_buf);
//}

// Why not write one export function which takes in the thing you want to export?
const uint8_t* MZone::exportAPNC()
{
	return (const uint8_t *)as->apNC.get();
}

const uint8_t* MZone::exportAPSC()
{
	return (const uint8_t *)as->apSC.get();
}

const uint8_t* MZone::exportAPBC()
{
	return (const uint8_t *)as->apBC.get();
}

const uint8_t* MZone::exportAPPC()
{
	return (const uint8_t *)as->apPC.get();
}

const uint8_t* MZone::exportAPIO()
{
	return (const uint8_t *)as->apIO.get();
}

const float* MZone::exportgBCPC()
{
	return (const float *)as->gBCPC.get();
}

const float* MZone::exportgPFPC()
{
	return (const float *)as->gPFPC.get();
}

const float* MZone::exportVmBC()
{
	return (const float *)as->vBC.get();
}

const float* MZone::exportVmPC()
{
	return (const float *)as->vPC.get();
}

const float* MZone::exportVmNC()
{
	return (const float *)as->vNC.get();
}

const float* MZone::exportVmIO()
{
	return (const float *)as->vIO.get();
}

const unsigned int* MZone::exportAPBufBC()
{
	return (const unsigned int *)as->apBufBC.get();
}

const uint32_t* MZone::exportAPBufPC()
{
	return (const uint32_t *)as->apBufPC.get();
}

const uint8_t* MZone::exportAPBufIO()
{
	return (const uint8_t *)as->apBufIO.get();
}

const uint32_t* MZone::exportAPBufNC()
{
	return (const uint32_t *)as->apBufNC.get();
}

void MZone::testReduction()
{
	cudaError_t error;
	cudaStream_t *sts = new cudaStream_t[numGPUs];

	float hostTestData[num_gr] = {0.0};
	float hostPCSum[num_pc] = {0.0};
	float hostBCSum[num_bc] = {0.0};
	float hostSCSum[num_sc] = {0.0};

	float gpuToHostPCSum[num_pc] = {0.0};
	float gpuToHostBCSum[num_bc] = {0.0};
	float gpuToHostSCSum[num_sc] = {0.0};

	// leaving these dynamic for now as i do not understand cuda oof
	float **gpuPCTestData = new float*[numGPUs];
	float **gpuBCTestData = new float*[numGPUs];
	float **gpuSCTestData = new float*[numGPUs];

	size_t *gpuPCP = new size_t[numGPUs];
	size_t *gpuBCP = new size_t[numGPUs];
	size_t *gpuSCP = new size_t[numGPUs];

	float **gpuPCSum = new float*[numGPUs];
	float **gpuBCSum = new float*[numGPUs];
	float **gpuSCSum = new float*[numGPUs];

	for (int i = 0; i < num_gr; i++)
	{
		hostTestData[i] = randGen->Random();
	}

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);

		cudaStreamCreate(&sts[i]);

		cudaMallocPitch(&gpuPCTestData[i], &gpuPCP[i],
				num_p_pc_from_gr_to_pc * sizeof(float), num_pc / numGPUs);
		cudaMallocPitch(&gpuBCTestData[i], &gpuBCP[i],
				num_p_bc_from_gr_to_bc * sizeof(float), num_bc / numGPUs);
		cudaMallocPitch(&gpuSCTestData[i], &gpuSCP[i],
				num_p_sc_from_gr_to_sc*sizeof(float), numSCPerGPU);

		cudaMalloc(&gpuPCSum[i], num_pc / numGPUs * sizeof(float));
		cudaMalloc(&gpuBCSum[i], num_bc / numGPUs * sizeof(float));
		cudaMalloc(&gpuSCSum[i], numSCPerGPU * sizeof(float));

		LOG_DEBUG("Allocating memory for gpu %d", i);
		LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));

		cudaDeviceSynchronize();
	}

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);

		for (int j = 0; j < num_pc / numGPUs; j++)
		{
		   cudaMemcpy(((char *)gpuPCTestData[i] + j * gpuPCP[i]),
				 &hostTestData[i * numGRPerGPU + j * num_p_pc_from_gr_to_pc],
				 num_p_pc_from_gr_to_pc * sizeof(float), cudaMemcpyHostToDevice);
		}

		for (int j = 0; j < num_bc / numGPUs; j++)
		{
		   cudaMemcpy(((char *)gpuBCTestData[i] + j * gpuBCP[i]),
				 &hostTestData[i * numGRPerGPU + j * num_p_bc_from_gr_to_bc],
				 num_p_bc_from_gr_to_bc * sizeof(float), cudaMemcpyHostToDevice);
		}

		for (int j = 0; j < numSCPerGPU; j++)
		{
		   cudaMemcpy(((char *)gpuSCTestData[i] + j * gpuSCP[i]),
				 &hostTestData[i * numGRPerGPU + j * num_p_sc_from_gr_to_sc],
				 num_p_sc_from_gr_to_sc * sizeof(float), cudaMemcpyHostToDevice);
		}

		LOG_DEBUG("Copying memory for gpu %d", i);
		LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));

		cudaDeviceSynchronize();
	}

	for (int i = 0; i < num_pc; i++)
	{
		hostPCSum[i] = 0;

		for (int j = 0; j < num_p_pc_from_gr_to_pc; j++)
		{
			hostPCSum[i] += hostTestData[i * num_p_pc_from_gr_to_pc + j];
		}
	}

	for (int i = 0; i < num_bc; i++)
	{
		hostBCSum[i] = 0;

		for (int j = 0; j < num_p_bc_from_gr_to_bc; j++)
		{
			hostBCSum[i] += hostTestData[i * num_p_bc_from_gr_to_bc + j];
		}
	}

	for (int i = 0; i < num_sc; i++)
	{
		hostSCSum[i] = 0;

		for (int j = 0; j < num_p_sc_from_gr_to_sc; j++)
		{
			hostSCSum[i] += hostTestData[i * num_p_sc_from_gr_to_sc + j];
		}
	}

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		callSumKernel<float, true, false>(sts[i], gpuPCTestData[i], gpuPCP[i],
				gpuPCSum[i], 1, num_pc / numGPUs, 1, num_p_pc_from_gr_to_pc);

		callSumKernel<float, true, false>(sts[i], gpuBCTestData[i], gpuBCP[i],
				gpuBCSum[i], 1, num_bc / numGPUs, 1, num_p_bc_from_gr_to_bc);

		callSumKernel<float, true, false>(sts[i], gpuSCTestData[i], gpuSCP[i],
				gpuSCSum[i], 1, numSCPerGPU, 1, num_p_sc_from_gr_to_sc);

		cudaDeviceSynchronize();

		LOG_DEBUG("Calling sum kernels for gpu %d", i);
		LOG_DEBUG("Last error: %s", cudaGetErrorString(cudaGetLastError()));
	}

	for (int i = 0; i<numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);

		cudaMemcpy(&gpuToHostPCSum[i * num_pc / numGPUs], gpuPCSum[i],
			num_pc / numGPUs * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(&gpuToHostBCSum[i * num_bc / numGPUs], gpuBCSum[i],
			num_bc / numGPUs * sizeof(float), cudaMemcpyDeviceToHost);
		cudaMemcpy(&gpuToHostSCSum[i * numSCPerGPU], gpuSCSum[i],
			numSCPerGPU * sizeof(float), cudaMemcpyDeviceToHost);

		cudaDeviceSynchronize();
	}

	LOG_DEBUG("NumPC per GPU: %d", num_pc / numGPUs);
	LOG_DEBUG("NumBC per GPU: %d", num_bc / numGPUs);
	LOG_DEBUG("NumSC per GPU: %d", numSCPerGPU);

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i+gpuIndStart);
		cudaFree(gpuPCTestData[i]);
		cudaFree(gpuBCTestData[i]);
		cudaFree(gpuSCTestData[i]);
		cudaFree(gpuPCSum[i]);
		cudaFree(gpuBCSum[i]);
		cudaFree(gpuSCSum[i]);
		
	}

	for (int i = 0; i < numGPUs; i++)
	{
		cudaSetDevice(i + gpuIndStart);
		cudaStreamDestroy(sts[i]);
	}
	delete[] sts;

	delete[] gpuPCTestData;
	delete[] gpuBCTestData;
	delete[] gpuSCTestData;

	delete[] gpuPCP;
	delete[] gpuBCP;
	delete[] gpuSCP;

	delete[] gpuPCSum;
	delete[] gpuSCSum;
	delete[] gpuBCSum;
}

