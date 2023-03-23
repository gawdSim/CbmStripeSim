/*
 * cbmsimcore.cpp
 *
 *  Created on: Dec 15, 2011
 *      Author: consciousness
 */

#include "logger.h"
#include "cbmsimcore.h"

//#define NO_ASYNC
//#define DISP_CUDA_ERR

CBMSimCore::CBMSimCore() {}

CBMSimCore::CBMSimCore(std::fstream &psth_file_buf, MZoneState *state,
	int gpuIndStart, int numGPUP2)
{
	CRandomSFMT0 randGen(time(0));
	int *mzoneRSeed = new int[state->getNumZones()];

	for (int i = 0; i < state->getNumZones(); i++)
	{
		mzoneRSeed[i] = randGen.IRandom(0, INT_MAX);
	}

	construct(psth_file_buf, state, mzoneRSeed, gpuIndStart, numGPUP2);

	delete[] mzoneRSeed;
}

CBMSimCore::~CBMSimCore()
{
	if (grs) delete grs;

	for (int i = 0; i < numZones; i++)
	{
		delete zones[i];
	}

	delete[] zones;

	for (int i = 0; i < numGPUs; i++)
	{
		// How could gpuIndStart ever not be 0,
		// given we're looping from 0 to numGPUs?
		cudaSetDevice(i + gpuIndStart);

		for (int j = 0; j < 8; j++)
		{
			cudaStreamDestroy(streams[i][j]);
		}
		delete[] streams[i];
	}

	delete[] streams;
}

// for speed
void CBMSimCore::writeToState()
{
	for (int i = 0; i < numZones; i++)
	{
		zones[i]->writeToState();
	}
}

void CBMSimCore::writeState(std::fstream& outfile)
{
	writeToState();
	simState->writeState(outfile); // using internal cp
}

void CBMSimCore::initCUDAStreams()
{
	cudaError_t error;

	int maxNumGPUs;
	// TODO: use assert, try, and catch for these types of errors
	error = cudaGetDeviceCount(&maxNumGPUs);

	LOG_DEBUG("CUDA max num devices: %d", maxNumGPUs);
	LOG_DEBUG("%s", cudaGetErrorString(error));
	LOG_DEBUG("CUDA num devices: %d, starting at GPU %d", numGPUs, gpuIndStart);

	streams = new cudaStream_t*[numGPUs];

	for (int i = 0; i < numGPUs; i++)
	{
		error = cudaSetDevice(i + gpuIndStart);
		LOG_DEBUG("Selecting device %d", i);
		LOG_DEBUG("%s", cudaGetErrorString(error));
		streams[i] = new cudaStream_t[8];
		LOG_DEBUG("Resetting device %d", i);
		LOG_DEBUG("%s", cudaGetErrorString(error));
		cudaDeviceSynchronize();

		for (int j = 0; j < 8; j++)
		{
			error = cudaStreamCreate(&streams[i][j]);
			LOG_DEBUG("Initializing stream %d for device %d",j, i);
			LOG_DEBUG("%s", cudaGetErrorString(error));
		}
		cudaDeviceSynchronize();
		error = cudaGetLastError();
		LOG_DEBUG("Cuda device %d", i);
		LOG_DEBUG("%s", cudaGetErrorString(error));
	}
}

void CBMSimCore::syncCUDA(std::string title)
{
	cudaError_t error;
	for (int i = 0; i < numGPUs; i++)
	{
		error = cudaSetDevice(i + gpuIndStart);
#ifdef DISP_CUDA_ERR
		LOG_TRACE("sync point  %s, switching to gpu %d", title.c_str(), i);
		LOG_TRACE("%s", cudaGetErrorString(error));
#endif
		error = cudaDeviceSynchronize();
#ifdef DISP_CUDA_ERR
		LOG_TRACE("sync point  %s, switching to gpu %d", title.c_str(), i);
		LOG_TRACE("%s", cudaGetErrorString(error));
#endif
	}
}

void CBMSimCore::calcActivity(enum plasticity pf_pc_plast, uint32_t ts)
{
	syncCUDA("1");

	grs->calcGRPoissActivity(ts, streams, 5);

	for (int i = 0; i < numZones; i++)
	{
		if (pf_pc_plast == GRADED)
		{
			zones[i]->runPFPCPlastCUDA(streams, 1, ts);
		}
		zones[i]->runUpdatePFPCOutCUDA(streams, i + 2);
		zones[i]->runSumPFPCCUDA(streams, i + 1);
		zones[i]->cpyPFPCSumCUDA(streams, i);

		zones[i]->runUpdatePFBCOutCUDA(streams, i+4);
		zones[i]->runSumPFBCCUDA(streams, 2);
		zones[i]->cpyPFBCSumCUDA(streams, i);
	
		zones[i]->runUpdatePFSCOutCUDA(streams, i+4);
		zones[i]->runSumPFSCCUDA(streams, 3);
		zones[i]->cpyPFSCSumCUDA(streams, i);

		zones[i]->calcBCActivities();
		zones[i]->calcSCActivities();

		zones[i]->updateBCPCOut();
		zones[i]->updateSCPCOut();

		zones[i]->calcPCActivities();
		zones[i]->updatePCOut();

		zones[i]->calcNCActivities();
		zones[i]->updateNCOut();

		zones[i]->calcIOActivities();
		zones[i]->updateIOOut();
	}
}

void CBMSimCore::updateErrDrive(unsigned int zoneN, float errDriveRelative)
{
	zones[zoneN]->setErrDrive(errDriveRelative);
}

MZone** CBMSimCore::getMZoneList()
{
	return (MZone **)zones;
}

void CBMSimCore::construct(std::fstream &psth_file_buf, MZoneState *state,
	int *mzoneRSeed, int gpuIndStart, int numGPUP2)
{
	int maxNumGPUs;

	numZones = state->getNumZones();

	cudaGetDeviceCount(&maxNumGPUs);

	if (gpuIndStart <= 0)
	{
		this->gpuIndStart = 0;
	}
	else if (gpuIndStart >= maxNumGPUs)
	{
		this->gpuIndStart = maxNumGPUs - 1;
	}
	else
	{
		this->gpuIndStart = gpuIndStart;
	}

	if (numGPUP2 < 0)
	{
		numGPUs = maxNumGPUs;
	}
	else
	{
		numGPUs = (unsigned int)numGPUP2;
	}

	if (this->gpuIndStart + numGPUs > maxNumGPUs)
	{
		numGPUs = 1;
	}
	LOG_DEBUG("Calculated (?) number of GPUs: %d", numGPUs);

	LOG_DEBUG("Initializing cuda streams...");
	initCUDAStreams();
	LOG_DEBUG("Finished initialzing cuda streams.");

	LOG_DEBUG("Initializing granule cells ...");
	grs = new PoissonRegenCells(psth_file_buf, streams);
	LOG_DEBUG("granule cells initialized.");

	LOG_DEBUG("Initializing mzones...");
	zones = new MZone*[numZones];
	for (int i = 0; i < numZones; i++)
	{
		// same thing for zones as with innet
		zones[i] = new MZone(state->getMZoneConStateInternal(i),
			state->getMZoneActStateInternal(i), grs->get_ap_buf_gr_gpu(),
			grs->get_ap_hist_gr_gpu(), mzoneRSeed[i], this->gpuIndStart, numGPUs);
	}
	LOG_DEBUG("Mzone construction complete");

	simState = state; // shallow copy
}

