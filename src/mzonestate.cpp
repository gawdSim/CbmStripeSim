/*
 * mzonestate.cpp
 *
 *  Created on: Jan 18, 2023
 *      Author: Self
 */

#include "logger.h"
#include "sfmt.h"
#include "mzonestate.h"

MZoneState::MZoneState() {}

MZoneState::MZoneState(unsigned int nZones) : numZones(nZones)
{
	LOG_DEBUG("Generating mzone state...");
	CRandomSFMT randGen(time(0));

	int *mzoneCRSeed = new int[nZones];
	int *mzoneARSeed = new int[nZones];

	mzoneConStates = new MZoneConnectivityState*[nZones];
	mzoneActStates = new MZoneActivityState*[nZones];
	for (int i = 0; i < nZones; i++)
	{
		mzoneCRSeed[i] = randGen.IRandom(0, INT_MAX);
		mzoneARSeed[i] = randGen.IRandom(0, INT_MAX);
		mzoneConStates[i] = new MZoneConnectivityState(mzoneCRSeed[i]);
		mzoneActStates[i] = new MZoneActivityState(mzoneARSeed[i]);
	}
	delete[] mzoneCRSeed;
	delete[] mzoneARSeed;
	LOG_DEBUG("Finished generating mzone state.");
}

MZoneState::MZoneState(unsigned int nZones, std::fstream &sim_file_buf) : numZones(nZones)
{
	LOG_DEBUG("Initializing cbm state from file...");
	mzoneConStates = new MZoneConnectivityState*[nZones];
	mzoneActStates = new MZoneActivityState*[nZones];

	for (int i = 0; i < nZones; i++)
	{
		mzoneConStates[i] = new MZoneConnectivityState(sim_file_buf);
		mzoneActStates[i] = new MZoneActivityState(sim_file_buf);
	}
	LOG_DEBUG("Finished initializing cbm state.");
}

MZoneState::~MZoneState()
{
	for (int i = 0; i < numZones; i++) 
	{
		delete mzoneConStates[i];
		delete mzoneActStates[i];
	}
	delete[] mzoneConStates;
	delete[] mzoneActStates;
}

void MZoneState::readState(std::fstream &infile)
{
	for (int i = 0; i < numZones; i++)
	{
		mzoneConStates[i]->readState(infile);
		mzoneActStates[i]->readState(infile);
	}
}

void MZoneState::writeState(std::fstream &outfile)
{
	for (int i = 0; i < numZones; i++)
	{
		mzoneConStates[i]->writeState(outfile);
		mzoneActStates[i]->writeState(outfile);
	}
}

uint32_t MZoneState::getNumZones() { return numZones; }

MZoneActivityState* MZoneState::getMZoneActStateInternal(unsigned int zoneN) { return mzoneActStates[zoneN]; }

MZoneConnectivityState* MZoneState::getMZoneConStateInternal(unsigned int zoneN) { return mzoneConStates[zoneN]; }

