/*
 * cbmstate.h
 *
 *  Created on: Dec 5, 2012
 *      Author: consciousness
 */

#ifndef M_ZONE_STATE_H_
#define M_ZONE_STATE_H_

#include <fstream>
#include <iostream>
#include <time.h>
#include <limits.h>

#include <cstdint>
#include "mzoneconnectivitystate.h"
#include "mzoneactivitystate.h"
#include "connectivityparams.h" // <-- added in 06/01/2022
#include "activityparams.h"

class MZoneState 
{
	public:
		MZoneState();
		MZoneState(unsigned int nZones);
		// TODO: make a choice which of two below constructors want to keep
		MZoneState(unsigned int nZones, std::fstream &sim_file_buf);
		~MZoneState();

		void readState(std::fstream &infile);
		void writeState(std::fstream &outfile);

		uint32_t getNumZones();

		MZoneConnectivityState *getMZoneConStateInternal(unsigned int zoneN);
		MZoneActivityState *getMZoneActStateInternal(unsigned int zoneN);

	private:
		uint32_t numZones;

		MZoneConnectivityState **mzoneConStates;
		MZoneActivityState **mzoneActStates;
};

#endif /* M_ZONE_STATE_H_ */

