/*
 * mzoneconnectivitystate.h
 *
 *  Created on: Nov 21, 2012
 *      Author: consciousness
 */

#ifndef MZONECONNECTIVITYSTATE_H_
#define MZONECONNECTIVITYSTATE_H_

#include <fstream>
#include <cstdint>

class MZoneConnectivityState
{
public:
	MZoneConnectivityState();
	MZoneConnectivityState(int randSeed);
	MZoneConnectivityState(std::fstream &infile);
	~MZoneConnectivityState();

	void readState(std::fstream &infile);
	void writeState(std::fstream &outfile);

	//granule cells
	uint32_t *pGRDelayMaskfromGRtoBSP;
	uint32_t *pGRfromGRtoPC; // index on gr side of what pc that gr connects to
	uint32_t *pGRfromGRtoBC; // index on gr side of what bc that gr connects to
	uint32_t *pGRfromGRtoSC; // index on gr side of what sc that gr connects to

	//basket cells
	uint32_t **pBCfromBCtoPC;
	uint32_t **pBCfromPCtoBC;

	//stellate cells
	uint32_t **pSCfromSCtoPC;

	//purkinje cells
	uint32_t **pPCfromBCtoPC;
	uint32_t **pPCfromPCtoBC;
	uint32_t **pPCfromSCtoPC;
	uint32_t **pPCfromPCtoNC;
	uint32_t *pPCfromIOtoPC;

	//nucleus cells
	uint32_t **pNCfromPCtoNC;
	uint32_t **pNCfromNCtoIO;

	//inferior olivary cells
	uint32_t **pIOfromIOtoPC;
	uint32_t **pIOfromNCtoIO;
	uint32_t **pIOInIOIO;
	uint32_t **pIOOutIOIO;

private:
	void allocateMemory();
	void initializeVals();
	void deallocMemory();
	void stateRW(bool read, std::fstream &file);

	void assignGRDelays();
	void connectGRtoPC();
	void connectGRtoBC();
	void connectGRtoSC();
	void connectBCtoPC();
	void connectPCtoBC();
	void connectSCtoPC();
	void connectPCtoNC(CRandomSFMT0 &randGen);
	void connectNCtoIO();
	void connectIOtoPC();
	void connectIOtoIO();
};

#endif /* MZONECONNECTIVITYSTATE_H_ */

