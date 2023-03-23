/*
 * mzoneconnectivitystate.cpp
 *
 *  Created on: Nov 21, 2012
 *      Author: consciousness
 */

#include <algorithm>
#include <iostream>
#include <limits.h>

#include "logger.h"
#include "array_util.h"
#include "dynamic2darray.h"
#include "sfmt.h"
#include "connectivityparams.h"
#include "mzoneconnectivitystate.h"

MZoneConnectivityState::MZoneConnectivityState(int randSeed)
{
	CRandomSFMT0 randGen(randSeed);
	LOG_DEBUG("Allocating and initializing mzone connectivity arrays...");
	allocateMemory();
	//initializeVals();
	LOG_DEBUG("Initializing mzone connections...");
	LOG_DEBUG("Assigning GR delays");
	assignGRDelays();
	LOG_DEBUG("Connecting GR to PC");
	connectGRtoPC();
	LOG_DEBUG("Connecting GR to BC");
	connectGRtoBC();
	LOG_DEBUG("Connecting BC to PC");
	connectBCtoPC();
	LOG_DEBUG("Connecting PC to BC");
	connectPCtoBC();
	LOG_DEBUG("Connecting SC and PC");
	connectSCtoPC();
	LOG_DEBUG("Connecting PC and NC");
	connectPCtoNC(randGen);
	LOG_DEBUG("Connecting NC and IO");
	connectNCtoIO();
	LOG_DEBUG("Connecting IO and PC");
	connectIOtoPC();
	LOG_DEBUG("Connecting IO and IO");
	connectIOtoIO();
	LOG_DEBUG("Finished making mzone connections.");
}

MZoneConnectivityState::MZoneConnectivityState(std::fstream &infile)
{
	allocateMemory();
	stateRW(true, infile);
}

MZoneConnectivityState::~MZoneConnectivityState() {deallocMemory();}

void MZoneConnectivityState::readState(std::fstream &infile)
{
	stateRW(true, infile);
}

void MZoneConnectivityState::writeState(std::fstream &outfile)
{
	stateRW(false, outfile);
}

void MZoneConnectivityState::allocateMemory()
{
	//granule cells
	pGRDelayMaskfromGRtoBSP = new uint32_t[num_gr];

	pGRfromGRtoPC = new uint32_t[num_gr]();
	pGRfromGRtoBC = new uint32_t[num_gr]();

	//basket cells
	pBCfromBCtoPC  = allocate2DArray<uint32_t>(num_bc, num_p_bc_from_bc_to_pc);
	pBCfromPCtoBC  = allocate2DArray<uint32_t>(num_bc, num_p_bc_from_pc_to_bc);

	//stellate cells
	pSCfromSCtoPC = allocate2DArray<uint32_t>(num_sc, num_p_sc_from_sc_to_pc);

	//purkinje cells
	pPCfromBCtoPC = allocate2DArray<uint32_t>(num_pc, num_p_pc_from_bc_to_pc);
	pPCfromPCtoBC = allocate2DArray<uint32_t>(num_pc, num_p_pc_from_pc_to_bc);
	pPCfromSCtoPC = allocate2DArray<uint32_t>(num_pc, num_p_pc_from_sc_to_pc);
	pPCfromPCtoNC = allocate2DArray<uint32_t>(num_pc, num_p_pc_from_pc_to_nc);
	pPCfromIOtoPC = new uint32_t[num_pc]();

	//nucleus cells
	pNCfromPCtoNC = allocate2DArray<uint32_t>(num_nc, num_p_nc_from_pc_to_nc);
	pNCfromNCtoIO = allocate2DArray<uint32_t>(num_nc, num_p_nc_from_nc_to_io);

	//inferior olivary cells
	pIOfromIOtoPC = allocate2DArray<uint32_t>(num_io, num_p_io_from_io_to_pc);
	pIOfromNCtoIO = allocate2DArray<uint32_t>(num_io, num_p_io_from_nc_to_io);
	pIOInIOIO = allocate2DArray<uint32_t>(num_io, num_p_io_in_io_to_io);
	pIOOutIOIO = allocate2DArray<uint32_t>(num_io, num_p_io_out_io_to_io);
}

//void MZoneConnectivityState::initializeVals()
//{
//	// granule cells
//	std::fill(pGRDelayMaskfromGRtoBSP, pGRDelayMaskfromGRtoBSP + num_gr, 0);
//
//	// basket cells
//	std::fill(pBCfromBCtoPC[0], pBCfromBCtoPC[0]
//			+ num_bc * num_p_bc_from_bc_to_pc, 0);
//	std::fill(pBCfromPCtoBC[0], pBCfromPCtoBC[0]
//			+ num_bc * num_p_bc_from_pc_to_bc, 0);
//
//	// stellate cells
//	std::fill(pSCfromSCtoPC[0], pSCfromSCtoPC[0]
//			+ num_sc * num_p_sc_from_sc_to_pc, 0);
//
//	// purkinje cells
//	std::fill(pPCfromBCtoPC[0], pPCfromBCtoPC[0]
//			+ num_pc * num_p_pc_from_bc_to_pc, 0);
//	std::fill(pPCfromPCtoBC[0], pPCfromPCtoBC[0]
//			+ num_pc * num_p_pc_from_pc_to_bc, 0);
//	std::fill(pPCfromSCtoPC[0], pPCfromSCtoPC[0]
//			+ num_pc * num_p_pc_from_sc_to_pc, 0);
//	std::fill(pPCfromPCtoNC[0], pPCfromPCtoNC[0]
//			+ num_pc * num_p_pc_from_pc_to_nc, 0);
//
//	// nucleus cells
//	std::fill(pNCfromPCtoNC[0], pNCfromPCtoNC[0]
//			+ num_nc * num_p_nc_from_pc_to_nc, 0);
//	std::fill(pNCfromNCtoIO[0], pNCfromNCtoIO[0]
//			+ num_nc * num_p_nc_from_nc_to_io, 0);
//
//	// inferior olivary cells
//	std::fill(pIOfromIOtoPC[0], pIOfromIOtoPC[0]
//			+ num_io * num_p_io_from_io_to_pc, 0);
//	std::fill(pIOfromNCtoIO[0], pIOfromNCtoIO[0]
//			+ num_io * num_p_io_from_nc_to_io, 0);
//	std::fill(pIOInIOIO[0], pIOInIOIO[0]
//			+ num_io * num_p_io_in_io_to_io, 0);
//	std::fill(pIOOutIOIO[0], pIOOutIOIO[0]
//			+ num_io * num_p_io_out_io_to_io, 0);
//}

void MZoneConnectivityState::deallocMemory()
{
	// granule cells
	delete[] pGRDelayMaskfromGRtoBSP;

	delete[] pGRfromGRtoPC;
	delete[] pGRfromGRtoBC;

	// basket cells
	delete2DArray<uint32_t>(pBCfromBCtoPC);
	delete2DArray<uint32_t>(pBCfromPCtoBC);

	// stellate cells
	delete2DArray<uint32_t>(pSCfromSCtoPC);

	// purkinje cells
	delete2DArray<uint32_t>(pPCfromBCtoPC);
	delete2DArray<uint32_t>(pPCfromPCtoBC);
	delete2DArray<uint32_t>(pPCfromSCtoPC);
	delete2DArray<uint32_t>(pPCfromPCtoNC);
	delete[] pPCfromIOtoPC;

	// nucleus cells
	delete2DArray<uint32_t>(pNCfromPCtoNC);
	delete2DArray<uint32_t>(pNCfromNCtoIO);

	// inferior olivary cells
	delete2DArray<uint32_t>(pIOfromIOtoPC);
	delete2DArray<uint32_t>(pIOfromNCtoIO);
	delete2DArray<uint32_t>(pIOInIOIO);
	delete2DArray<uint32_t>(pIOOutIOIO);
}

void MZoneConnectivityState::stateRW(bool read, std::fstream &file)
{
	// granule cells
	rawBytesRW((char *)pGRDelayMaskfromGRtoBSP, num_gr * sizeof(uint32_t), read, file);

	// basket cells
	rawBytesRW((char *)pBCfromBCtoPC[0], num_bc * num_p_bc_from_bc_to_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pBCfromPCtoBC[0], num_bc * num_p_bc_from_pc_to_bc * sizeof(uint32_t), read, file);

	// stellate cells
	rawBytesRW((char *)pSCfromSCtoPC[0], num_sc * num_p_sc_from_sc_to_pc * sizeof(uint32_t), read, file);

	// purkinje cells
	rawBytesRW((char *)pPCfromBCtoPC[0], num_pc * num_p_pc_from_bc_to_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pPCfromPCtoBC[0], num_pc * num_p_pc_from_pc_to_bc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pPCfromSCtoPC[0], num_pc * num_p_pc_from_sc_to_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pPCfromPCtoNC[0], num_pc * num_p_pc_from_pc_to_nc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pPCfromIOtoPC, num_pc * sizeof(uint32_t), read, file);

	// nucleus cells
	rawBytesRW((char *)pNCfromPCtoNC[0], num_nc * num_p_nc_from_pc_to_nc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pNCfromNCtoIO[0], num_nc * num_p_nc_from_nc_to_io * sizeof(uint32_t), read, file);

	// inferior olivary cells
	rawBytesRW((char *)pIOfromIOtoPC[0], num_io * num_p_io_from_io_to_pc * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pIOfromNCtoIO[0], num_io * num_p_io_from_nc_to_io * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pIOInIOIO[0], num_io * num_p_io_in_io_to_io * sizeof(uint32_t), read, file);
	rawBytesRW((char *)pIOOutIOIO[0], num_io * num_p_io_out_io_to_io * sizeof(uint32_t), read, file);
}

void MZoneConnectivityState::assignGRDelays()
{
	for (int i = 0; i < num_gr; i++)
	{
		// calculate x coordinate of GR position
		int grPosX = i % gr_x;

		// calculate distance of GR (assume soma) to BC, PC, and SC (aa + pf distance)
		// and assign time delay.
		int grBCPCSCDist = abs((int)(gr_x / 2 - grPosX));
		pGRDelayMaskfromGRtoBSP[i] = 0x1 << (int)((grBCPCSCDist / gr_pf_vel_in_gr_x_per_t_step
			+ gr_af_delay_in_t_step) / msPerTimeStep);
		
	}
}

//TODO: TEST
void MZoneConnectivityState::connectGRtoPC() {
	for (uint32_t i = 0; i < num_gr; i++) {
		pGRfromGRtoPC[i] = i % num_pc;
	}
	fisher_yates_shuffle<uint32_t>(pGRfromGRtoPC, num_gr);
}

void MZoneConnectivityState::connectGRtoBC() {
	for (uint32_t i = 0; i < num_gr; i++) {
		pGRfromGRtoBC[i] = i % num_bc;
	}
	fisher_yates_shuffle<uint32_t>(pGRfromGRtoBC, num_gr);
}


void MZoneConnectivityState::connectBCtoPC()
{
	int bcToPCRatio = num_bc / num_pc;

	for (int i = 0; i < num_pc; i++)
	{
		pBCfromBCtoPC[i * bcToPCRatio][0] = ((i + 1) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio][1] = ((i - 1) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio][2] = ((i + 2) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio][3] = ((i - 2) % num_pc + num_pc) % num_pc;

		pBCfromBCtoPC[i * bcToPCRatio+1][0] = ((i + 1) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+1][1] = ((i - 1) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+1][2] = ((i + 3) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+1][3] = ((i - 3) % num_pc + num_pc) % num_pc;

		pBCfromBCtoPC[i * bcToPCRatio+2][0] = ((i + 3) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+2][1] = ((i - 3) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+2][2] = ((i + 6) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+2][3] = ((i - 6) % num_pc + num_pc) % num_pc;

		pBCfromBCtoPC[i * bcToPCRatio+3][0] = ((i + 4) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+3][1] = ((i - 4) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+3][2] = ((i + 9) % num_pc + num_pc) % num_pc;
		pBCfromBCtoPC[i * bcToPCRatio+3][3] = ((i - 9) % num_pc + num_pc) % num_pc;
	}
}

void MZoneConnectivityState::connectPCtoBC()
{
	int bcToPCRatio = num_bc / num_pc;

	for (int i = 0; i < num_pc; i++)
	{
		for (int j = 0; j < num_p_pc_from_pc_to_bc; j++)
		{
			int bcInd = i * bcToPCRatio - 6 + j;
			pPCfromPCtoBC[i][j] = (bcInd % num_bc + num_bc) % num_bc;
		}
	}
}

void MZoneConnectivityState::connectSCtoPC()
{
	for (int i = 0; i < num_sc; i++)
	{
		for (int j = 0; j < num_p_sc_from_sc_to_pc; j++)
		{
			int pcInd = i / num_p_pc_from_sc_to_pc;
			pSCfromSCtoPC[i][j] = pcInd;
			pPCfromSCtoPC[pcInd][i % num_p_pc_from_sc_to_pc] = i;
		}
	}
}

void MZoneConnectivityState::connectPCtoNC(CRandomSFMT0 &randGen)
{
	int pcNumConnected[num_pc] = {0};
	std::fill(pcNumConnected, pcNumConnected + num_pc, 1);
	
	for (int i = 0; i < num_nc; i++)
	{
		for (int j = 0; j < num_pc / num_nc; j++)
		{
			int pcInd = i * (num_pc / num_nc) + j;
			pNCfromPCtoNC[i][j] = pcInd;
			pPCfromPCtoNC[pcInd][0] = i;
		}
	}

	for (int i = 0; i < num_nc - 1; i++)
	{
		for (int j = num_p_nc_from_pc_to_nc / 3; j < num_p_nc_from_pc_to_nc; j++)
		{
			int countPCNC = 0;
			int pcInd;
			while(true)
			{
				bool connect = true;

				pcInd = randGen.IRandomX(0, num_pc - 1);

				if (pcNumConnected[pcInd] >= num_p_pc_from_pc_to_nc) continue;
				for (int k = 0; k < pcNumConnected[pcInd]; k++)
				{
					if (pPCfromPCtoNC[pcInd][k] == i)
					{
						connect = false;
						break;
					}
				}
				if (connect || countPCNC > 100) break;
				countPCNC++;
			}

			pNCfromPCtoNC[i][j] = pcInd;
			pPCfromPCtoNC[pcInd][pcNumConnected[pcInd]] = i;

			pcNumConnected[pcInd]++;
		}
	}

	// static cast?
	unsigned int numSyn = num_pc / num_nc;

	for (int h = 1; h < num_p_pc_from_pc_to_nc; h++)
	{
		for(int i = 0; i < num_pc; i++)
		{
			if (pcNumConnected[i] < num_p_pc_from_pc_to_nc)
			{
				pNCfromPCtoNC[num_nc - 1][numSyn] = i;
				pPCfromPCtoNC[i][pcNumConnected[i]] = num_nc - 1;

				pcNumConnected[i]++;
				numSyn++;
			}
		}
	}
}

void MZoneConnectivityState::connectNCtoIO()
{
	for (int i = 0; i < num_io; i++)
	{
		for (int j = 0; j < num_p_io_from_nc_to_io; j++)
		{
			pIOfromNCtoIO[i][j] = j;
			pNCfromNCtoIO[j][i] = i;
		}
	}
}

void MZoneConnectivityState::connectIOtoPC()
{
	for (int i = 0; i < num_io; i++)
	{
		for (int j = 0; j < num_p_io_from_io_to_pc; j++)
		{
			int pcInd = i * num_p_io_from_io_to_pc + j;

			pIOfromIOtoPC[i][j]  = pcInd;
			pPCfromIOtoPC[pcInd] = i;
		}
	}
}

void MZoneConnectivityState::connectIOtoIO()
{
	for (int i = 0; i < num_io; i++)
	{
		int inInd = 0;
		for (int j = 0; j < num_p_io_in_io_to_io; j++)
		{
			if (inInd == i) inInd++;
			pIOInIOIO[i][j] = inInd;
			inInd++;
		}
	}
}

