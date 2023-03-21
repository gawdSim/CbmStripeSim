/*
 * activityparams.h
 *
 *  Created on: Oct 11, 2012
 *      Author: varicella
 */

#ifndef ACTIVITYPARAMS_H_
#define ACTIVITYPARAMS_H_

#include <math.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include "file_parse.h"
#include "stdint.h"

extern bool act_params_populated;

/* raw params */
extern float coupleRiRjRatioIO;
extern float eBCtoPC;
extern float eNCtoIO; 
extern float ePCtoBC;
extern float ePCtoNC;
extern float eSCtoPC;
extern float gDecTauBCtoPC;
extern float gIncBCtoPC;
extern float gDecT0ofNCtoIO;
extern float gDecTSofNCtoIO;
extern float gDecTTofNCtoIO;
extern float gIncNCtoIO;
extern float gIncTauNCtoIO;
extern float gDecTauPCtoBC;
extern float gDecTauPCtoNC;
extern float gIncAvgPCtoNC;
extern float gDecTauGRtoBC;
extern float gDecTauGRtoPC;
extern float gDecTauGRtoSC;
extern float gIncGRtoPC;
extern float gDecTauSCtoPC;
extern float gIncSCtoPC;
extern float synLTDStepSizeGRtoPC;
extern float synLTPStepSizeGRtoPC;
extern float maxExtIncVIO;
extern float msLTDDurationIO;
extern float msLTDStartAPIO;
extern float msLTPEndAPIO;
extern float msLTPStartAPIO;
extern float msPerHistBinGR;
extern float relPDecT0ofNCtoIO;
extern float relPDecTSofNCtoIO;
extern float relPDecTTofNCtoIO;
extern float relPIncNCtoIO;
extern float relPIncTauNCtoIO;
extern float gIncPCtoBC;
extern float gIncGRtoBC;
extern float gIncGRtoSC;
extern float rawGLeakBC;
extern float rawGLeakGR;
extern float rawGLeakIO;
extern float rawGLeakNC;
extern float rawGLeakPC;
extern float rawGLeakSC;
extern float threshDecTauBC;
//extern float threshDecTauGR;
extern float threshDecTauIO;
extern float threshDecTauNC;
extern float threshDecTauPC;
extern float threshDecTauSC;
extern float threshMaxBC;
//extern float threshMaxGR;
extern float threshMaxIO;
extern float threshMaxNC;
extern float threshMaxPC;
extern float threshMaxSC;
extern float weightScale;

/* derived act params */
//extern float threshDecGR;
extern float tsPerHistBinGR;
extern float gLeakSC;
extern float gDecGRtoSC;
extern float threshDecSC;
extern float gDecGRtoBC;
extern float gDecPCtoBC;
extern float threshDecBC;
extern float threshDecPC;
extern float gLeakPC;
extern float gDecGRtoPC;
extern float gDecBCtoPC;
extern float gDecSCtoPC;
extern float tsPopHistPC; /* used for updating MFNC syn plasticity */
extern float tsPerPopHistBinPC; /* used for updating MFNC syn plasticity */ 
// extern float numPopHistBinsPC; /* used for updating MFNC syn plasticity */ 
extern float gLeakIO;
extern float threshDecIO;
extern float tsLTDDurationIO;
extern float tsLTDstartAPIO;
extern float tsLTPstartAPIO;
extern float tsLTPEndAPIO; 
extern float grPCHistCheckBinIO; /* used in PFPC syn plasticity */
extern float gDecPCtoNC;
extern float gLeakNC;
extern float threshDecNC;
extern float gLeakBC;

void populate_act_params(parsed_sess_file &s_file);
void read_act_params(std::fstream &in_param_buf);
void write_act_params(std::fstream &out_param_buf);

#endif /* ACTIVITYPARAMS_H_ */

