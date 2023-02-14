/*
 * activityparams.cpp
 *
 *  Created on: Oct 11, 2012
 *      Author: varicella
 */

#include <fstream>
#include <cstring>
#include <assert.h>

#include "file_utility.h"
#include "connectivityparams.h"
#include "activityparams.h"

bool act_params_populated = false;

float coupleRiRjRatioIO          = 0.0;
float eBCtoPC                    = 0.0;
float eNCtoIO                    = 0.0;
float ePCtoBC                    = 0.0;
float ePCtoNC                    = 0.0;
float eSCtoPC                    = 0.0;
float gDecTauBCtoPC              = 0.0;
float gIncBCtoPC                 = 0.0;
float gDecT0ofNCtoIO             = 0.0;
float gDecTSofNCtoIO             = 0.0;
float gDecTTofNCtoIO             = 0.0;
float gIncNCtoIO                 = 0.0;
float gIncTauNCtoIO              = 0.0;
float gDecTauPCtoBC              = 0.0;
float gDecTauPCtoNC              = 0.0;
float gIncAvgPCtoNC              = 0.0;
float gDecTauGRtoBC              = 0.0;
float gDecTauGRtoPC              = 0.0;
float gDecTauGRtoSC              = 0.0;
float gIncGRtoPC                 = 0.0;
float gDecTauSCtoPC              = 0.0;
float gIncSCtoPC                 = 0.0;
float synLTDStepSizeGRtoPC       = 0.0; 
float synLTPStepSizeGRtoPC       = 0.0; 
float maxExtIncVIO               = 0.0;
float msLTDDurationIO            = 0.0;
float msLTDStartAPIO             = 0.0;
float msLTPEndAPIO               = 0.0;
float msLTPStartAPIO             = 0.0;
float msPerHistBinGR             = 0.0;
float relPDecT0ofNCtoIO          = 0.0;
float relPDecTSofNCtoIO          = 0.0; 
float relPDecTTofNCtoIO          = 0.0; 
float relPIncNCtoIO              = 0.0;
float relPIncTauNCtoIO           = 0.0; 
float gIncPCtoBC                 = 0.0;
float gIncGRtoBC                 = 0.0;
float gIncGRtoSC                 = 0.0;
float rawGLeakBC                 = 0.0;
float rawGLeakGO                 = 0.0;
float rawGLeakGR                 = 0.0;
float rawGLeakIO                 = 0.0;
float rawGLeakNC                 = 0.0;
float rawGLeakPC                 = 0.0;
float rawGLeakSC                 = 0.0;
float threshDecTauBC             = 0.0;
//float threshDecTauGR             = 0.0;
float threshDecTauIO             = 0.0;
float threshDecTauNC             = 0.0;
float threshDecTauPC             = 0.0;
float threshDecTauSC             = 0.0;
float threshMaxBC                = 0.0;
//float threshMaxGR                = 0.0;
float threshMaxIO                = 0.0;
float threshMaxNC                = 0.0;
float threshMaxPC                = 0.0;
float threshMaxSC                = 0.0;
float weightScale                = 0.0; 

/* derived act params */
//float threshDecGR         = 0.0; 
float tsPerHistBinGR      = 0.0; 
float gLeakSC             = 0.0; 
float gDecGRtoSC          = 0.0; 
float threshDecSC         = 0.0; 
float gDecGRtoBC          = 0.0; 
float gDecPCtoBC          = 0.0; 
float threshDecBC         = 0.0; 
float threshDecPC         = 0.0; 
float gLeakPC             = 0.0; 
float gDecGRtoPC          = 0.0; 
float gDecBCtoPC          = 0.0; 
float gDecSCtoPC          = 0.0; 
float tsPopHistPC         = 0.0; 
float tsPerPopHistBinPC   = 0.0; 
//float numPopHistBinsPC    = 0.0; 
float gLeakIO             = 0.0; 
float threshDecIO         = 0.0; 
float tsLTDDurationIO     = 0.0; 
float tsLTDstartAPIO      = 0.0; 
float tsLTPstartAPIO      = 0.0; 
float tsLTPEndAPIO        = 0.0; 
float grPCHistCheckBinIO  = 0.0; 
float gDecPCtoNC          = 0.0; 
float gLeakNC             = 0.0; 
float threshDecNC         = 0.0; 
float gLeakBC             = 0.0; 

void populate_act_params(parsed_sess_file &s_file)
{
	coupleRiRjRatioIO          = std::stof(s_file.parsed_var_sections["activity"].param_map["coupleRiRjRatioIO"].value); 
	eBCtoPC                    = std::stof(s_file.parsed_var_sections["activity"].param_map["eBCtoPC"].value); 
	eNCtoIO                    = std::stof(s_file.parsed_var_sections["activity"].param_map["eNCtoIO"].value); 
	ePCtoBC                    = std::stof(s_file.parsed_var_sections["activity"].param_map["ePCtoBC"].value); 
	ePCtoNC                    = std::stof(s_file.parsed_var_sections["activity"].param_map["ePCtoNC"].value); 
	eSCtoPC                    = std::stof(s_file.parsed_var_sections["activity"].param_map["eSCtoPC"].value); 
	gDecTauBCtoPC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauBCtoPC"].value); 
	gIncBCtoPC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncBCtoPC"].value); 
	gDecT0ofNCtoIO             = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecT0ofNCtoIO"].value); 
	gDecTSofNCtoIO             = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTSofNCtoIO"].value); 
	gDecTTofNCtoIO             = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTTofNCtoIO"].value); 
	gIncNCtoIO                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncNCtoIO"].value); 
	gIncTauNCtoIO              = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncTauNCtoIO"].value); 
	gDecTauPCtoBC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauPCtoBC"].value); 
	gDecTauPCtoNC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauPCtoNC"].value); 
	gIncAvgPCtoNC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncAvgPCtoNC"].value); 
	gDecTauGRtoBC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauGRtoBC"].value); 
	gDecTauGRtoPC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauGRtoPC"].value); 
	gDecTauGRtoSC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauGRtoSC"].value); 
	gIncGRtoPC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncGRtoPC"].value); 
	gDecTauSCtoPC              = std::stof(s_file.parsed_var_sections["activity"].param_map["gDecTauSCtoPC"].value); 
	gIncSCtoPC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncSCtoPC"].value);  
	synLTDStepSizeGRtoPC       = std::stof(s_file.parsed_var_sections["activity"].param_map["synLTDStepSizeGRtoPC"].value); 
	synLTPStepSizeGRtoPC       = std::stof(s_file.parsed_var_sections["activity"].param_map["synLTPStepSizeGRtoPC"].value); 
	maxExtIncVIO               = std::stof(s_file.parsed_var_sections["activity"].param_map["maxExtIncVIO"].value); 
	msLTDDurationIO            = std::stof(s_file.parsed_var_sections["activity"].param_map["msLTDDurationIO"].value); 
	msLTDStartAPIO             = std::stof(s_file.parsed_var_sections["activity"].param_map["msLTDStartAPIO"].value); 
	msLTPEndAPIO               = std::stof(s_file.parsed_var_sections["activity"].param_map["msLTPEndAPIO"].value); 
	msLTPStartAPIO             = std::stof(s_file.parsed_var_sections["activity"].param_map["msLTPStartAPIO"].value); 
	msPerHistBinGR             = std::stof(s_file.parsed_var_sections["activity"].param_map["msPerHistBinGR"].value); 
	relPDecT0ofNCtoIO          = std::stof(s_file.parsed_var_sections["activity"].param_map["relPDecT0ofNCtoIO"].value); 
	relPDecTSofNCtoIO          = std::stof(s_file.parsed_var_sections["activity"].param_map["relPDecTSofNCtoIO"].value); 
	relPDecTTofNCtoIO          = std::stof(s_file.parsed_var_sections["activity"].param_map["relPDecTTofNCtoIO"].value); 
	relPIncNCtoIO              = std::stof(s_file.parsed_var_sections["activity"].param_map["relPIncNCtoIO"].value); 
	relPIncTauNCtoIO           = std::stof(s_file.parsed_var_sections["activity"].param_map["relPIncTauNCtoIO"].value); 
	gIncPCtoBC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncPCtoBC"].value); 
	gIncGRtoBC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncGRtoBC"].value); 
	gIncGRtoSC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["gIncGRtoSC"].value); 
	rawGLeakBC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["rawGLeakBC"].value); 
	rawGLeakGR                 = std::stof(s_file.parsed_var_sections["activity"].param_map["rawGLeakGR"].value); 
	rawGLeakIO                 = std::stof(s_file.parsed_var_sections["activity"].param_map["rawGLeakIO"].value); 
	rawGLeakNC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["rawGLeakNC"].value); 
	rawGLeakPC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["rawGLeakPC"].value); 
	rawGLeakSC                 = std::stof(s_file.parsed_var_sections["activity"].param_map["rawGLeakSC"].value); 
	threshDecTauBC             = std::stof(s_file.parsed_var_sections["activity"].param_map["threshDecTauBC"].value); 
	//threshDecTauGR             = std::stof(s_file.parsed_var_sections["activity"].param_map["threshDecTauGR"].value); 
	threshDecTauIO             = std::stof(s_file.parsed_var_sections["activity"].param_map["threshDecTauIO"].value); 
	threshDecTauNC             = std::stof(s_file.parsed_var_sections["activity"].param_map["threshDecTauNC"].value); 
	threshDecTauPC             = std::stof(s_file.parsed_var_sections["activity"].param_map["threshDecTauPC"].value); 
	threshDecTauSC             = std::stof(s_file.parsed_var_sections["activity"].param_map["threshDecTauSC"].value); 
	threshMaxBC                = std::stof(s_file.parsed_var_sections["activity"].param_map["threshMaxBC"].value); 
	//threshMaxGR                = std::stof(s_file.parsed_var_sections["activity"].param_map["threshMaxGR"].value); 
	threshMaxIO                = std::stof(s_file.parsed_var_sections["activity"].param_map["threshMaxIO"].value); 
	threshMaxNC                = std::stof(s_file.parsed_var_sections["activity"].param_map["threshMaxNC"].value); 
	threshMaxPC                = std::stof(s_file.parsed_var_sections["activity"].param_map["threshMaxPC"].value); 
	threshMaxSC                = std::stof(s_file.parsed_var_sections["activity"].param_map["threshMaxSC"].value); 
	//weightScale                = std::stof(s_file.parsed_var_sections["activity"].param_map["weightScale"].value); 

	/* derived act params */
	// assume that we have initialized conparams already, else will get float error
	//threshDecGR         = 1 - exp(-msPerTimeStep / threshDecTauGR);
	tsPerHistBinGR      = msPerHistBinGR / msPerTimeStep;
	gLeakSC             = rawGLeakSC / (6 - msPerTimeStep);
	gDecGRtoSC          = exp(-msPerTimeStep / gDecTauGRtoSC);
	threshDecSC         = 1 - exp(-msPerTimeStep / threshDecTauSC);
	gDecGRtoBC          = exp(-msPerTimeStep / gDecTauGRtoBC);
	gDecPCtoBC          = exp(-msPerTimeStep / gDecTauPCtoBC);
	threshDecBC         = 1 - exp(-msPerTimeStep / threshDecTauBC);
	threshDecPC         = 1 - exp(-msPerTimeStep / threshDecTauPC);
	gLeakPC             = rawGLeakPC / (6 - msPerTimeStep);
	gDecGRtoPC          = exp(-msPerTimeStep / gDecTauGRtoPC);
	gDecBCtoPC          = exp(-msPerTimeStep / gDecTauBCtoPC);
	gDecSCtoPC          = exp(-msPerTimeStep / gDecTauSCtoPC);
	tsPopHistPC         = 40 / msPerTimeStep;
	tsPerPopHistBinPC   =  5 / msPerTimeStep, 
	//numPopHistBinsPC    =  8.0; tsPopHistPC / tsPerPopHistBinPC
	gLeakIO             = rawGLeakIO / (6 - msPerTimeStep);
	threshDecIO         = 1 - exp(-msPerTimeStep / threshDecTauIO);
	tsLTDDurationIO     = msLTDDurationIO / msPerTimeStep;
	tsLTDstartAPIO      = msLTDStartAPIO  / msPerTimeStep;
	tsLTPstartAPIO      = msLTPStartAPIO  / msPerTimeStep;
	tsLTPEndAPIO        = msLTPEndAPIO    / msPerTimeStep;
	grPCHistCheckBinIO  = abs(msLTPEndAPIO / msPerHistBinGR);
	gDecPCtoNC          = exp(-msPerTimeStep / gDecTauPCtoNC);
	gLeakNC             = rawGLeakNC / (6 - msPerTimeStep);
	threshDecNC         = 1 - exp(-msPerTimeStep / threshDecTauNC);
	gLeakBC             = rawGLeakBC;

	act_params_populated = true;
}

void read_act_params(std::fstream &in_param_buf)
{
	in_param_buf.read((char *)&coupleRiRjRatioIO, sizeof(float));
	in_param_buf.read((char *)&eBCtoPC, sizeof(float));
	in_param_buf.read((char *)&eNCtoIO, sizeof(float)); 
	in_param_buf.read((char *)&ePCtoBC, sizeof(float));
	in_param_buf.read((char *)&ePCtoNC, sizeof(float));
	in_param_buf.read((char *)&eSCtoPC, sizeof(float));
	in_param_buf.read((char *)&gDecTauBCtoPC, sizeof(float));
	in_param_buf.read((char *)&gIncBCtoPC, sizeof(float));
	in_param_buf.read((char *)&gDecT0ofNCtoIO, sizeof(float));
	in_param_buf.read((char *)&gDecTSofNCtoIO, sizeof(float));
	in_param_buf.read((char *)&gDecTTofNCtoIO, sizeof(float));
	in_param_buf.read((char *)&gIncNCtoIO, sizeof(float));
	in_param_buf.read((char *)&gIncTauNCtoIO, sizeof(float));
	in_param_buf.read((char *)&gDecTauPCtoBC, sizeof(float));
	in_param_buf.read((char *)&gDecTauPCtoNC, sizeof(float));
	in_param_buf.read((char *)&gIncAvgPCtoNC, sizeof(float));
	in_param_buf.read((char *)&gDecTauGRtoBC, sizeof(float));
	in_param_buf.read((char *)&gDecTauGRtoPC, sizeof(float));
	in_param_buf.read((char *)&gDecTauGRtoSC, sizeof(float));
	in_param_buf.read((char *)&gIncGRtoPC, sizeof(float));
	in_param_buf.read((char *)&gDecTauSCtoPC, sizeof(float));
	in_param_buf.read((char *)&gIncSCtoPC, sizeof(float));
	in_param_buf.read((char *)&synLTDStepSizeGRtoPC, sizeof(float));
	in_param_buf.read((char *)&synLTPStepSizeGRtoPC, sizeof(float));
	in_param_buf.read((char *)&maxExtIncVIO, sizeof(float));
	in_param_buf.read((char *)&msLTDDurationIO, sizeof(float));
	in_param_buf.read((char *)&msLTDStartAPIO, sizeof(float));
	in_param_buf.read((char *)&msLTPEndAPIO, sizeof(float));
	in_param_buf.read((char *)&msLTPStartAPIO, sizeof(float));
	in_param_buf.read((char *)&msPerHistBinGR, sizeof(float));
	in_param_buf.read((char *)&relPDecT0ofNCtoIO, sizeof(float));
	in_param_buf.read((char *)&relPDecTSofNCtoIO, sizeof(float));
	in_param_buf.read((char *)&relPDecTTofNCtoIO, sizeof(float));
	in_param_buf.read((char *)&relPIncNCtoIO, sizeof(float));
	in_param_buf.read((char *)&relPIncTauNCtoIO, sizeof(float));
	in_param_buf.read((char *)&gIncPCtoBC, sizeof(float));
	in_param_buf.read((char *)&gIncGRtoBC, sizeof(float));
	in_param_buf.read((char *)&gIncGRtoSC, sizeof(float));
	in_param_buf.read((char *)&rawGLeakBC, sizeof(float));
	in_param_buf.read((char *)&rawGLeakGR, sizeof(float));
	in_param_buf.read((char *)&rawGLeakIO, sizeof(float));
	in_param_buf.read((char *)&rawGLeakNC, sizeof(float));
	in_param_buf.read((char *)&rawGLeakPC, sizeof(float));
	in_param_buf.read((char *)&rawGLeakSC, sizeof(float));
	in_param_buf.read((char *)&threshDecTauBC, sizeof(float));
	//in_param_buf.read((char *)&threshDecTauGR, sizeof(float));
	in_param_buf.read((char *)&threshDecTauIO, sizeof(float));
	in_param_buf.read((char *)&threshDecTauNC, sizeof(float));
	in_param_buf.read((char *)&threshDecTauPC, sizeof(float));
	in_param_buf.read((char *)&threshDecTauSC, sizeof(float));
	in_param_buf.read((char *)&threshMaxBC, sizeof(float));
	//in_param_buf.read((char *)&threshMaxGR, sizeof(float));
	in_param_buf.read((char *)&threshMaxIO, sizeof(float));
	in_param_buf.read((char *)&threshMaxNC, sizeof(float));
	in_param_buf.read((char *)&threshMaxPC, sizeof(float));
	in_param_buf.read((char *)&threshMaxSC, sizeof(float));
	//in_param_buf.read((char *)&weightScale, sizeof(float));

	/* derived params */
	//in_param_buf.read((char *)&threshDecGR, sizeof(float));
	in_param_buf.read((char *)&tsPerHistBinGR, sizeof(float));
	in_param_buf.read((char *)&gLeakSC, sizeof(float));
	in_param_buf.read((char *)&gDecGRtoSC, sizeof(float));
	in_param_buf.read((char *)&threshDecSC, sizeof(float));
	in_param_buf.read((char *)&gDecGRtoBC, sizeof(float));
	in_param_buf.read((char *)&gDecPCtoBC, sizeof(float));
	in_param_buf.read((char *)&threshDecBC, sizeof(float));
	in_param_buf.read((char *)&threshDecPC, sizeof(float));
	in_param_buf.read((char *)&gLeakPC, sizeof(float));
	in_param_buf.read((char *)&gDecGRtoPC, sizeof(float));
	in_param_buf.read((char *)&gDecBCtoPC, sizeof(float));
	in_param_buf.read((char *)&gDecSCtoPC, sizeof(float));
	in_param_buf.read((char *)&tsPopHistPC, sizeof(float));
	in_param_buf.read((char *)&tsPerPopHistBinPC, sizeof(float));
	//in_param_buf.read((char *)&numPopHistBinsPC, sizeof(float));
	in_param_buf.read((char *)&gLeakIO, sizeof(float));
	in_param_buf.read((char *)&threshDecIO, sizeof(float));
	in_param_buf.read((char *)&tsLTDDurationIO, sizeof(float));
	in_param_buf.read((char *)&tsLTDstartAPIO, sizeof(float));
	in_param_buf.read((char *)&tsLTPstartAPIO, sizeof(float));
	in_param_buf.read((char *)&tsLTPEndAPIO, sizeof(float)); 
	in_param_buf.read((char *)&grPCHistCheckBinIO, sizeof(float)); 
	in_param_buf.read((char *)&gDecPCtoNC, sizeof(float));
	in_param_buf.read((char *)&gLeakNC, sizeof(float));
	in_param_buf.read((char *)&threshDecNC, sizeof(float));

	act_params_populated = true;
}

void write_act_params(std::fstream &out_param_buf)
{
	out_param_buf.write((char *)&coupleRiRjRatioIO, sizeof(float));
	out_param_buf.write((char *)&eBCtoPC, sizeof(float));
	out_param_buf.write((char *)&eNCtoIO, sizeof(float)); 
	out_param_buf.write((char *)&ePCtoBC, sizeof(float));
	out_param_buf.write((char *)&ePCtoNC, sizeof(float));
	out_param_buf.write((char *)&eSCtoPC, sizeof(float));
	out_param_buf.write((char *)&gDecTauBCtoPC, sizeof(float));
	out_param_buf.write((char *)&gIncBCtoPC, sizeof(float));
	out_param_buf.write((char *)&gDecT0ofNCtoIO, sizeof(float));
	out_param_buf.write((char *)&gDecTSofNCtoIO, sizeof(float));
	out_param_buf.write((char *)&gDecTTofNCtoIO, sizeof(float));
	out_param_buf.write((char *)&gIncNCtoIO, sizeof(float));
	out_param_buf.write((char *)&gIncTauNCtoIO, sizeof(float));
	out_param_buf.write((char *)&gDecTauPCtoBC, sizeof(float));
	out_param_buf.write((char *)&gDecTauPCtoNC, sizeof(float));
	out_param_buf.write((char *)&gIncAvgPCtoNC, sizeof(float));
	out_param_buf.write((char *)&gDecTauGRtoBC, sizeof(float));
	out_param_buf.write((char *)&gDecTauGRtoPC, sizeof(float));
	out_param_buf.write((char *)&gDecTauGRtoSC, sizeof(float));
	out_param_buf.write((char *)&gIncGRtoPC, sizeof(float));
	out_param_buf.write((char *)&gDecTauSCtoPC, sizeof(float));
	out_param_buf.write((char *)&gIncSCtoPC, sizeof(float));
	out_param_buf.write((char *)&synLTDStepSizeGRtoPC, sizeof(float));
	out_param_buf.write((char *)&synLTPStepSizeGRtoPC, sizeof(float));
	out_param_buf.write((char *)&maxExtIncVIO, sizeof(float));
	out_param_buf.write((char *)&msLTDDurationIO, sizeof(float));
	out_param_buf.write((char *)&msLTDStartAPIO, sizeof(float));
	out_param_buf.write((char *)&msLTPEndAPIO, sizeof(float));
	out_param_buf.write((char *)&msLTPStartAPIO, sizeof(float));
	out_param_buf.write((char *)&msPerHistBinGR, sizeof(float));
	out_param_buf.write((char *)&relPDecT0ofNCtoIO, sizeof(float));
	out_param_buf.write((char *)&relPDecTSofNCtoIO, sizeof(float));
	out_param_buf.write((char *)&relPDecTTofNCtoIO, sizeof(float));
	out_param_buf.write((char *)&relPIncNCtoIO, sizeof(float));
	out_param_buf.write((char *)&relPIncTauNCtoIO, sizeof(float));
	out_param_buf.write((char *)&gIncPCtoBC, sizeof(float));
	out_param_buf.write((char *)&gIncGRtoBC, sizeof(float));
	out_param_buf.write((char *)&gIncGRtoSC, sizeof(float));
	out_param_buf.write((char *)&rawGLeakBC, sizeof(float));
	out_param_buf.write((char *)&rawGLeakGR, sizeof(float));
	out_param_buf.write((char *)&rawGLeakIO, sizeof(float));
	out_param_buf.write((char *)&rawGLeakNC, sizeof(float));
	out_param_buf.write((char *)&rawGLeakPC, sizeof(float));
	out_param_buf.write((char *)&rawGLeakSC, sizeof(float));
	out_param_buf.write((char *)&threshDecTauBC, sizeof(float));
	//out_param_buf.write((char *)&threshDecTauGR, sizeof(float));
	out_param_buf.write((char *)&threshDecTauIO, sizeof(float));
	out_param_buf.write((char *)&threshDecTauNC, sizeof(float));
	out_param_buf.write((char *)&threshDecTauPC, sizeof(float));
	out_param_buf.write((char *)&threshDecTauSC, sizeof(float));
	out_param_buf.write((char *)&threshMaxBC, sizeof(float));
	//out_param_buf.write((char *)&threshMaxGR, sizeof(float));
	out_param_buf.write((char *)&threshMaxIO, sizeof(float));
	out_param_buf.write((char *)&threshMaxNC, sizeof(float));
	out_param_buf.write((char *)&threshMaxPC, sizeof(float));
	out_param_buf.write((char *)&threshMaxSC, sizeof(float));
	//out_param_buf.write((char *)&weightScale, sizeof(float));

	/* derived params */
	//out_param_buf.write((char *)&threshDecGR, sizeof(float));
	out_param_buf.write((char *)&tsPerHistBinGR, sizeof(float));
	out_param_buf.write((char *)&gLeakSC, sizeof(float));
	out_param_buf.write((char *)&gDecGRtoSC, sizeof(float));
	out_param_buf.write((char *)&threshDecSC, sizeof(float));
	out_param_buf.write((char *)&gDecGRtoBC, sizeof(float));
	out_param_buf.write((char *)&gDecPCtoBC, sizeof(float));
	out_param_buf.write((char *)&threshDecBC, sizeof(float));
	out_param_buf.write((char *)&threshDecPC, sizeof(float));
	out_param_buf.write((char *)&gLeakPC, sizeof(float));
	out_param_buf.write((char *)&gDecGRtoPC, sizeof(float));
	out_param_buf.write((char *)&gDecBCtoPC, sizeof(float));
	out_param_buf.write((char *)&gDecSCtoPC, sizeof(float));
	out_param_buf.write((char *)&tsPopHistPC, sizeof(float));
	out_param_buf.write((char *)&tsPerPopHistBinPC, sizeof(float));
	//out_param_buf.write((char *)&numPopHistBinsPC, sizeof(float));
	out_param_buf.write((char *)&gLeakIO, sizeof(float));
	out_param_buf.write((char *)&threshDecIO, sizeof(float));
	out_param_buf.write((char *)&tsLTDDurationIO, sizeof(float));
	out_param_buf.write((char *)&tsLTDstartAPIO, sizeof(float));
	out_param_buf.write((char *)&tsLTPstartAPIO, sizeof(float));
	out_param_buf.write((char *)&tsLTPEndAPIO, sizeof(float)); 
	out_param_buf.write((char *)&grPCHistCheckBinIO, sizeof(float)); 
	out_param_buf.write((char *)&gDecPCtoNC, sizeof(float));
	out_param_buf.write((char *)&gLeakNC, sizeof(float));
	out_param_buf.write((char *)&threshDecNC, sizeof(float));

}

