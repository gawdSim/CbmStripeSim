 /*
 * kernels.cu
 *
 *  Created on: Jun 6, 2011
 *      Author: consciousness
 */

#include <curand_kernel.h>
#include "kernels.h"

 extern __shared__ uint32_t sharedIOBufGR[];
 extern __shared__ float  sharedIOBufGRfloat[];

__global__ void testKernel(float *a, float *b, float *c)
 {
 	int i = blockIdx.x * blockDim.x + threadIdx.x;
 	c[i] = a[i] + b[i];
 }

//**---------------random kernels-------------**

template <typename randState>
__global__ void curandSetupKernel(randState *state, uint32_t seed)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	/* every thread gets same seed, different sequence number,
	   no offset */
	curand_init(seed, id, 0, &state[id]);
}

//FIXME: for each ts, rng produces same number
template <typename randState>
__global__ void curandGenerateUniformsKernel(randState *state, float *randoms, uint32_t rand_offset)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	curandStateMRG32k3a localState = state[i];
	randoms[i + rand_offset] = curand_uniform(&localState);
	state[i] = localState;
}

// NOTE: only call this once per thread!!! else curand_uniform
//       will spit out the same numbers!
__device__ float getRandFloat(uint64_t seed, int tid, uint32_t *threadCallCount)
{
	curandStateMRG32k3a d_state;
	curand_init(seed + tid + threadCallCount[tid], 0, 0, &d_state);
	threadCallCount[tid]++;
	return curand_uniform(&d_state);
}

//**---------------end random kernels---------**

//**---------------end random kernels-------------**

//**-----------------GR Kernels------------------**

__global__ void calcActivityGRGPU(float *threshGPU, uint8_t *apGPU, uint32_t *apBufGPU, uint64_t *apHistGPU, float *randoms,
	  float *gr_templateGPU, size_t gr_template_pitchGPU, size_t num_gr_old, size_t ts, float s_per_ts, uint32_t ap_buf_hist_mask, 
	  float threshBase, float threshMax, float threshInc)
{
	int tix = blockIdx.x * blockDim.x + threadIdx.x;
	float *gr_template_row = (float *)((char *)gr_templateGPU + ts * gr_template_pitchGPU);
	//threshGPU[tix] += (threshMax - threshGPU[tix]) * threshInc;
	apGPU[tix] = randoms[tix] < (gr_template_row[tix] * s_per_ts); // * threshGPU[tix]);
	apBufGPU[tix] = (apBufGPU[tix] << 1) | apGPU[tix];
	apHistGPU[tix] <<= 1;
	apHistGPU[tix] |= ((apBufGPU[tix] & ap_buf_hist_mask) > 0) * 0x00000001;

	//threshGPU[tix] = apGPU[tix] * threshBase + (apGPU[tix] - 1) * threshGPU[tix];

}

__global__ void updateGRHistory(uint32_t *apBuf, uint64_t *apHist, uint32_t bufTestMask)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	uint64_t tempHist=apHist[i]<<1;
	apHist[i]=tempHist|((apBuf[i]&bufTestMask)>0)*0x00000001; 
}

__global__ void updatePFBCOutGPU(uint32_t *apBuf, uint32_t *delay,
		uint32_t *pfBC, size_t pfBCPitch, unsigned int numPFInPerBC, unsigned int numPFInPerBCP2)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int *pfBCRow = (uint32_t *)((char *)pfBC + (tid >> numPFInPerBCP2) * pfBCPitch);

	pfBCRow[tid & (numPFInPerBC-1)] = (apBuf[tid] & delay[tid]) > 0;
}


__global__ void updatePFSCOutGPU(uint32_t *apBuf, uint32_t *delay,
		uint32_t *pfSC, size_t pfSCPitch, unsigned int numPFInPerSC, unsigned int numPFInPerSCP2)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int *pfSCRow = (uint32_t *)((char *)pfSC + (tid >> numPFInPerSCP2) * pfSCPitch);

	pfSCRow[tid & (numPFInPerSC-1)] = (apBuf[tid] & delay[tid]) > 0;
}

// keep in mind: grConOut has len numGRPerGPU -> contains pc ids that these gr connect to
// also that now pfPC now has dims (numBlocks, num_pc)
__global__ void updatePFPCOutGPU(uint32_t *apBuf, uint32_t *grConOut, uint32_t *delay,
		float *synWeight, float *pfPC, size_t pfPCPitch, uint32_t nWrites)
{
	int tid=threadIdx.x;
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	uint32_t *pcRow=(uint32_t *)((char *)pfPC+blockIdx.x*pfPCPitch);

	uint32_t tempOut;

	for (uint32_t i=0; i<nWrites; i++)
	{
		sharedIOBufGR[tid+i*blockDim.x]=0;
	}
	__syncthreads();

	tempOut = (apBuf[index] & delay[index]) > 0;
	if (tempOut > 0)
	{
		atomicAdd(&sharedIOBufGR[grConOut[index]], synWeight[index]);
	}
	__syncthreads();

	for(int i=0; i<nWrites; i++)
	{
		pcRow[tid+i*blockDim.x]=sharedIOBufGR[tid+i*blockDim.x];
	}


	//unsigned int tempOut;
	//float *pfPCRow = (float *)((char *)pfPC + (tid >> numPFInPerPCP2) * pfPCPitch);

	//tempOut = (apBuf[tid] & delay[tid]) > 0;

	//pfPCRow[tid & (numPFInPerPC-1)] = synWeight[tid] * tempOut;
}

__global__ void sumPFPCOutGPU(unsigned int nRows, float *pcOut, size_t pcOutPitch, float *pcOutSum)
{
	unsigned int *pcOutRow;
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	unsigned int tempSum;

	tempSum=0;
	for(int i=0; i<nRows; i++)
	{
		pcOutRow=(unsigned int *)((char *)pcOut+i*pcOutPitch);

		tempSum+=pcOutRow[index];
	}
	pcOutSum[index]=tempSum;
}



//**---------------end GR Kernels-------------------**

//**---------------SC Kernels-------------------**

__global__ void calcSCActivity(float *vSC, float *gPFSC, float *threshSC, uint8_t *apSC,
	uint32_t *apBufSC, uint32_t *gInputSumPFSC, float eLeakSC, float gLeakSC, float gIncPFSC,
	float gDecPFSC, float threshRestSC, float threshMaxSC, float threshDecaySC)
{
	int i = blockIdx.x*blockDim.x+threadIdx.x;

	float tempV        = vSC[i];
	float tempThresh   = threshSC[i];
	uint8_t tempAP     = apSC[i];
	uint32_t tempAPBuf = apBufSC[i];

	gPFSC[i] += gInputSumPFSC[i] * gIncPFSC;
	gPFSC[i] *= gDecPFSC;

	tempV += gLeakSC * (eLeakSC - tempV) - gPFSC[i] * tempV;

	tempThresh += threshDecaySC * (threshRestSC - tempThresh);
	tempAP = tempV > tempThresh;

	tempAPBuf = (tempAPBuf << 1) | (tempAP * 0x00000001);
	tempThresh = tempAP * threshMaxSC + (1 - tempAP) * tempThresh;

	vSC[i]      = tempV;
	threshSC[i]  = tempThresh;
	apSC[i]      = tempAP;
	apBufSC[i]   = tempAPBuf;
}
//**---------------end SC Kernels-------------------**

//**---------------PC kernels-----------------**

__global__ void calcPCActivity(float *vPC, float *gPFPC, float *gSCPC, float *gBCPC, float *threshPC,
	uint8_t *apPC, uint32_t *apBufPC, float *gInputSumPFPC, float *gInputSumSCPC, float *gInputSumBCPC,
	float eLeakPC, float gLeakPC, float gIncPFPC, float gIncSCPC, float gIncBCPC, float eSCtoPC, float eBCtoPC,
	float gDecPFPC, float gDecSCPC, float gDecBCPC, float threshRestPC, float threshMaxPC, float threshDecayPC)
{
	int i = blockIdx.x*blockDim.x+threadIdx.x;

	float tempV        = vPC[i];
	float tempThresh   = threshPC[i];
	uint8_t tempAP     = apPC[i];
	uint32_t tempAPBuf = apBufPC[i];

	gPFPC[i] += gInputSumPFPC[i] * gIncPFPC;
	gPFPC[i] *= gDecPFPC;

	gSCPC[i] += gInputSumSCPC[i] * gIncSCPC;
	gSCPC[i] *= gDecSCPC;

	gBCPC[i] += gInputSumBCPC[i] * gIncBCPC;
	gBCPC[i] *= gDecBCPC;

	tempV    += gLeakPC * (eLeakPC - tempV)
			  - gPFPC[i] * tempV 
			  + gSCPC[i] * (eSCtoPC - tempV)
			  + gBCPC[i] * (eBCtoPC - tempV);

	tempThresh +=  threshDecayPC * (threshRestPC - tempThresh);
	tempAP = tempV > tempThresh;
	
	tempAPBuf = (tempAPBuf << 1) | (tempAP * 0x00000001);
	tempThresh = tempAP * threshMaxPC + (1 - tempAP) * tempThresh;

	vPC[i]      = tempV;
	threshPC[i]  = tempThresh;
	apPC[i]      = tempAP;
	apBufPC[i]   = tempAPBuf;
}

//**---------------end PC kernels-------------**

//**---------------BC kernels-----------------**

__global__ void calcBCActivity(float *vBC, float *gPFBC, float *gPCBC, float *threshBC,
	uint8_t *apBC, uint32_t *apBufBC, uint32_t *gInputSumPFBC, uint32_t *gInputSumPCBC, 
	float eLeakBC, float gLeakBC, float gIncPFBC, float gIncPCBC, float ePCtoBC,
	float gDecPFBC, float gDecPCBC, float threshRestBC, float threshMaxBC, float threshDecayBC)
{
	int i = blockIdx.x*blockDim.x+threadIdx.x;

	float tempV        = vBC[i];
	float tempThresh   = threshBC[i];
	uint8_t tempAP     = apBC[i];
	uint32_t tempAPBuf = apBufBC[i];

	gPFBC[i] += gInputSumPFBC[i] * gIncPFBC;
	gPFBC[i] *= gDecPFBC;

	gPCBC[i] += gInputSumPCBC[i] * gIncPCBC;
	gPCBC[i] *= gDecPCBC;

	tempV    += gLeakBC * (eLeakBC - tempV)
			  - gPFBC[i] * tempV
			  + gPCBC[i] * (ePCtoBC - tempV);

	tempThresh +=  threshDecayBC * (threshRestBC - tempThresh);
	tempAP = tempV > tempThresh;
	
	tempAPBuf = (tempAPBuf << 1) | (tempAP * 0x00000001);
	tempThresh = tempAP * threshMaxBC + (1 - tempAP) * tempThresh;

	vBC[i]      = tempV;
	threshBC[i]  = tempThresh;
	apBC[i]      = tempAP;
	apBufBC[i]   = tempAPBuf;
}

//**---------------end BC kernels-------------**

//**---------------IO kernels-----------------**

__global__ void updatePFPCSynIO(float *synWPFPC, uint64_t *historyGR, uint64_t plastCheckMask,
		unsigned int offset, float plastStep)
{

	int i=blockIdx.x*blockDim.x+threadIdx.x+offset;
	synWPFPC[i]=synWPFPC[i]+((historyGR[i]&plastCheckMask)>0)*plastStep;

	synWPFPC[i]=(synWPFPC[i]>0)*synWPFPC[i];
	synWPFPC[i]=(synWPFPC[i]>1)+(synWPFPC[i]<=1)*synWPFPC[i];
}

// TODO: compute the gInputSumNCIO separately
// NOTE: 'threadCallCounts' is an array of size num_io
__global__ void calcIOActivity(uint64_t seed, uint32_t *threadCallCounts, float *vIO,
	float *threshIO, uint8_t *apIO, uint32_t *apBufIO, float *gInputSumNCIO,
	float *vCoupleIO, float eLeakIO, float gLeakIO, float eNCtoIO, float *errDrive,
	float threshRestIO, float threshMaxIO, float threshDecayIO)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	float gNoise = getRandFloat(seed, tid, threadCallCounts);
	gNoise = 2.0 * (gNoise - 0.5);

	float tempGInputSum = 1.5 * gInputSumNCIO[tid] / 3.1;
	float tempV         = vIO[tid];
	float tempThresh    = threshIO[tid];
	uint8_t tempAP      = apIO[tid];
	uint32_t tempAPBuf  = apBufIO[tid];

	tempV  += gLeakIO * (eLeakIO - tempV)
			+ tempGInputSum * (eNCtoIO - tempV)
			+ vCoupleIO[tid]
			+ *errDrive + gNoise;

	tempThresh +=  threshDecayIO * (threshRestIO - tempThresh);
	tempAP = tempV > tempThresh;

	tempAPBuf = (tempAPBuf << 1) | (tempAP * 0x00000001);
	tempThresh = tempAP * threshMaxIO + (1 - tempAP) * tempThresh;

	vIO[tid]     = tempV;
	threshIO[tid] = tempThresh;
	apIO[tid]     = tempAP;
	apBufIO[tid]  = tempAPBuf;

	*errDrive = 0.0; // for now we are going to assume that we want some global error drive...
}

//**---------------end IO kernels-------------**

//**---------------NC kernels-----------------**

// TODO: precompute the input sums elsewhere: update the conductance delays here
__global__ void calcNCActivity(float *vNC, float *threshNC,
	uint8_t *apNC, uint32_t *apBufNC, float *gInputSumPCNC, float *gInputSumMFNCNMDA,
	float *gInputSumMFNCAMPA, float eLeakNC, float gLeakNC, float ePCtoNC,
	float threshRestNC, float threshMaxNC, float threshDecayNC)
{
	int i = blockIdx.x*blockDim.x+threadIdx.x;

	float tempV         = vNC[i];
	float tempThresh    = threshNC[i];
	uint8_t tempAP      = apNC[i];
	uint32_t tempAPBuf  = apBufNC[i];

	tempV  += gLeakNC * (eLeakNC - tempV)
			- (gInputSumMFNCNMDA[i] + gInputSumMFNCAMPA[i]) * tempV
			+ gInputSumPCNC[i] * (ePCtoNC - tempV);

	tempThresh +=  threshDecayNC * (threshRestNC - tempThresh);
	tempAP = tempV > tempThresh;

	tempAPBuf = (tempAPBuf << 1) | (tempAP * 0x00000001);
	tempThresh = tempAP * threshMaxNC + (1 - tempAP) * tempThresh;

	vNC[i] = tempV;
	threshNC[i] = tempThresh;
	apNC[i] = tempAP;
	apBufNC[i] = tempAPBuf;
}

//**---------------end NC kernels-------------**


//**---------------common kernels-------------**

template <typename Type, unsigned int blockSize, unsigned int sDataSize, bool inMultiPitch, bool outMultiPitch>
__global__ void sumInputsNew(Type *input, unsigned int inputPitch,
		Type *output, unsigned int outputPitch, unsigned int rowLength)
{
		__shared__ Type sData[sDataSize];

	int tid=threadIdx.x;
	int index=blockIdx.x*(blockSize*2)+tid;
	int gridSize=blockSize*2*gridDim.x;
	Type *inputRow;

	Type tempSum=0;

	if(inMultiPitch)
	{
		inputRow=(Type *)((char *)input+blockIdx.y*inputPitch);
	}
	else
	{
		inputRow=input+blockIdx.y;
	}

	while(index<rowLength)
	{
		tempSum+=inputRow[index]+inputRow[index+blockSize];
		index+=gridSize;
	}
	sData[tid]=tempSum;
	__syncthreads();

	if(blockSize>=512)
	{
		if(tid<256)
			sData[tid]+=sData[tid+256];
		__syncthreads();
	}

	if(blockSize>=256)
	{
		if(tid<128)
			sData[tid]+=sData[tid+128];
		__syncthreads();
	}

	if(blockSize>=128)
	{
		if(tid<64)
			sData[tid]+=sData[tid+64];
		__syncthreads();
	}

	if(tid<32)
	{
		volatile Type* sMem = sData;
		if(blockSize>=64)
			sMem[tid]+=sMem[tid+32];
		if(blockSize>=32)
			sMem[tid]+=sMem[tid+16];
		if(blockSize>=16)
			sMem[tid]+=sMem[tid+8];
		if(blockSize>=8)
			sMem[tid]+=sMem[tid+4];
		if(blockSize>=4)
			sMem[tid]+=sMem[tid+2];
		if(blockSize>=2)
			sMem[tid]+=sMem[tid+1];
	}
	if(tid==0)
	{
		Type *outputRow;
		if(outMultiPitch)
		{
			outputRow=(Type *)((char *)output+blockIdx.y*outputPitch);
		}
		else
		{
			outputRow=output+blockIdx.y;
		}
		outputRow[blockIdx.x]=sData[0];
	}
}

template<typename Type>
__global__ void broadcastValue(Type *val, Type *outArr)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;

	outArr[i]=*val;
}

//**---------------end common kernels---------**


//**---------------kernel calls---------------**

void callTestKernel(cudaStream_t &st, float *a, float *b, float *c)
{
	testKernel<<<1, 128>>>(a, b, c);
}

template <typename randState, typename blockDims, typename threadDims>
void callCurandSetupKernel(cudaStream_t &st, randState *state, uint32_t seed,
						   blockDims &block_dim, threadDims &thread_dim)
{
	curandSetupKernel<randState><<<block_dim, thread_dim, 0, st>>>(state, seed);
}

template <typename randState>
void callCurandGenerateUniformKernel(cudaStream_t &st, randState *state, uint32_t block_dim,
	  uint32_t thread_dim, float *randoms, size_t rand_offset)
{
	curandGenerateUniformsKernel<randState><<<block_dim, thread_dim, 0, st>>>(state, randoms, rand_offset);
}

void callGRActKernel(cudaStream_t &st, uint32_t numBlocks, uint32_t numGRPerBlock,
    float *threshGPU, uint8_t *apGPU, uint32_t *apBufGPU, uint64_t *apHistGPU, float *randoms, float *gr_templateGPU,
	size_t gr_template_pitchGPU, size_t num_gr_old, size_t ts, float s_per_ts, uint32_t ap_buf_hist_mask, float threshBase,
	float threshMax, float threshInc)
{
	calcActivityGRGPU<<<numBlocks, numGRPerBlock, 0, st>>>(threshGPU, apGPU, apBufGPU, apHistGPU, randoms, gr_templateGPU,
      gr_template_pitchGPU, num_gr_old, ts, s_per_ts, ap_buf_hist_mask, threshBase, threshMax, threshInc);
}

void callPCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numPCPerBlock,
	float *vPC, float *gPFPC, float *gSCPC, float *gBCPC, float *threshPC, uint8_t *apPC, uint32_t *apBufPC,
	float *gInputSumPFPC, float *gInputSumSCPC, float *gInputSumBCPC, float eLeakPC, float gLeakPC, float gIncPFPC,
	float gIncSCPC, float gIncBCPC, float eSCtoPC, float eBCtoPC, float gDecPFPC, float gDecSCPC, float gDecBCPC,
	float threshRestPC, float threshMaxPC, float threshDecayPC)
{
	calcPCActivity<<<numBlocks, numPCPerBlock, 0, st>>>(vPC, gPFPC, gSCPC, gBCPC, threshPC, apPC, apBufPC,
		gInputSumPFPC, gInputSumSCPC, gInputSumBCPC, eLeakPC, gLeakPC, gIncPFPC, gIncSCPC, gIncBCPC,
		eSCtoPC, eBCtoPC, gDecPFPC, gDecSCPC, gDecBCPC, threshRestPC, threshMaxPC, threshDecayPC);
}

void callSCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numSCPerBlock, 
	float *vSC, float *gPFSC, float *threshSC, uint8_t *apSC, uint32_t *apBufSC, uint32_t *gInputSumPFSC,
	float eLeakSC, float gLeakSC, float gIncPFSC, float gDecPFSC, float threshRestSC, float threshMaxSC,
	float threshDecaySC)
{
	calcSCActivity<<<numBlocks, numSCPerBlock, 0, st>>>(vSC, gPFSC, threshSC, apSC, apBufSC, gInputSumPFSC,
		eLeakSC, gLeakSC, gIncPFSC, gDecPFSC, threshRestSC, threshMaxSC, threshDecaySC);
}

void callBCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numBCPerBlock, 
	float *vBC, float *gPFBC, float *gPCBC, float *threshBC, uint8_t *apBC, uint32_t *apBufBC,
	uint32_t *gInputSumPFBC, uint32_t *gInputSumPCBC, float eLeakBC, float gLeakBC, float gIncPFBC,
	float gIncPCBC, float ePCtoBC, float gDecPFBC, float gDecPCBC, float threshRestBC,
	float threshMaxBC, float threshDecayBC)
{
	calcBCActivity<<<numBlocks, numBCPerBlock, 0, st>>>(vBC, gPFBC, gPCBC, threshBC, apBC, apBufBC,
		gInputSumPFBC, gInputSumPCBC, eLeakBC, gLeakBC, gIncPFBC, gIncPCBC, ePCtoBC, gDecPFBC,
		gDecPCBC, threshRestBC, threshMaxBC, threshDecayBC);
}

void callNCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numNCPerBlock,
	float *vNC, float *threshNC, uint8_t *apNC, uint32_t *apBufNC, float *gInputSumPCNC,
	float *gInputSumMFNCNMDA, float *gInputSumMFNCAMPA, float eLeakNC, float gLeakNC,
	float ePCtoNC, float threshRestNC, float threshMaxNC, float threshDecayNC)
{
	calcNCActivity<<<numBlocks, numNCPerBlock, 0, st>>>(vNC, threshNC, apNC, apBufNC,
		gInputSumPCNC, gInputSumMFNCNMDA, gInputSumMFNCAMPA, eLeakNC, gLeakNC, ePCtoNC,
		threshRestNC, threshMaxNC, threshDecayNC);
}

void callIOActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numIOPerBlock, 
	uint64_t seed, uint32_t *threadCallCounts, float *vIO, float *threshIO, uint8_t *apIO,
	uint32_t *apBufIO, float *gInputSumNCIO, float *vCoupleIO, float eLeakIO, float gLeakIO,
	float eNCtoIO, float *errDrive, float threshRestIO, float threshMaxIO, float threshDecayIO)
{
	calcIOActivity<<<numBlocks, numIOPerBlock, 0, st>>>(seed, threadCallCounts, vIO, threshIO,
		apIO, apBufIO, gInputSumNCIO, vCoupleIO, eLeakIO, gLeakIO, eNCtoIO, errDrive,
		threshRestIO, threshMaxIO, threshDecayIO);
}

template<typename Type, bool inMultiP, bool outMultiP>
void callSumKernel(cudaStream_t &st, Type *inGPU, size_t inGPUP, Type *outSumGPU, size_t outSumGPUP,
		unsigned int nOutCells, unsigned int nOutCols, unsigned int rowLength)
{
	unsigned int numElementsPerBlock;
	dim3 dimGrid(nOutCols, nOutCells);

	numElementsPerBlock=rowLength/nOutCols;

	if(numElementsPerBlock>=2048)
	{
		sumInputsNew<Type, 512, 512, inMultiP, outMultiP><<<dimGrid, 512, 0, st>>>
				(inGPU, inGPUP, outSumGPU, outSumGPUP, rowLength);
	}
	else if(numElementsPerBlock>=512)
	{
		sumInputsNew<Type, 128, 128, inMultiP, outMultiP><<<dimGrid, 128, 0, st>>>
				(inGPU, inGPUP, outSumGPU, outSumGPUP, rowLength);
	}
	else if(numElementsPerBlock>=128)
	{
		sumInputsNew<Type, 32, 64, inMultiP, outMultiP><<<dimGrid, 32, 0, st>>>
				(inGPU, inGPUP, outSumGPU, outSumGPUP, rowLength);
	}
	else if(numElementsPerBlock>=32)
	{
		sumInputsNew<Type, 8, 64, inMultiP, outMultiP><<<dimGrid, 8, 0, st>>>
				(inGPU, inGPUP, outSumGPU, outSumGPUP, rowLength);
	}
	else
	{
		sumInputsNew<Type, 2, 64, inMultiP, outMultiP><<<dimGrid, 2, 0, st>>>
				(inGPU, inGPUP, outSumGPU, outSumGPUP, rowLength);
	}
}

template<typename Type>
void callBroadcastKernel(cudaStream_t &st, Type *broadcastVal, Type *outArray,
		unsigned int nBlocks, unsigned int rowLength)
{
	broadcastValue<Type><<<nBlocks, rowLength/nBlocks, 0, st>>>(broadcastVal, outArray);
}

void callUpdatePFBCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint32_t *delayMaskGPU, uint32_t *inPFBCGPU, size_t inPFBCGPUPitch,
		unsigned int numPFInPerBCP2)
{
	updatePFBCOutGPU<<<numBlocks, numGRPerBlock, 0, st>>>(apBufGPU, delayMaskGPU,
			inPFBCGPU, inPFBCGPUPitch, 1<<numPFInPerBCP2, numPFInPerBCP2);
}

void callUpdatePFSCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint32_t *delayMaskGPU, uint32_t *inPFSCGPU, size_t inPFSCGPUPitch,
		unsigned int numPFInPerSCP2)
{
	updatePFSCOutGPU<<<numBlocks, numGRPerBlock, 0, st>>>(apBufGPU, delayMaskGPU,
			inPFSCGPU, inPFSCGPUPitch, 1<<numPFInPerSCP2, numPFInPerSCP2);
}

void callUpdatePFPCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t numPC, uint32_t *apBufGPU, uint32_t *grOutConGPU, uint32_t *delayMaskGPU, float *pfPCSynWGPU,
		float *inPFPCGPU, size_t inPFPCGPUPitch)
{
	updatePFPCOutGPU<<<numBlocks, numGRPerBlock, numPC * sizeof(uint32_t), st>>>(apBufGPU, grOutConGPU, delayMaskGPU,
		pfPCSynWGPU, inPFPCGPU, inPFPCGPUPitch, numPC / numGRPerBlock);
}

void callSumPFPCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numPCPerBlock,
		unsigned int numPFPCOutRows, float *grInPCGPU,  size_t pfInPCGPUPitch, float *pfInPCSGPU)
{
	sumPFPCOutGPU<<<numBlocks, numPCPerBlock, 0, st>>>(numPFPCOutRows, grInPCGPU, pfInPCGPUPitch, pfInPCSGPU);
}

void callUpdateGRHistKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint64_t *historyGPU, uint32_t apBufGRHistMask)
{
		updateGRHistory<<<numBlocks, numGRPerBlock, 0, st>>>(apBufGPU, historyGPU, apBufGRHistMask);
}

void callUpdatePFPCPlasticityIOKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		float *synWeightGPU, uint64_t *historyGPU, unsigned int pastBinNToCheck,
		int offSet, float pfPCPlastStep)
{
	uint64_t mask = ((uint64_t)1)<<(pastBinNToCheck-1);
		updatePFPCSynIO<<<numBlocks, numGRPerBlock, 0, st>>>(synWeightGPU, historyGPU,
				mask, offSet, pfPCPlastStep);
}

//**---------------end kernel calls------------**

// template initializations

template void callCurandSetupKernel<curandStateMRG32k3a, dim3, dim3>
(cudaStream_t &st, curandStateMRG32k3a *state, uint32_t seed, dim3 &block_dim, dim3 &thread_dim);

template void callCurandGenerateUniformKernel<curandStateMRG32k3a>(cudaStream_t &st, curandStateMRG32k3a *state,
		uint32_t block_dim, uint32_t thread_dim, float *randoms, size_t rand_offset);

template void callSumKernel<float, true, false>
		(cudaStream_t &st, float *inPFGPU, size_t inPFGPUP, float *outPFSumGPU, size_t outPFSumGPUP,
		unsigned int nOutCells, unsigned int nOutCols, unsigned int rowLength);

template void callSumKernel<uint32_t, true, false>
		(cudaStream_t &st, uint32_t *inPFGPU, size_t inPFGPUP, uint32_t *outPFSumGPU, size_t outPFSumGPUP,
		unsigned int nOutCells, unsigned int nOutCols, unsigned int rowLength);

template void callSumKernel<uint32_t, false, false>
		(cudaStream_t &st, uint32_t *inPFGPU, size_t inPFGPUP, uint32_t *outPFSumGPU, size_t outPFSumGPUP,
		unsigned int nOutCells, unsigned int nOutCols, unsigned int rowLength);

template void callBroadcastKernel<uint32_t>
(cudaStream_t &st, uint32_t *broadCastVal, uint32_t *outArray, unsigned int nBlocks, unsigned int rowLength);

