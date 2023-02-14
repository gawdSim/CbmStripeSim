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

// NOTE: only call this once per thread!!! else curand_uniform
//       will spit out the same numbers!
__device__ float getRandFloat(uint64_t seed, int tid, uint32_t *threadCallCount)
{
	curandStateMRG32k3a d_state;
	curand_init(seed + tid + threadCallCount[tid], 0, 0, &d_state);
	threadCallCount[tid]++;
	return curand_uniform(&d_state);
}

//**---------------end random kernels-------------**

//**-----------------GR Kernels------------------**

__global__ void calcActivityGRGPU(float *vm, float *gKCa, float *gLeak, float *gNMDA, float *gNMDAInc,
	float *thresh, uint32_t *apBuf, uint8_t *apOutGR, uint32_t *apGR, int *apMFtoGR,
	float *gESum, float *gISum, float eLeak, float eGOIn, float gAMPAInc, 
	float threshBase, float threshMax, float threshDecay)
{
	float tempThresh;
	unsigned int tempAP;

	int i = blockIdx.x * blockDim.x + threadIdx.x;

	float tempGKCa = gKCa[i];
	float tempV = vm[i];

	gLeak[i] = 0.0000001021370733 * tempV * tempV * tempV * tempV
	   		 + 0.00001636462 * tempV * tempV * tempV
			 + 0.00113971219 * tempV * tempV
			 + 0.038772 * tempV
			 + 0.6234929;
	
	
	gNMDAInc[i] = 0.00000011969 * tempV * tempV * tempV
	   			+ 0.000089369 * tempV * tempV
				+ 0.0151 * tempV
				+ 0.7713;

	gNMDA[i] = gNMDAInc[i] * gAMPAInc * apMFtoGR[i] + gNMDA[i] * 0.9672;


	tempV = tempV + gLeak[i] * (eLeak - tempV) - gESum[i] * tempV 
	   	  - gNMDA[i] * tempV + gISum[i] * (eGOIn - tempV); 

	if (tempV > threshMax) tempV = threshMax;

	tempThresh = thresh[i] + (threshBase - thresh[i]) * threshDecay;
	tempAP 	   = tempV > tempThresh;
	thresh[i]  = tempAP * threshMax + (!tempAP) * tempThresh;

	tempGKCa = tempGKCa * 0.9999f; 
	gKCa[i] = tempAP * (tempGKCa + 0.000f) + (!tempAP) * tempGKCa;

	apBuf[i]   = (apBuf[i] << 1) | tempAP;
	apOutGR[i] = tempAP;
	apGR[i]    = tempAP;
	vm[i]      = tempV;
}

__global__ void updateGRBCOutGPU(uint32_t *apBuf,
		uint32_t *bcOut, size_t bcOutPitch,
		uint32_t *delay, size_t delayPitch,
		uint32_t *con, size_t conPitch,
		int32_t *numSyn, int nWrites)
{
	int tid=threadIdx.x;
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	unsigned int *conRow;
	unsigned int *delayRow;
	unsigned int *bcRow=(unsigned int *)((char *)bcOut+blockIdx.x*bcOutPitch);

	int tempNS=numSyn[index];
	unsigned int tempOut;

	for(int i=0; i<nWrites; i++)
	{
		sharedIOBufGR[tid+i*blockDim.x]=0;
	}

	__syncthreads();
	for(int i=0; i<tempNS; i++)
	{
		conRow=(uint32_t *)((char *)con+i*conPitch);
		delayRow=(uint32_t *)((char *)delay+i*delayPitch);

		tempOut=(apBuf[index]&delayRow[index])>0;

		if(tempOut>0)
		{
			atomicAdd(&sharedIOBufGR[conRow[index]], 1);
		}
	}
	__syncthreads();
	for(int i=0; i<nWrites; i++)
	{
		bcRow[tid+i*blockDim.x]=sharedIOBufGR[tid+i*blockDim.x];
	}
}

__global__ void sumGRBCOutGPU(unsigned int nRows, uint32_t *bcOut, size_t bcOutPitch, uint32_t *bcOutSum)
{
	unsigned int *bcOutRow;
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	unsigned int tempSum;

	tempSum=0;
	for(int i=0; i<nRows; i++)
	{
		bcOutRow=(unsigned int *)((char *)bcOut+i*bcOutPitch);

		tempSum+=bcOutRow[index];
	}

	bcOutSum[index]=tempSum;
	//bcOutSum[index]=1;
}

__global__ void updateGRHistory(uint32_t *apBuf, uint64_t *apHist, uint32_t bufTestMask)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	uint64_t tempHist=apHist[i]<<1;
	apHist[i]=tempHist|((apBuf[i]&bufTestMask)>0)*0x00000001; 
}

__global__ void updatePFBCSCOutGPU(uint32_t *apBuf, uint32_t *delay,
		uint32_t *pfBC, size_t pfBCPitch, unsigned int numPFInPerBC, unsigned int numPFInPerBCP2,
		uint32_t *pfSC, size_t pfSCPitch, unsigned int numPFInPerSC, unsigned int numPFInPerSCP2)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	uint32_t tempOut;
	unsigned int *pfBCRow=(uint32_t *)((char *)pfBC+(index>>numPFInPerBCP2)*pfBCPitch);
	unsigned int *pfSCRow=(uint32_t *)((char *)pfSC+(index>>numPFInPerSCP2)*pfSCPitch);

	tempOut=(apBuf[index]&delay[index])>0;

	pfBCRow[index&(numPFInPerBC-1)]=tempOut;
	pfSCRow[index&(numPFInPerSC-1)]=tempOut;
}

__global__ void updatePFPCOutGPU(uint32_t *apBuf, uint32_t *delay,
		float *synWeight, float *pfPC, size_t pfPCPitch, unsigned int numPFInPerPC, unsigned int numPFInPerPCP2)
{
	int index=blockIdx.x*blockDim.x+threadIdx.x;
	unsigned int tempOut;
	float *pfPCRow=(float *)((char *)pfPC+(index>>numPFInPerPCP2)*pfPCPitch);

	tempOut=(apBuf[index]&delay[index])>0;

	pfPCRow[index&(numPFInPerPC-1)]=synWeight[index]*tempOut;
}

//**---------------end GR Kernels-------------------**

//**---------------SC Kernels-------------------**

__global__ void calcSCActivity(float *vSC, float *gPFSC, float *threshSC, uint8_t *apSC,
	uint32_t *apBufSC,  float *gInputSumPFSC, float eLeakSC, float gLeakSC, float gIncPFSC,
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
	uint8_t *apBC, uint32_t *apBufBC, float *gInputSumPFBC, float *gInputSumPCBC, 
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

// TODO: update as a gamma generator
void callGRActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		float *vGPU, float *gKCaGPU, float *gLeakGPU, float *gNMDAGRGPU, float *gNMDAIncGRGPU,
		float *threshGPU, uint32_t *apBufGPU, uint8_t *apOutGRGPU, uint32_t *apGRGPU,
		int *apMFtoGRGPU, float *gESumGPU, float *gISumGPU, float eLeak,
		float eGOIn, float gAMPAInc, float threshBase, float threshMax, float threshDecay)
{
	calcActivityGRGPU<<<numBlocks, numGRPerBlock, 0, st>>>(vGPU, gKCaGPU, gLeakGPU, gNMDAGRGPU,
		  gNMDAIncGRGPU, threshGPU, apBufGPU, apOutGRGPU, apGRGPU, apMFtoGRGPU, gESumGPU, 
		  gISumGPU, eLeak, eGOIn, gAMPAInc, threshBase, threshMax, threshDecay);
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
	float *vSC, float *gPFSC, float *threshSC, uint8_t *apSC, uint32_t *apBufSC, float *gInputSumPFSC,
	float eLeakSC, float gLeakSC, float gIncPFSC, float gDecPFSC, float threshRestSC, float threshMaxSC,
	float threshDecaySC)
{
	calcSCActivity<<<numBlocks, numSCPerBlock, 0, st>>>(vSC, gPFSC, threshSC, apSC, apBufSC, gInputSumPFSC,
		eLeakSC, gLeakSC, gIncPFSC, gDecPFSC, threshRestSC, threshMaxSC, threshDecaySC);
}

void callBCActKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numBCPerBlock, 
	float *vBC, float *gPFBC, float *gPCBC, float *threshBC, uint8_t *apBC, uint32_t *apBufBC,
	float *gInputSumPFBC, float *gInputSumPCBC, float eLeakBC, float gLeakBC, float gIncPFBC,
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

void callSumGRBCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numBCPerBlock,
		unsigned int numGROutRows, uint32_t *grInBCGPU,  size_t grInBCGPUPitch, uint32_t *grInBCSGPU)
{
	sumGRBCOutGPU<<<numBlocks, numBCPerBlock, 0, st>>>(numGROutRows, grInBCGPU, grInBCGPUPitch, grInBCSGPU);
}

void callUpdatePFBCSCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint32_t *delayMaskGPU,
		uint32_t *inPFBCGPU, size_t inPFBCGPUPitch, unsigned int numPFInPerBCP2,
		uint32_t *inPFSCGPU, size_t inPFSCGPUPitch, unsigned int numPFInPerSCP2)
{
	updatePFBCSCOutGPU<<<numBlocks, numGRPerBlock, 0, st>>>(apBufGPU, delayMaskGPU,
			inPFBCGPU, inPFBCGPUPitch, 1<<numPFInPerBCP2, numPFInPerBCP2,
			inPFSCGPU, inPFSCGPUPitch, 1<<numPFInPerSCP2, numPFInPerSCP2);
}

void callUpdatePFPCOutKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock,
		uint32_t *apBufGPU, uint32_t *delayMaskGPU, float *pfPCSynWGPU, float *inPFPCGPU,
		size_t inPFPCGPUPitch, unsigned int numPFInPerPCP2)
{
	updatePFPCOutGPU<<<numBlocks, numGRPerBlock, 0, st>>>(apBufGPU, delayMaskGPU, pfPCSynWGPU,
			inPFPCGPU, inPFPCGPUPitch, 1<<numPFInPerPCP2, numPFInPerPCP2);
}

void callUpdateGROutBCKernel(cudaStream_t &st, unsigned int numBlocks, unsigned int numGRPerBlock, unsigned int numBC,
		uint32_t *apBufGPU, uint32_t *grInBCGPU, uint32_t grInBCGPUPitch, uint32_t *delayMasksGPU, uint32_t delayMasksGPUPitch,
		uint32_t *conGRtoBCGPU, size_t conGRtoBCGPUPitch, int32_t *numBCPerGRGPU)
{
	updateGRBCOutGPU<<<numBlocks, numGRPerBlock, numBC*sizeof(uint32_t), st>>>(apBufGPU, grInBCGPU, grInBCGPUPitch,
			delayMasksGPU, delayMasksGPUPitch, conGRtoBCGPU, conGRtoBCGPUPitch, numBCPerGRGPU, numBC/numGRPerBlock);
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

