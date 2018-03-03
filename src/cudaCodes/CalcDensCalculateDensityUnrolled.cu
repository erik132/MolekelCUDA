#include "CalcDensCalculateDensityUnrolled.cuh"

#include "molekelHelpFunctions/CalcChiCalculateDensityUnrolled.cu"

#include "gputimer.h"


#if __CUDA_ARCH__ < 600 
__device__ double atomicAddLegacy(double* address, double val) { 
	unsigned long long int* address_as_ull = (unsigned long long int*)address; 
	unsigned long long int old = *address_as_ull, assumed; 
	do { 
		assumed = old; 
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed))); 
		// Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN) 
	} while (assumed != old); 
	return __longlong_as_double(old); 
} 
#endif


__global__ void calculateDensityUnrolled(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, float *densities, double *results){
	double result = 0;
	int indexZ = threadIdx.z + (blockDim.z*blockIdx.z);
	int	indexY = threadIdx.y + (blockDim.y*blockIdx.y);
	int	indexX = threadIdx.x + (blockDim.x*(blockIdx.x + internalData.offsetx));
	int realX, rowNr;
	float x,y,z;

	realX = indexX / molecule->nBasisFunctions;
	rowNr = indexX - (molecule->nBasisFunctions * realX);

	if(realX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
			
		x = internalData.dim0 + realX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim4 + indexZ*internalData.dz;
		
		result = calcChiCalculateDensityUnrolled(densities,orbital, molecule, x, y, z, internalData.densityLength,rowNr);
		atomicAddLegacy(&results[realX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)], result);
		//results[realX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] += result;
	}
}


vtkImageData* CalcDensCalculateDensityUnrolled::runComputation(){
	char buffer[100];
	cudaError_t status;
	int i, j, k, counter=0;
	vtkImageData* imageData;
	int originalx =0;

	dim3 blockSize = this->getBlockSize();
	dim3 gridSize = this->getGridSize(blockSize);

	sprintf(buffer, "original grid size %d %d %d block size %d %d %d", gridSize.x,gridSize.y,gridSize.z, blockSize.x,blockSize.y,blockSize.z);
	this->esl->logMessage(buffer);

	originalx = gridSize.x;
	gridSize = this->limitBlocks(10000,gridSize);


	sprintf(buffer, "grid size %d %d %d block size %d %d %d", gridSize.x,gridSize.y,gridSize.z, blockSize.x,blockSize.y,blockSize.z);
	this->esl->logMessage(buffer);
	
	for(this->calcData.offsetx = 0; this->calcData.offsetx < originalx; this->calcData.offsetx += gridSize.x){
		calculateDensityUnrolled<<<gridSize, blockSize>>>(this->deviceMolecule,this->calcData,this->deviceOrbital,this->deviceDensityMatrix, this->deviceResults);

		status = cudaGetLastError();
		if(status != cudaSuccess){
			sprintf(buffer, "kernel launch failed, errorcode %s", cudaGetErrorString(status));
			this->esl->logMessage(buffer);
			return NULL;
		}

		status = cudaDeviceSynchronize();
		if(status != cudaSuccess){
			sprintf(buffer, "failed to sync devices %s", cudaGetErrorString(status));
			this->esl->logMessage(buffer);
			return NULL;
		}
	}


	

	status = cudaMemcpy(results,deviceResults,resultsLength*sizeof(double),cudaMemcpyDeviceToHost);
	if(status != cudaSuccess){
		sprintf(buffer, "results copy back to host failed, errorcode %s", cudaGetErrorString(status));
		this->esl->logMessage(buffer);
		return NULL;
	}

	imageData = initImageData();
	counter = 0;
	for (i=0; i<calcData.ncub2; i++) {
		for (j=0; j<calcData.ncub1; j++) {
			for (k=0; k<calcData.ncub0; k++) {
				imageData->SetScalarComponentFromDouble( k, j, i, 0, results[counter] );
				counter++;
			}
		}
	}

	return imageData;
}

CalcDensCalculateDensityUnrolled::CalcDensCalculateDensityUnrolled(CalcDensDataPack *data): CalcDensCalculateDensity(data){
	
}


dim3 CalcDensCalculateDensityUnrolled::getBlockSize(){
	dim3 result(400,1,1);
	return result;
}

dim3 CalcDensCalculateDensityUnrolled::getGridSize(dim3 blockSize){
	int total, xsize;
	dim3 result(1,1,1);

	total = this->mol->nBasisFunctions * this->calcData.ncub0;
	xsize = total/blockSize.x;
	if(total > xsize * blockSize.x){
		xsize++;
	}
	
	result.x = xsize;
	result.y = this->calcData.ncub1;
	result.z = this->calcData.ncub2;
	return result;
}


