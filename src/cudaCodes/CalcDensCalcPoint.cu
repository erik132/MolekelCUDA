#include "CalcDensCalcPoint.cuh"

#include "molekelHelpFunctions/CalcChiCalcPoint.cu"
#include "gputimer.h"

__global__ void calcPoint(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, double *results){

	double result = 0;
	int indexZ = threadIdx.z + (blockDim.z*blockIdx.z);
	int	indexY = threadIdx.y + (blockDim.y*blockIdx.y);
	int	indexX = threadIdx.x + (blockDim.x*blockIdx.x);
	float x,y,z;

	if(indexX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
		
		x = internalData.dim0 + indexX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim4 + indexZ*internalData.dz;
		
		result = calcChiCalcPoint(orbital, molecule, x, y, z);
		results[indexX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] = result;
	}
	
}

cudaError_t CalcDensCalcPoint::initData(){
	ESLogger esl("CudaCalcDensCalcPoint.txt");
	char buffer[100];
	cudaError_t status;

	resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	results = new double[resultsLength];
	
	status = CalcDensCalcPoint::moleculeToDevice();
	if(status != cudaSuccess){
		esl.logMessage("Molecule copy failed");
		return status;
	}

	status = CalcDensCalcPoint::orbitalToDevice();
	if(status != cudaSuccess){
		esl.logMessage("Orbital copy failed");
		return status;
	}

	status=cudaMalloc((void**)&deviceResults, sizeof(double)*resultsLength);
	if(status != cudaSuccess){
		sprintf(buffer, "memory allocation on device failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
		return status;
	}

	return cudaSuccess;
}

vtkImageData* CalcDensCalcPoint::runComputation(){
	
	ESLogger esl("CudaCalcDensCalcPoint.txt");
	cudaError_t status;
	int i, j, k, counter;
	char buffer[100];
	vtkImageData* imageData;
	
	dim3 blockSize(BLOCK_DIM,BLOCK_DIM,BLOCK_DIM);
	dim3 gridSize = getGridSize();

	
	calcPoint<<<gridSize, blockSize>>>(deviceMolecule, calcData, deviceOrbital, deviceResults);

	status = cudaGetLastError();
	if(status != cudaSuccess){
		sprintf(buffer, "Kernel launch failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
	}
	status = cudaDeviceSynchronize();
	if(status != cudaSuccess){
		sprintf(buffer, "Device synchronization failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
	}

	status = cudaMemcpy(results, deviceResults, sizeof(double)*resultsLength, cudaMemcpyDeviceToHost);
	if(status != cudaSuccess){
		sprintf(buffer, "memcpy from device failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
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

void CalcDensCalcPoint::cleanupData(){
	CalcDensCalcPoint::deleteDeviceMoleculeData();
	CalcDensCalcPoint::deleteDeviceOrbitalData();
	cudaFree(deviceResults);
	delete[] results;
}

CalcDensCalcPoint::CalcDensCalcPoint(CalcDensDataPack *data): CalcDensCudaFunction(data){
	
}


