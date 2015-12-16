#include "CalcDensCalculateDensity.cuh"

#include "molekelHelpFunctions/CalcChi.cu"

#include "gputimer.h"

__global__ void checkDensityMatrix(float* densities, float* resultArray){
	int idx = threadIdx.x;
	int i;
	float result = 0;
	densities = densities + (idx*(idx+1))/2;

	for(i=0; i<=idx; i++){
		result += densities[i];
	}
	resultArray[idx] = result;
}

__global__ void calculateDensity(/*float *densityMatrix, double * results, CudaMolecule *molecule, 
								 CalcDensInternalData internalData, dim3 maxes, size_t *bufferPointers, size_t ranblocks*/){

	/*int indexX;
	int indexY;
	int indexZ;
	int totalBlockId = ranblocks + blockIdx.x;
	double *bufferBlock;
	float x,y,z;

	bufferBlock = (double*)bufferPointers[blockIdx.x];

	indexZ = totalBlockId / (maxes.x*maxes.y);
	indexY = (totalBlockId - (indexZ * maxes.x * maxes.y)) /  maxes.x;
	indexX = totalBlockId - (indexZ * maxes.x * maxes.y) - (indexY * maxes.x);

	bufferBlock += (threadIdx.x + (blockDim.x*threadIdx.y) + (blockDim.x*blockDim.y*threadIdx.z)) * molecule->nBasisFunctions;

	indexZ = threadIdx.z + (blockDim.z*indexZ);
	indexY = threadIdx.y + (blockDim.y*indexY);
	indexX = threadIdx.x + (blockDim.x*indexX);

	if(indexX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
		x = internalData.dim0 + indexX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim4 + indexZ*internalData.dz;

		calcChi(bufferBlock, molecule, x, y, z);

		results[indexX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] = bufferBlock[500];
	}*/


}

vtkImageData* CalcDensCalculateDensity::calcImageData(){
	ESLogger esl("CalcDensCalculateDensity.txt");
	esl.logMessage("function started");

	const int resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	cudaError_t status;
	int k, intBufferBlockCount;
	char buffer[1000];
	size_t i,j;
	size_t CHECK_SIZE = 200;
	float checkSum = 0;
	float *deviceCheckArray;
	cudaDeviceProp devProps;
	float *checkArray;
	size_t totalBufferSize, blockBufferSize, bufferBlockCount, *buffers, *deviceBuffers;
	dim3 gridSize = getGridSize(), blockSize(BLOCK_DIM, BLOCK_DIM, BLOCK_DIM);
	size_t totalBlocks = gridSize.x * gridSize.y * gridSize.z;
	double *deviceBufferBlock;
	double *results, *deviceResults;
	
	results = new double[resultsLength];
	
	/* calculate all memory data */
	//if(CHECK_SIZE > mol->nBasisFunctions) CHECK_SIZE = mol->nBasisFunctions;
	cudaGetDeviceProperties(&devProps, 0);
	totalBufferSize = devProps.totalGlobalMem * 0.7;
	blockBufferSize = cudaMolecule.nBasisFunctions * sizeof(double) * BLOCK_DIM * BLOCK_DIM * BLOCK_DIM;
	bufferBlockCount = totalBufferSize / blockBufferSize;

	if(bufferBlockCount > totalBlocks){
		bufferBlockCount = totalBlocks;
	}

	buffers = new size_t[bufferBlockCount];
	

	checkArray = (float*)malloc(CHECK_SIZE*sizeof(float));
	
	sprintf(buffer, "resultsLength is %d" , resultsLength);
	esl.logMessage(buffer);

	status = moleculeToDevice();
	
	if(status == cudaSuccess){
		esl.logMessage("Molecule copied with success");
	}else{
		esl.logMessage("Molecule copy failed");
	}

	if(this->createDensityMatrix()){
		esl.logMessage("density matrix success");
	}else{
		esl.logMessage("density matrix failed");
	}

	status = densityMatrixToDevice();

	if(status == cudaSuccess){
		esl.logMessage("Density matrix copied with success");
	}else{
		esl.logMessage("Density matrix copy failed");
	}
	
	/*for(i=0; i<mol->nBasisFunctions; i++){
		for(j=0; j<=i; j++){
			sprintf(buffer, "matrix point at %d %d value: %f", i, j, densityMatrix[i][j]);
			esl.logMessage(buffer);
		}
	}*/

	cudaMalloc(&deviceCheckArray, CHECK_SIZE*sizeof(float));
	cudaMalloc(&deviceResults, resultsLength);

	for(i=0; i<bufferBlockCount; i++){
		status = cudaMalloc(&deviceBufferBlock, blockBufferSize);
		if(status != cudaSuccess){
			esl.logMessage("Cuda block allocation interrupted");
			bufferBlockCount = i;
			break;
		}
		buffers[i] = (size_t)deviceBufferBlock;
		deviceBufferBlock = NULL;
	}

	if(bufferBlockCount <= 0){
		return NULL;
	}

	cudaMalloc(&deviceBuffers, sizeof(size_t)*bufferBlockCount);
	cudaMemcpy(deviceBuffers, buffers, sizeof(size_t)*bufferBlockCount, cudaMemcpyDeviceToHost);

	intBufferBlockCount = (int)bufferBlockCount;
	dim3 realGridSize(intBufferBlockCount,1,1);
	sprintf(buffer, "Buffer block count value %d", intBufferBlockCount);
	esl.logMessage(buffer);
	for(i=0; i<totalBlocks; i += bufferBlockCount){
		calculateDensity<<<1,1>>>(/*deviceDensityMatrix, deviceResults, deviceMolecule, 
															calcData, gridSize, deviceBuffers, i*/);
		status = cudaGetLastError();
		if(status == cudaSuccess){
			esl.logMessage("calculate density ran");
		}else{
			sprintf(buffer, "calculate density failed. Error is: %s", cudaGetErrorString(status));
			esl.logMessage(buffer);
		}
		cudaThreadSynchronize();

	}

	/*checkDensityMatrix<<<1, CHECK_SIZE>>>(deviceDensityMatrix,deviceCheckArray);
	cudaThreadSynchronize();
	status = cudaGetLastError();

	if(status == cudaSuccess){
		esl.logMessage("density matrix check ran");
	}else{
		esl.logMessage("Density matrix check failed");
	}*/

	cudaMemcpy(checkArray, deviceCheckArray, CHECK_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(results, deviceResults, resultsLength*sizeof(double), cudaMemcpyDeviceToHost);
	/*for(i=0; i<CHECK_SIZE; i++){
		sprintf(buffer, "CUDA check sum at %d value: %f", i, checkArray[i]);
		esl.logMessage(buffer);
	}

	for(i=0; i<CHECK_SIZE; i++){
		checkSum = 0;
		for(j=0; j<=i; j++){
			checkSum += densityMatrix[i][j];
		}
		sprintf(buffer, "check sum at %d value: %f", i, checkSum);
		esl.logMessage(buffer);
	}*/

	for(k=0; k<resultsLength; k++){
		sprintf(buffer, "result check %d value: %.15f", k, results[k]);
		esl.logMessage(buffer);
	}

	/*for(i=0; i<CHECK_SIZE; i++){
		sprintf(buffer, "check sum at %d value: %f", i, densityMatrix[0][i]);
		esl.logMessage(buffer);
	}*/

	
	deleteDeviceMoleculeData();
	deleteDeviceDensityMatrix();
	this->deleteDensityMatrix();
	for(i=0; i<bufferBlockCount; i++){
		deviceBufferBlock = (double*)buffers[i];
		cudaFree(deviceBufferBlock);
	}

	delete[] buffers;
	delete[] results;
	cudaFree(deviceCheckArray);
	cudaFree(deviceBuffers);
	cudaFree(deviceResults);
	free(checkArray);
	return NULL;
}

CalcDensCalculateDensity::CalcDensCalculateDensity(CalcDensDataPack *data): CalcDensCudaFunction(data){

}