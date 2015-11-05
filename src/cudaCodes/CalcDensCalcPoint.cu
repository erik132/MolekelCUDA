#include "CalcDensCalcPoint.cuh"

#include "molekelHelpFunctions/CalcChiCalcPoint.cu"
#include "molekelHelpFunctions/CalcChi.cu"

#include "gputimer.h"

__global__ void calcPoint(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, double *results){

	double result = 1;
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


vtkImageData* CalcDensCalcPoint::calcImageData(){
	
	ESLogger esl("CudaCalcDensCalcPoint.txt");
	esl.logMessage("function started");
	cudaError_t status;
	double *results, *deviceResults;
	const int resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	int i, j, k, counter;
	char buffer[100];
	vtkImageData* imageData;
	GpuTimer gputimer;
	
	sprintf(buffer, "resultsLength is %d" , resultsLength);
	esl.logMessage(buffer);
	

	results = new double[resultsLength];

	dim3 blockSize(BLOCK_DIM,BLOCK_DIM,BLOCK_DIM);
	dim3 gridSize = getGridSize();

	status = CalcDensCalcPoint::moleculeToDevice();

	
	
	if(status == cudaSuccess){
		esl.logMessage("molecule copied with success");
	}else{
		esl.logMessage("Molecule copy failed");
	}

	sprintf(buffer, "ncub0 %d, ncub1 %d, ncub2 %d", calcData.ncub0, calcData.ncub1, calcData.ncub2);
	esl.logMessage(buffer);

	status=cudaMalloc((void**)&deviceResults, sizeof(double)*resultsLength);

	if(status == cudaSuccess){
		esl.logMessage("memory allocation on device success");
	}else{
		sprintf(buffer, "memory allocation on device failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
	}
	
	gputimer.Start();
	calcPoint<<<gridSize, blockSize>>>(deviceMolecule, calcData, deviceOrbital, deviceResults);
	gputimer.Stop();
	sprintf(buffer, "kernel run time: %f", gputimer.Elapsed());
	esl.logMessage(buffer);
	status = cudaGetLastError();

	if(status == cudaSuccess){
		esl.logMessage("kernel launch success");
	}else{
		sprintf(buffer, "Kernel launch failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
	}
	status = cudaDeviceSynchronize();

	status = cudaMemcpy(results, deviceResults, sizeof(double)*resultsLength, cudaMemcpyDeviceToHost);

	if(status == cudaSuccess){
		esl.logMessage("memcpy from device success");
	}else{
		esl.logMessage("memcpy from device failed");
	}
	cudaFree(deviceResults);
	
	/*for(i=0; i<BLOCK_DIM; i++){
		for(j=0; j<BLOCK_DIM; j++){
			for(k=0; k<BLOCK_DIM; k++){
				sprintf(buffer, "x: %d y: %d z: %d value is %.15f", k, j, i, results[k + (BLOCK_DIM*j) + (BLOCK_DIM*BLOCK_DIM*i)]);
				esl.logMessage(buffer);
			}
		}
	}*/
	imageData = initImageData();
	counter = 0;
	for (i=0; i<calcData.ncub2; i++) {
		for (j=0; j<calcData.ncub1; j++) {
			for (k=0; k<calcData.ncub0; k++) {
				
				/*if(i < 5 && j <5 && k<5){
					sprintf(buffer, "x: %d y: %d z: %d value is %.15f", k, j, i, results[counter]);
					esl.logMessage(buffer);
				}*/
				//imageData->SetScalarComponentFromDouble( k, j, i, 0, results[k + (calcData.ncub0*j) + (calcData.ncub0*calcData.ncub1*i)] );
				imageData->SetScalarComponentFromDouble( k, j, i, 0, results[counter] );
				counter++;

			}
		}
	}

	/*for(i=0; i<cudaMolecule.nBasisFunctions; i++){
		sprintf(buffer, "nr: %d value is %.15f",  i, hostChi[i]);
		esl.logMessage(buffer);
	}*/

	/*for(i=0; i<resultsLength; i++){
		sprintf(buffer, "nr %d value is %.15f", i, results[i]);
		esl.logMessage(buffer);
	}*/


	CalcDensCalcPoint::deleteDeviceData();

	delete[] results;
	
	return imageData;
}

CalcDensCalcPoint::CalcDensCalcPoint(CalcDensDataPack *data): CalcDensCudaFunction(data){
	
}


