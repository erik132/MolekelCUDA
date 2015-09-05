#include "CalcDensCalcPoint.cuh"

#include "molekelHelpFunctions/CalcChi.cu"

__global__ void calcPoint(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, double *results){

	double result = 123456;
	int indexZ = threadIdx.z + blockDim.z*BLOCK_DIM;
	int	indexY = threadIdx.y + blockDim.y*BLOCK_DIM;
	int	indexX = threadIdx.x + blockDim.x*BLOCK_DIM;
	float x,y,z;
	const int basisFunctions = molecule->nBasisFunctions;
	int i;
	


	if(indexX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
		

		x = internalData.dim0 + indexX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim3 + indexZ*internalData.dz;

		for(i = 0; i<molecule->nBasisFunctions; i++){
			memLocation[i] = 0;
		}
		//this needs to be redone. You can not put 1MB or more on a 64KB chip!
		//calcChi(memLocation, molecule, x, y, z);

		for(i=0; i<basisFunctions; i++){
			//result += memLocation[i] + orbital->deviceCoefficients[i];
		}

		
	}
	results[threadIdx.x] = result;
	
}


vtkImageData* CalcDensCalcPoint::calcImageData(){
	
	//CudaMolecule molecule;
	ESLogger esl("CudaCalcDensCalcPoint.txt");
	esl.logMessage("function started");
	cudaError_t status;
	double *results, *deviceResults;
	const int resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	int i, j, k;
	char buffer[100];
	
	sprintf(buffer, "resultsLength is %d" , resultsLength);
	esl.logMessage(buffer);
	/*if(calcData.mol == NULL){
		esl.logMessage("mol pointer was null");
	}else{
		esl.logMessage("mol pointer was not null");
	}*/

	results = new double[resultsLength];

	for(i=0; i<resultsLength; i++){
		results[i] = 2;
	}

	dim3 blockSize(BLOCK_DIM,BLOCK_DIM,BLOCK_DIM);
	dim3 gridSize(1);

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

	calcPoint<<<gridSize, blockSize>>>(deviceMolecule, calcData, deviceOrbital, deviceResults);
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
	
	for(i=0; i<BLOCK_DIM; i++){
		for(j=0; j<BLOCK_DIM; j++){
			for(k=0; k<BLOCK_DIM; k++){
				sprintf(buffer, "x: %d y: %d z: %d value is %.15f", k, j, i, results[k + (BLOCK_DIM*j) + (BLOCK_DIM*BLOCK_DIM*i)]);
				esl.logMessage(buffer);
			}
		}
	}


	CalcDensCalcPoint::deleteDeviceData();

	delete[] results;
	
	return NULL;
}

CalcDensCalcPoint::CalcDensCalcPoint(CalcDensDataPack *data): CalcDensCudaFunction(data){
	//CalcDensCudaFunction::CalcDensCudaFunction(Molecule *mol, float *dim, int *ncubes);
}


