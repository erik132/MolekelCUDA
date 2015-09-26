#include "CalcDensCalcPoint.cuh"

#include "molekelHelpFunctions/CalcChiCalcPoint.cu"
#include "molekelHelpFunctions/CalcChi.cu"

__global__ void calcPoint(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, double *results, double *chiArray){

	double result = 1;
	int indexZ = threadIdx.z + (blockDim.z*blockIdx.z);
	int	indexY = threadIdx.y + (blockDim.y*blockIdx.y);
	int	indexX = threadIdx.x + (blockDim.x*blockIdx.x);
	float x,y,z;
	const int basisFunctions = molecule->nBasisFunctions;
	int i;
	


	if(indexX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
		
		x = internalData.dim0 + indexX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim4 + indexZ*internalData.dz;

		//this needs to be redone. You can not put 1MB or more on a 64KB chip!
		//calcChi(memLocation, molecule, x, y, z);

		/*for(i=0; i<basisFunctions; i++){
			result += memLocation[i] + orbital->deviceCoefficients[i];
		}*/

		results[indexX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] = calcChiCalcPoint(orbital, molecule, x, y, z);
		//results[indexX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] = 123;
	}
	
	if(indexX ==0 && indexY == 0 && indexZ == 0){
		for(i=0; i<molecule->nBasisFunctions; i++)
			chiArray[i] = 0;
		calcChi(chiArray, molecule, x, y, z);
	}
	
}


vtkImageData* CalcDensCalcPoint::calcImageData(){
	
	//CudaMolecule molecule;
	ESLogger esl("CudaCalcDensCalcPoint.txt");
	esl.logMessage("function started");
	cudaError_t status;
	double *results, *deviceResults, *hostChi, *deviceChi;
	const int resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	int i, j, k, counter;
	char buffer[100];
	float x,y,z;
	
	sprintf(buffer, "resultsLength is %d" , resultsLength);
	esl.logMessage(buffer);
	

	results = new double[resultsLength];

	for(i=0; i<resultsLength; i++){
		results[i] = 2;
	}

	dim3 blockSize(BLOCK_DIM,BLOCK_DIM,BLOCK_DIM);
	dim3 gridSize = getGridSize();
	//dim3 gridSize(1);

	status = CalcDensCalcPoint::moleculeToDevice();

	
	
	if(status == cudaSuccess){
		esl.logMessage("molecule copied with success");
	}else{
		esl.logMessage("Molecule copy failed");
	}

	sprintf(buffer, "ncub0 %d, ncub1 %d, ncub2 %d", calcData.ncub0, calcData.ncub1, calcData.ncub2);
	esl.logMessage(buffer);

	status=cudaMalloc((void**)&deviceResults, sizeof(double)*resultsLength);
	status=cudaMalloc((void**)&deviceChi, sizeof(double)*cudaMolecule.nBasisFunctions);
	hostChi = new double[cudaMolecule.nBasisFunctions];

	if(status == cudaSuccess){
		esl.logMessage("memory allocation on device success");
	}else{
		sprintf(buffer, "memory allocation on device failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
	}

	calcPoint<<<gridSize, blockSize>>>(deviceMolecule, calcData, deviceOrbital, deviceResults, deviceChi);
	status = cudaGetLastError();

	if(status == cudaSuccess){
		esl.logMessage("kernel launch success");
	}else{
		sprintf(buffer, "Kernel launch failed, errorcode %s", cudaGetErrorString(status));
		esl.logMessage(buffer);
	}
	status = cudaDeviceSynchronize();

	status = cudaMemcpy(results, deviceResults, sizeof(double)*resultsLength, cudaMemcpyDeviceToHost);
	status = cudaMemcpy(hostChi, deviceChi, sizeof(double)*cudaMolecule.nBasisFunctions, cudaMemcpyDeviceToHost);

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
	/*counter = 0;
	for (i=0, z=calcData.dim4; i<calcData.ncub2; i++, z += calcData.dz) {
		for (j=0, y=calcData.dim2; j<calcData.ncub1; j++, y += calcData.dy) {
			for (k=0, x=calcData.dim0; k<calcData.ncub0; k++, x += calcData.dx) {
				
				if(i < 5 && j <5 && k<5){
					sprintf(buffer, "x: %d y: %d z: %d value is %.15f", k, j, i, results[counter]);
					esl.logMessage(buffer);
				}
				counter++;

			}
		}
	}*/

	/*for(i=0; i<cudaMolecule.nBasisFunctions; i++){
		sprintf(buffer, "nr: %d value is %.15f",  i, hostChi[i]);
		esl.logMessage(buffer);
	}*/

	for(i=0; i<resultsLength; i++){
		sprintf(buffer, "nr %d value is %.15f", i, results[i]);
		esl.logMessage(buffer);
	}


	CalcDensCalcPoint::deleteDeviceData();

	delete[] results;
	delete[] hostChi;
	cudaFree(deviceChi);
	
	return NULL;
}

CalcDensCalcPoint::CalcDensCalcPoint(CalcDensDataPack *data): CalcDensCudaFunction(data){
	//CalcDensCudaFunction::CalcDensCudaFunction(Molecule *mol, float *dim, int *ncubes);
}


