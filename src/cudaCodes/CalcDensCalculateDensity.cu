#include "CalcDensCalculateDensity.cuh"

__global__ void checkDensityMatrix(float* densities, float* resultArray, int size){
	int idx = threadIdx.x;
	int i;
	float result = 0;
	densities = densities + (idx*(idx+1))/2;

	for(i=0; i<=idx; i++){
		result += densities[i];
	}
	resultArray[idx] = result;
}

vtkImageData* CalcDensCalculateDensity::calcImageData(){
	ESLogger esl("CalcDensCalculateDensity.txt");
	esl.logMessage("function started");

	const int resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	cudaError_t status;
	char buffer[1000];
	int i, j;
	int CHECK_SIZE = 200;
	float checkSum = 0;
	float *deviceCheckArray;
	if(CHECK_SIZE > mol->nBasisFunctions) CHECK_SIZE = mol->nBasisFunctions;
	float *checkArray;
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
	checkDensityMatrix<<<1, CHECK_SIZE>>>(deviceDensityMatrix,deviceCheckArray, sizeof(float));
	cudaThreadSynchronize();
	status = cudaGetLastError();

	if(status == cudaSuccess){
		esl.logMessage("density matrix check ran");
	}else{
		esl.logMessage("Density matrix check failed");
	}

	cudaMemcpy(checkArray, deviceCheckArray, CHECK_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
	
	for(i=0; i<CHECK_SIZE; i++){
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
	}

	/*for(i=0; i<CHECK_SIZE; i++){
		sprintf(buffer, "check sum at %d value: %f", i, densityMatrix[0][i]);
		esl.logMessage(buffer);
	}*/

	this->deleteDensityMatrix();
	deleteDeviceMoleculeData();
	deleteDeviceDensityMatrix();
	cudaFree(deviceCheckArray);
	free(checkArray);
	return NULL;
}

CalcDensCalculateDensity::CalcDensCalculateDensity(CalcDensDataPack *data): CalcDensCudaFunction(data){

}