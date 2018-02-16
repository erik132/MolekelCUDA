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



cudaError_t CalcDensCalculateDensity::initData(){
	char buffer[100];
	cudaError_t status;

	resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	results = new double[resultsLength];

	this->esl->logMessage("starting to init data");
	
	status = this->moleculeToDevice();
	if(status != cudaSuccess){
		esl->logMessage("Molecule copy failed");
		return status;
	}

	status = this->orbitalToDevice();
	if(status != cudaSuccess){
		esl->logMessage("Orbital copy failed");
		return status;
	}

	status=cudaMalloc((void**)&deviceResults, sizeof(double)*resultsLength);
	if(status != cudaSuccess){
		sprintf(buffer, "memory allocation on device failed, errorcode %s", cudaGetErrorString(status));
		esl->logMessage(buffer);
		return status;
	}
	
	if(!this->createDensityMatrix()){
		this->esl->logMessage("can not create density matrix.");
		return cudaErrorMemoryAllocation;
	}

	status = this->densityMatrixToDevice();
	if(status != cudaSuccess){
		sprintf(buffer, "density matrix allocation to device failed, errorcode %s", cudaGetErrorString(status));
		this->esl->logMessage(buffer);
	}

	this->esl->logMessage("all operations on initdata finished as a success");
	return cudaSuccess;
}

vtkImageData* CalcDensCalculateDensity::runComputation(){

	return NULL;
}

void CalcDensCalculateDensity::cleanupData(){
	this->esl->logMessage("starting to clean up.");

	this->deleteDeviceMoleculeData();
	this->deleteDeviceOrbitalData();
	this->deleteDeviceDensityMatrix();
	this->deleteDensityMatrix();
	cudaFree(deviceResults);
	delete[] results;

	this->esl->logMessage("cleanup complete.");
}

CalcDensCalculateDensity::CalcDensCalculateDensity(CalcDensDataPack *data): CalcDensCudaFunction(data){
	this->esl = new ESLogger("CalcDensCalculateDensity.txt");
}

CalcDensCalculateDensity::~CalcDensCalculateDensity(){
	delete this->esl;
}