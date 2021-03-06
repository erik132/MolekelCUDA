#include "CalcDensCalculateDensity.cuh"

#include "molekelHelpFunctions/CalcChiCalculateDensity.cu"

__global__ void calculateDensity(CudaMolecule *molecule, CalcDensInternalData internalData, CudaMolecularOrbital *orbital, float *densities, double *results, int offsetx){
	double result = 0;
	int indexZ = threadIdx.z + (blockDim.z*blockIdx.z);
	int	indexY = threadIdx.y + (blockDim.y*blockIdx.y);
	int	indexX = threadIdx.x + (blockDim.x*(blockIdx.x + offsetx));
	float x,y,z;

	if(indexX < internalData.ncub0 && indexY < internalData.ncub1 && indexZ < internalData.ncub2){
		
		x = internalData.dim0 + indexX*internalData.dx;
		y = internalData.dim2 + indexY*internalData.dy;
		z = internalData.dim4 + indexZ*internalData.dz;
		
		result = calcChiCalculateDensity(densities,orbital, molecule, x, y, z, internalData.densityLength);
		results[indexX + (internalData.ncub0*indexY) + (internalData.ncub0*internalData.ncub1*indexZ)] = result;
	}
}



cudaError_t CalcDensCalculateDensity::initData(){
	char buffer[100];
	cudaError_t status;

	resultsLength = calcData.ncub0*calcData.ncub1*calcData.ncub2;
	results = new double[resultsLength];

	
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

	status = cudaMemset(deviceResults,0,sizeof(double)*resultsLength);
	if(status != cudaSuccess){
		sprintf(buffer, "can not memset device results to 0, errorcode %s", cudaGetErrorString(status));
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

	return cudaSuccess;
}

vtkImageData* CalcDensCalculateDensity::runComputation(){
	char buffer[100];
	cudaError_t status;
	int i, j, k, counter=0;
	int originalGridx;
	int offsetx=0;
	vtkImageData* imageData;


	dim3 blockSize(BLOCK_DIM,BLOCK_DIM,BLOCK_DIM);
	dim3 gridSize = getGridSize();
	
	originalGridx = gridSize.x;
	gridSize = this->limitThreads(10000,gridSize);

	
	for(offsetx=0; offsetx<originalGridx; offsetx += gridSize.x){
		calculateDensity<<<gridSize, blockSize>>>(this->deviceMolecule,this->calcData,this->deviceOrbital,this->deviceDensityMatrix, this->deviceResults,offsetx);
		status = cudaDeviceSynchronize();
		if(status != cudaSuccess){
			sprintf(buffer, "failed to sync devices during loop %s", cudaGetErrorString(status));
			this->esl->logMessage(buffer);
			return NULL;
		}
	}

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

void CalcDensCalculateDensity::cleanupData(){

	this->deleteDeviceMoleculeData();
	this->deleteDeviceOrbitalData();
	this->deleteDeviceDensityMatrix();
	this->deleteDensityMatrix();
	cudaFree(deviceResults);
	delete[] results;

}

CalcDensCalculateDensity::CalcDensCalculateDensity(CalcDensDataPack *data): CalcDensCudaFunction(data){
	this->esl = new ESLogger("CalcDensCalculateDensity.txt");
}

CalcDensCalculateDensity::~CalcDensCalculateDensity(){
	delete this->esl;
}