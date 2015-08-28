#include "CalcDensCalcPoint.cuh"

#include "molekelHelpFunctions/CalcChi.cu"

__global__ void calcPoint(CudaMolecule *molecule, float dx, float dy, float dz, CalcDensInternalData specialData){


	extern __shared__ double chi[];
	int indexZ = threadIdx.z + blockDim.z*BLOCK_DIM;
	int	indexY = threadIdx.y + blockDim.y*BLOCK_DIM;
	int	indexX = threadIdx.x + blockDim.x*BLOCK_DIM;
	float x,y,z;


	if(indexX < specialData.ncub0 && indexY < specialData.ncub1 && indexZ < specialData.ncub2){
		x = specialData.dim0 + indexX*specialData.dx;
		y = specialData.dim2 + indexY*specialData.dy;
		z = specialData.dim3 + indexZ*specialData.dz;

		calcChi(chi, molecule, x, y, z);
	}
	
	
	
}


vtkImageData* CalcDensCalcPoint::calcImageData(){
	
	//CudaMolecule molecule;
	ESLogger esl("CudaCalcDensCalcPoint.txt");
	esl.logMessage("function started");
	cudaError_t status;
	
	if(calcData->mol == NULL){
		esl.logMessage("mol pointer was null");
	}else{
		esl.logMessage("mol pointer was not null");
	}

	dim3 blockSize(5,5,5);
	int memSize = 5*5*5*sizeof(double)*cudaMolecule.nBasisFunctions;
	dim3 gridSize(1);

	
	status = CalcDensCalcPoint::moleculeToDevice();

	calcPoint<<<gridSize, blockSize, memSize>>>(deviceMolecule, dx, dy, dz);
	
	if(status == cudaSuccess){
		esl.logMessage("molecule copied with success");
	}else{
		esl.logMessage("Molecule copy failed");
	}

	CalcDensCalcPoint::deleteDeviceData();
	
	return NULL;
}

CalcDensCalcPoint::CalcDensCalcPoint(CalcDensDataPack *data): CalcDensCudaFunction(data){
	//CalcDensCudaFunction::CalcDensCudaFunction(Molecule *mol, float *dim, int *ncubes);
}


