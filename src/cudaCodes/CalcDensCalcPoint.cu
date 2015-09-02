#include "CalcDensCalcPoint.cuh"

#include "molekelHelpFunctions/CalcChi.cu"

__global__ void calcPoint(CudaMolecule *molecule, CalcDensInternalData internal, MolecularOrbital *orbital){

	extern __shared__ double chi[];
	int indexZ = threadIdx.z + blockDim.z*BLOCK_DIM;
	int	indexY = threadIdx.y + blockDim.y*BLOCK_DIM;
	int	indexX = threadIdx.x + blockDim.x*BLOCK_DIM;
	float x,y,z;
	double *memLocation;
	const int memNr = threadIdx.x + (threadIdx.y*BLOCK_DIM) + (threadIdx.z*BLOCK_DIM*BLOCK_DIM);
	const int basisFunctions = molecule->nBasisFunctions;
	int i, value=0;


	if(indexX < internal.ncub0 && indexY < internal.ncub1 && indexZ < internal.ncub2){
		memLocation = chi + (memNr*molecule->nBasisFunctions);

		x = internal.dim0 + indexX*internal.dx;
		y = internal.dim2 + indexY*internal.dy;
		z = internal.dim3 + indexZ*internal.dz;

		for(i = 0; i<molecule->nBasisFunctions; i++){
			memLocation[i] = 0;
		}

		calcChi(memLocation, molecule, x, y, z);

		for(i=0; i<basisFunctions; i++){
			value += memLocation[i] + orbital[i];
		}
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

	dim3 blockSize(BLOCK_DIM,BLOCK_DIM,BLOCK_DIM);
	int memSize = BLOCK_DIM*BLOCK_DIM*BLOCK_DIM*sizeof(double)*cudaMolecule.nBasisFunctions;
	dim3 gridSize(1);

	
	status = CalcDensCalcPoint::moleculeToDevice();

	calcPoint<<<gridSize, blockSize, memSize>>>(deviceMolecule, internal, deviceOrbital);
	
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


