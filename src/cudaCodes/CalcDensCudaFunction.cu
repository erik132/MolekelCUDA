#include "CalcDensCudaFunction.cuh"


vtkImageData* CalcDensCudaFunction::initImageData(){

	vtkSmartPointer< vtkImageData > image( vtkImageData::New() );

	image->SetDimensions( calcData.ncub0, calcData.ncub1, calcData.ncub2 );
	image->SetOrigin( calcData.dim0,calcData.dim2,calcData.dim4 );
	image->SetSpacing( calcData.dx, calcData.dy, calcData.dz );

	return image;
}

CalcDensCudaFunction::CalcDensCudaFunction(CalcDensDataPack *data){
	
	/*ESLogger esl("CalcDensCudaFunction.txt");
	esl.logMessage("function started");*/
	BLOCK_DIM = 5;

	densityMatrix = NULL;
	calcData.minValue = data->minValue;
	calcData.maxValue = data->maxValue;

	calcData.ncub0 = *data->ncubes++;
	calcData.ncub1 = *data->ncubes++;
	calcData.ncub2 = *data->ncubes++;

	calcData.dim0 = data->dim[0];
	calcData.dim2 = data->dim[2];
	calcData.dim4 = data->dim[4];
	//esl.logMessage("starting to cpy datasource");
	calcData.datasource = data->datasource;

	calcData.dx = (data->dim[1] - data->dim[0]) / (calcData.ncub0 - 1);
	calcData.dy = (data->dim[3] - data->dim[2]) / (calcData.ncub1 - 1);
	calcData.dz = (data->dim[5] - data->dim[4]) / (calcData.ncub2 - 1);
	//esl.logMessage("starting to cpy mol and key");
	this->mol = data->mol;
	this->key = data->key;
	//esl.logMessage("mol and key cpy success");

	if(data->orbital != NULL){
		cudaOrbital.setProperties(data->orbital,data->mol);
	}
	cudaMolecule.setProperties(data->mol);
	calcData.densityLength = 0;

	
}



cudaError_t CalcDensCudaFunction::moleculeToDevice(){
	
	cudaError_t status;

	status = cudaMolecule.cpyInternalPointers();

	if(status!=cudaSuccess){
		return status;
	}

	status = cudaMalloc((void**)&deviceMolecule, sizeof(CudaMolecule));

	if(status!=cudaSuccess){
		return status;
	}
	status = cudaMemcpy(deviceMolecule, &cudaMolecule, sizeof(CudaMolecule),cudaMemcpyHostToDevice);

	return status;
}

void CalcDensCudaFunction::deleteDeviceMoleculeData(){

	cudaMolecule.clearCudaData();
	cudaFree(deviceMolecule);
	

}

void CalcDensCudaFunction::deleteDeviceOrbitalData(){
	cudaOrbital.clearCudaData();
	cudaFree(deviceOrbital);
}

cudaError_t CalcDensCudaFunction::orbitalToDevice(){

	cudaError_t status = cudaOrbital.cpyInternalPointers();
	if(status != cudaSuccess){
		return status;
	}

	status = cudaMalloc((void**)&deviceOrbital, sizeof(CudaMolecularOrbital));

	if(status != cudaSuccess){
		return status;
	}
	status = cudaMemcpy(deviceOrbital, &cudaOrbital, sizeof(CudaMolecularOrbital),cudaMemcpyHostToDevice);

	return status;
}

dim3 CalcDensCudaFunction::getGridSize(){
	
	int x,y,z; //size nr for each direction
	
	x= getSingleGridSize(calcData.ncub0, BLOCK_DIM);
	y= getSingleGridSize(calcData.ncub1, BLOCK_DIM);
	z= getSingleGridSize(calcData.ncub2, BLOCK_DIM);

	dim3 gridSize(x,y,z);
	return gridSize;

}

int CalcDensCudaFunction::getSingleGridSize(int elements, int blockSize){
	int result;

	result = elements/blockSize;
	if(result*blockSize < elements){
		result++;
	}

	return result;
}

int CalcDensCudaFunction::createDensityMatrix(){
	register short i, j, k;
	float adder;

	if(!mol->alphaDensity && calcData.datasource == USE_MATRICES) return 0;

	if((densityMatrix = (float **)alloc_trimat(mol->nBasisFunctions, sizeof(float))) == NULL){
		fprintf(stderr, "can't allocate density-matrix\n");
		return 0;
	}

	/* use the alpha (and beta) density matrices from the gaussian output file */
	if(calcData.datasource == USE_MATRICES){
	if(key == EL_DENS){
	  if(mol->betaDensity){
		for(i=0; i<mol->nBasisFunctions; i++){
		  for(j=0; j<=i; j++)
			densityMatrix[i][j] = mol->alphaDensity[i][j] + mol->betaDensity[i][j];
		}
	  }
	  else {
		for(i=0; i<mol->nBasisFunctions; i++){
		  for(j=0; j<=i; j++) densityMatrix[i][j] = mol->alphaDensity[i][j];
		}
	  }
	}
	if(key == SPIN_DENS){
	  if(!mol->betaDensity) return 0;
	  for(i=0; i<mol->nBasisFunctions; i++){
		for(j=0; j<=i; j++)
		  densityMatrix[i][j]  = mol->alphaDensity[i][j] - mol->betaDensity[i][j];
	  }
	}
	}

	/* generate the density matrices with the coefficients */
	else if(calcData.datasource == USE_COEFFS) {
	if(key == EL_DENS){
	  for(i=0; i<mol->nBasisFunctions; i++){
		for(j=0; j<=i; j++){
		  densityMatrix[i][j] = 0;
		  for(k=0; k<mol->nBeta; k++){
			adder = mol->alphaOrbital[k].coefficient[i] *
				  mol->alphaOrbital[k].coefficient[j];
			if(mol->alphaBeta)
			  adder += mol->betaOrbital[k].coefficient[i] *
					mol->betaOrbital[k].coefficient[j];
			else adder *= 2.0;
			densityMatrix[i][j] += adder;
		  }
		  for(; k<mol->nAlpha; k++)
			densityMatrix[i][j] += mol->alphaOrbital[k].coefficient[i] *
						mol->alphaOrbital[k].coefficient[j];
		}
	  }
	}
	if(key == SPIN_DENS){
	  if(!mol->alphaBeta) {
		for(i=0; i<mol->nBasisFunctions; i++){
		  for(j=0; j<=i; j++){
			densityMatrix[i][j] = 0;
			for(k=mol->nBeta; k<mol->nAlpha; k++)
			  densityMatrix[i][j] += mol->alphaOrbital[k].coefficient[i] *
						  mol->alphaOrbital[k].coefficient[j];
		  }
		}
	  }
	  else {
		for(i=0; i<mol->nBasisFunctions; i++){
		  for(j=0; j<=i; j++){
			densityMatrix[i][j] = 0;
			for(k=0; k<mol->nAlpha; k++)
			  densityMatrix[i][j] += mol->alphaOrbital[k].coefficient[i] *
						  mol->alphaOrbital[k].coefficient[j];
			for(k=0; k<mol->nBeta; k++)
			  densityMatrix[i][j] -= mol->betaOrbital[k].coefficient[i] *
						  mol->betaOrbital[k].coefficient[j];
		  }
		}
	  }
	}
	}
	else return 0;

	return 1;
}

void CalcDensCudaFunction::deleteDeviceDensityMatrix(){
	cudaFree(deviceDensityMatrix);
}

/* take straight from old/utilites.cpp */
void *CalcDensCudaFunction::alloc_trimat(int n, size_t size)
{
	void **pointerarray;
	char *array;
	register short i;

	if((array = (char*) malloc((n*(n+1))/2*size)) == NULL) return NULL;
	  /* array will hold the data */
	if((pointerarray = (void**) malloc(n*sizeof(char *))) == NULL) return NULL;
	  /* pointerarray will hold the pointers to the rows of data */
	for(i=0; i<n; i++) pointerarray[i] = array + (i*(i+1))/2*size;

	return pointerarray;
}

void CalcDensCudaFunction::deleteDensityMatrix(){
	if(densityMatrix){
		free(densityMatrix[0]);
		free(densityMatrix);
		densityMatrix = NULL;
	}
}

cudaError_t CalcDensCudaFunction::densityMatrixToDevice(){
	
	densityMatrixLength = (mol->nBasisFunctions*(mol->nBasisFunctions+1))/2*sizeof(float);
	cudaError_t status;

	status = cudaMalloc((void**)&deviceDensityMatrix, densityMatrixLength);
	if(status != cudaSuccess) return status;
	status = cudaMemcpy(deviceDensityMatrix, densityMatrix[0], densityMatrixLength, cudaMemcpyHostToDevice);
	calcData.densityLength = densityMatrixLength/sizeof(float);
	return status;

}

vtkImageData* CalcDensCudaFunction::calcImageData(){
	vtkImageData* result = NULL;

	if(initData() == cudaSuccess){
		result = runComputation();
	}
	cleanupData();
	return result;
}