#include "CalcDensCudaFunction.cuh"


vtkImageData* CalcDensCudaFunction::initImageData(){

	vtkSmartPointer< vtkImageData > image( vtkImageData::New() );

	image->SetDimensions( calcData.ncub0, calcData.ncub1, calcData.ncub2 );
	image->SetOrigin( calcData.dim0,calcData.dim2,calcData.dim4 );
	image->SetSpacing( calcData.dx, calcData.dy, calcData.dz );

	return image;
}

CalcDensCudaFunction::CalcDensCudaFunction(CalcDensDataPack *data){
	
	ESLogger esl("CalcDensCudaFunction.txt");
	esl.logMessage("function started");
	

	calcData.ncub0 = *data->ncubes++;
	calcData.ncub1 = *data->ncubes++;
	calcData.ncub2 = *data->ncubes++;

	calcData.dim0 = data->dim[0];
	calcData.dim2 = data->dim[2];
	calcData.dim4 = data->dim[4];

	calcData.dx = (data->dim[1] - data->dim[0]) / (calcData.ncub0 - 1);
	calcData.dy = (data->dim[3] - data->dim[2]) / (calcData.ncub1 - 1);
	calcData.dz = (data->dim[5] - data->dim[4]) / (calcData.ncub2 - 1);
	
	cudaMolecule.setProperties(data->mol);
	cudaOrbital.setProperties(data->orbital,data->mol);

	

	if(data->orbital == NULL){
		esl.logMessage("data pack molecularOrbital is NULL");
	}else{
		esl.logMessage("data pac molecularOrbital is not NULL");
	}
}



cudaError_t CalcDensCudaFunction::moleculeToDevice(){
	
	cudaError_t status;

	status = cpyOrbital();

	if(status != cudaSuccess){
		return status;
	}


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

void CalcDensCudaFunction::deleteDeviceData(){

	cudaMolecule.clearCudaData();
	cudaFree(deviceMolecule);
	cudaOrbital.clearCudaData();
	cudaFree(deviceOrbital);

}

cudaError_t CalcDensCudaFunction::cpyOrbital(){

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