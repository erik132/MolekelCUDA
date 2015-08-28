#include "CalcDensCudaFunction.cuh"


vtkImageData* CalcDensCudaFunction::initImageData(){

	vtkSmartPointer< vtkImageData > image( vtkImageData::New() );

	image->SetDimensions( ncub[0], ncub[1], ncub[2] );
	image->SetOrigin( calcData->dim[0],calcData->dim[2],calcData->dim[4] );
	image->SetSpacing( dx, dy, dz );

	return image;
}

CalcDensCudaFunction::CalcDensCudaFunction(CalcDensDataPack *data){
	
	ESLogger esl("CalcDensCudaFunction.txt");
	esl.logMessage("function started");
	

	ncub[0] = *data->ncubes++;
	ncub[1] = *data->ncubes++;
	ncub[2] = *data->ncubes++;

	dx = (data->dim[1]-data->dim[0])/(ncub[0]-1);
	dy = (data->dim[3]-data->dim[2])/(ncub[1]-1);
	dz = (data->dim[5]-data->dim[4])/(ncub[2]-1);
	
	cudaMolecule.setProperties(data->mol);

	calcData = data;

	

	if(CalcDensCudaFunction::calcData->mol == NULL){
		esl.logMessage("data pack molecule is NULL");
	}else{
		esl.logMessage("data pac molecule is not NULL");
	}
}

CalcDensDataPack *CalcDensCudaFunction::getDataPack(){
	return calcData;
}



cudaError_t CalcDensCudaFunction::moleculeToDevice(){
	
	cudaError_t status;

	status = cudaMolecule.cpyInternalPointers();

	if(status!=cudaSuccess){
		return status;
	}

	status = cudaMalloc(&deviceMolecule, sizeof(CudaMolecule));

	if(status!=cudaSuccess){
		return status;
	}
	status = cudaMemcpy(deviceMolecule, &cudaMolecule, sizeof(CudaMolecule),cudaMemcpyHostToDevice);

	return status;
}

void CalcDensCudaFunction::deleteDeviceData(){

	cudaMolecule.clearCudaData();
	cudaFree(deviceMolecule);

}